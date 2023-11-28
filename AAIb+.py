#!/usr/bin/env python

__version__ = '1.0.1'


import gzip
import os
import shutil
import sys
import subprocess as sp
from argparse import ArgumentParser
from collections import Counter
from glob import glob
from multiprocessing import cpu_count
from tempfile import mkdtemp
from time import strftime

from Bio import SeqIO
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy import mean, std

def parseArgs():
	parser = ArgumentParser(description='Computes the average amino acid '
			'identity (AAI) between two protein sets', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-1', '--set1', required=True, metavar='FILE',
		help='first input FastA or GenBank file '
		'(optionally gunzip compressed)')
	req.add_argument('-2', '--set2', required=True, metavar='FILE',
		help='second input FastA or GenBank file '
		'(optionally gunzip compressed)')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--cpus', type=require_int_nonnegative,
		metavar='INT', default='0', help='number of CPUs [all]')
	opt.add_argument('-f', '--fraction', type=float, metavar='FLOAT',
		default=70.0, help='minimum alignment length percentage [70.0]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-i', '--identity', type=require_float_0_to_100,
		metavar='FLOAT', default=30.0, help='minimum percent identity [30.0]')
	opt.add_argument('-l', '--length', type=int, metavar='INT',
		default=0, help='minimum alignment character length (sum of all '
		'aligned segments and all gaps) [0]')
	opt.add_argument('-o', '--outpath', metavar='PATH',
		default=None, help='output directory [AAI--<date>_<time>]')
	opt.add_argument('-s', '--bitscore', type=float, metavar='FLOAT',
		default=0.0, help='minimum alignment Bit score [0.0]')
	opt.add_argument('--aligner', choices=['blastp', 'diamond'],
		default='blastp', help='local protein aligner [blastp]')
	opt.add_argument('--decimal-places', type=require_int_nonnegative,
		metavar='INT', default=3, help='precision decimal points to round '
		'AAI values [3]')
	opt.add_argument('--max-ACGT', type=require_float_0_to_1,
		metavar='FLOAT', default=0.9, help='input sequence sets above this '
		'fraction of ACGT will be assumed to be nucleotides and AAI '
		'skipped [0.9]')
	opt.add_argument('--min-aln-frac', type=require_float_0_to_1,
		metavar='FLOAT', default=0, help='minimum two-way alignment fraction '
		'of each set [0.0]')
	opt.add_argument('--min-aln-len', type=require_int_nonnegative,
		metavar='INT', default=0, help='minimum two-way alignment length of '
		'each set [0]')
	opt.add_argument('--name1', type=str, metavar='STR', default=None,
		help='identifier used for set1 output data [basename set1]')
	opt.add_argument('--name2', type=str, metavar='STR', default=None,
		help='identifier used for set2 output data [basename set2]')
	opt.add_argument('--refilter', default=False, action='store_true',
		help='skip alignment and re-filter previous search data with '
		'different cutoff parameters; overwrites output [off]')
	return parser.parse_args()

def require_int_nonnegative(x):
	try:
		if int(x) < 0 or '.' in str(x):
			sys.stderr.write('ERROR: {} must be a non-negative integer\n'.\
				format(x))
			sys.exit(1)
	except ValueError:
		sys.stderr.write('ERROR: {} must be an integer\n'.format(x))
		sys.exit(1)
	return int(x)

def require_float_0_to_100(x):
	try:
		x = float(x)
		if 0 < x > 100:
			sys.stderr.write('ERROR: {} must be in 0, 100 range\n'.format(x))
			sys.exit(1)
	except ValueError:
		sys.stderr.write('ERROR: {} must be a float\n'.format(x))
		sys.exit(1)
	return x

def require_float_0_to_1(x):
	try:
		x = float(x)
		if 0 < x > 1:
			sys.stderr.write('ERROR: {} must be in 0, 1.0 range\n'.format(x))
			sys.exit(1)
	except ValueError:
		sys.stderr.write('ERROR: {} must be a float\n'.format(x))
		sys.exit(1)
	return x

def require_dependency(dep):
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, dep)) and \
		not os.path.isdir(os.path.join(path, dep)):
			return True
	sys.stderr.write('ERROR: {} unavailable; not in $PATH\n'.format(dep))
	sys.exit(1)

def decompress_file(infile, outdir):
	uncompressed_file = os.path.basename(infile).rstrip('.gz')
	outfile = os.path.join(outdir, uncompressed_file)
	with gzip.open(infile, 'rb') as ifh, open(outfile, 'wb') as ofh:
		shutil.copyfileobj(ifh, ofh)
	return outfile

def genbank_to_faa(infile, outdir):
	faa_file = os.path.basename(infile).rsplit('.', 1)[0] + '.faa'
	outfile = os.path.join(outdir, faa_file)
	records = []
	for rec in SeqIO.parse(infile, 'genbank'):
		for feat in rec.features:
			if feat.type == 'CDS' and 'translation' in feat.qualifiers:
				if len(feat.qualifiers['translation']) == 1:
					new_rec = SeqRecord(Seq(feat.qualifiers['translation'][0],
						IUPAC.protein), id=feat.qualifiers['locus_tag'][0],
						description=rec.id)
					records.append(new_rec)
	with open(outfile, 'w') as ofh:
		SeqIO.write(records, ofh, 'fasta')
	return outfile

def get_seqlen(infile):
	seqlen = 0
	for rec in SeqIO.parse(infile, 'fasta'):
		seqlen += len(rec.seq)
	return seqlen

def verify_file_exists_and_nonempty(infile):
	if os.path.exists(infile):
		if os.stat(infile).st_size == 0:
			sys.stderr.write('ERROR: {} file empty\n'.format(infile))
			sys.exit(1)
	else:
		sys.stderr.write('ERROR: {} file absent\n'.format(infile))
		sys.exit(1)

def verify_protein_fasta_is_not_nucleotides(infile, max_acgt):
	c, total_length = Counter(), 0
	for rec in SeqIO.parse(infile, 'fasta'):
		c += Counter(rec.seq.upper())
		total_length += len(rec.seq)
	observed_acgt_fraction = 1. * (c['A']+c['C']+c['G']+c['T']) / total_length
	if observed_acgt_fraction > max_acgt:
		sys.stderr.write('ERROR: {} exceeds {} maximum allowed A, C, G, T '
			'composition. If {} is truly amino acids and not nucleotides, '
			'adjust the --max-ACGT value to compute AAI.\n'.format(
				observed_acgt_fraction, max_acgt, infile))
		sys.exit(1)

def count_fasta_records(infile):
	return sum([1 for ln in open(infile).readlines() if ln.startswith('>')])

def main():
	opts = parseArgs()
	require_dependency(opts.aligner)
	if opts.aligner == 'blastp':
		require_dependency('makeblastdb')

	# I/O handling
	if opts.cpus < 1:
		cpus = str(cpu_count())
	else:
		cpus = str(opts.cpus)
	dec_pts = opts.decimal_places
	psize, aln_frc = (0, 0), (0, 0, (0, 0))
	set1 = os.path.realpath(os.path.expanduser(opts.set1))
	set2 = os.path.realpath(os.path.expanduser(opts.set2))
	if opts.name1 is not None:
		b1 = opts.name1
	else:
		b1 = os.path.basename(opts.set1).rstrip('.gz').rsplit('.', 1)[0]
	if opts.name2 is not None:
		b2 = opts.name2
	else:
		b2 = os.path.basename(opts.set2).rstrip('.gz').rsplit('.', 1)[0]
	tmp = mkdtemp()
	if opts.outpath is not None:
		outpath = os.path.realpath(os.path.expanduser(opts.outpath))
	else:
		autocreated_dirname = 'AAI--' + strftime('%d%b%Y_%-I:%M:%S%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)
	if not os.path.exists(outpath):
		os.mkdir(outpath)

	if not opts.refilter:
		# Automatically handle gunzip compression and GenBank format input
		if set1.endswith('.gz'):
			set1 = decompress_file(set1, tmp)
		if set2.endswith('.gz'):
			set2 = decompress_file(set2, tmp)
		if set1.endswith(('.gbff', '.gbf', '.gbk', '.gb')):
			set1 = genbank_to_faa(set1, tmp)
		else:
			verify_protein_fasta_is_not_nucleotides(set1, opts.max_ACGT)
		if set2.endswith(('.gbff', '.gbf', '.gbk', '.gb')):
			set2 = genbank_to_faa(set2, tmp)
		else:
			verify_protein_fasta_is_not_nucleotides(set2, opts.max_ACGT)

		# Get input sequence lengths of proteomes
		psize = (get_seqlen(set1), get_seqlen(set2))

		# Execute bidirectional blast of proteins
		for s1, s2, sb1, sb2 in [(set1, set2, b1, b2), (set2, set1, b2, b1)]:
			aln = os.path.join(outpath, opts.aligner+'.' +sb1+','+sb2+ '.tab')
			ref_db = os.path.join(tmp, sb1)
			if opts.aligner == 'blastp':
				c1 = ['makeblastdb', '-in', s1, '-out', ref_db, '-dbtype',
					'prot']
				c2 = ['blastp', '-db', ref_db, '-query', s2, '-max_hsps', '1',
					'-max_target_seqs', '1', '-num_threads', cpus,
					'-out', aln, '-outfmt', '6 qseqid sseqid pident length'
					' mismatch gapopen qstart qend sstart send evalue'
					' bitscore qcovhsp gaps']
			elif opts.aligner == 'diamond':
				c1 = ['diamond', 'makedb', '--in', s1, '--db', ref_db]
				c2 = ['diamond', 'blastp', '--db', ref_db, '--query', s2,
					'--max-hsps', '1', '--max-target-seqs', '1',
					'--threads', cpus, '--out', aln, '--outfmt'] + \
					('6 qseqid sseqid pident length mismatch gapopen qstart'
					' qend sstart send evalue bitscore qcovhsp gaps').split()
			for cmd in (c1, c2):
				process = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
				_, err = process.communicate()
				if process.returncode != 0:
					sys.stderr.write(err)
					sys.stderr.write('ERROR: failed system call: {}\n'.\
						format(' '.join(cmd)))
					sys.exit(1)

	# Parse alignment output
	aai, tot, filt = [[] for _ in range(3)], Counter(), [[] for _ in range(3)]
	# idx0 set1 qry, idx1 set2 qry, idx2 twoway
	d = {} #filtered aln data
	# d keys are (qseqid, sseqid), vals are (%ident, gapless_aln_len)

	set1qry_base = os.path.join(outpath, opts.aligner + '.' + b2 + ',' + b1)
	set2qry_base = os.path.join(outpath, opts.aligner + '.' + b1 + ',' + b2)
	verify_file_exists_and_nonempty(set1qry_base + '.tab')
	verify_file_exists_and_nonempty(set2qry_base + '.tab')

	with open(set1qry_base + '.tab') as dat, \
	open(set1qry_base + '.filt.tab', 'w') as passed:
		for line in dat:
			tot[0] += 1
			l = line.split('\t')
			if float(l[2])>=opts.identity and int(l[3])>=opts.length \
			and float(l[11])>=opts.bitscore and float(l[12])>=opts.fraction:
				aai[0].append(float(l[2]))
				aln_len = int(l[3]) - int(l[13])
				d[(str(l[0]), str(l[1]))] = float(l[2]), aln_len
				filt[0].append(aln_len)
				passed.write(line)
	with open(set2qry_base + '.tab') as dat, \
	open(set2qry_base + '.filt.tab', 'w') as passed, \
	open(set2qry_base + '.filt.two-way.tab', 'w') as passed2:
		for line in dat:
			tot[1] += 1
			l = line.split('\t')
			if float(l[2])>=opts.identity and int(l[3])>=opts.length \
			and float(l[11])>=opts.bitscore and float(l[12])>=opts.fraction:
				aai[1].append(float(l[2]))
				aln_len = int(l[3]) - int(l[13])
				filt[1].append(aln_len)
				passed.write(line)
				pair = str(l[1]), str(l[0])
				if pair in d:
					aai[2].extend([float(l[2]), float(d[pair][0])])
					filt[2].append(aln_len)
					passed2.write(line)
	if any(len(l) == 0 for l in aai + filt):
		sys.stderr.write('ERROR: no alignments between sets.\n')
		sys.exit(1)

	# Calculate overall alignment (absolute and relative) values
	if not opts.refilter:
		set1_alnfrc = sum(filt[0]) / 1. / psize[0]
		set2_alnfrc = sum(filt[1]) / 1. / psize[1]
		bidirectional_alnfrc = (sum(filt[2]) / 1. / psize[0],
			sum(filt[2]) / 1. / psize[1])
		aln_frc = (set1_alnfrc, set2_alnfrc, bidirectional_alnfrc)
		tot = count_fasta_records(set1), count_fasta_records(set2)
	else:
		sys.stderr.write('WARNING: Total protein counts, which are the '
			'denominators in the second column output reported in the AAI '
			'stats tab file, were estimated from the blast alignment data '
			'files rather than from tallied input fragments. These values '
			'are likely a underestimation.\n'
			'WARNING: Proteome sizes and alignment fractions are reported as '
			'0 and 0.0. Turn off --refilter to get usable values.\n')

	# Optionally quit if low (absolute or relative) sequence alignment
	if opts.min_aln_len > sum(filt[2]):
		sys.stderr.write('ERROR: {} bi-directional (two-way) gapless '
			'alignment length is less than {}\n'.format(sum(filt[2]),
			opts.min_aln_len))
		sys.exit(1)
	if any(x < opts.min_aln_frac for x in aln_frc[2]):
		sys.stderr.write('ERROR: {} and/or {} bi-directional (two-way) '
			'gapless alignment fraction (relative to each input sequence '
			'length) is less than {}\n'.format(aln_frc[2][0], aln_frc[2][1],
			opts.min_aln_frac))
		sys.exit(1)

	# Write AAI data
	with open(os.path.join(outpath, 'aai.'+b1+','+b2+'.stats.tab'), 'w') \
	as ofh:
		ofh.write('Query_Sample(s)\tFiltered_Proteins\tAAI\tStDev\t'
			'Gapless_Aln\tInput_Size\tAln_Fraction\n')
		ofh.write('{}\t{}/{}\t{:.{p}f}%\t{:.{p}f}%\t{}\t{}\t{:.{p}f}\n'.\
			format(b1, len(filt[0]), tot[0], mean(aai[0]), std(aai[0]),
			sum(filt[0]), psize[0], aln_frc[0], p=dec_pts))
		ofh.write('{}\t{}/{}\t{:.{p}f}%\t{:.{p}f}%\t{}\t{}\t{:.{p}f}\n'.\
			format(b2, len(filt[1]), tot[1], mean(aai[1]), std(aai[1]),
			sum(filt[1]), psize[1], aln_frc[1], p=dec_pts))
		ofh.write('{},{}\t{}/{},{}/{}\t{:.{p}f}%\t{:.{p}f}%\t{},{}\t{},{}\t'
			'{:.{p}f},{:.{p}f}\n'.format(b1, b2, len(aai[2]) / 2, tot[0],
			len(aai[2]) / 2, tot[1], mean(aai[2]), std(aai[2]),
			sum(filt[2]), sum(filt[2]), psize[0], psize[1],
			aln_frc[2][0], aln_frc[2][1], p=dec_pts))
	shutil.rmtree(tmp)

if __name__ == '__main__':
	main()
