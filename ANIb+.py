#!/usr/bin/env python

__version__ = '1.1.1b'


import gzip
import os
import shutil
import sys
import subprocess as sp
from argparse import ArgumentParser
from collections import deque, Counter
from itertools import islice
from glob import glob
from multiprocessing import cpu_count
from tempfile import mkdtemp
from time import strftime

from Bio import SeqIO
from Bio.SeqUtils import GC
from numpy import mean, std

def parseArgs():
	parser = ArgumentParser(description='Computes the average nucleotide '
		'identity (ANI) between two nucleic acid sequence sets with BLASTn',
		add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-1', '--set1', required=True, metavar='FILE',
		help='first input FastA or GenBank file '
		'(optionally gunzip compressed)')
	req.add_argument('-2', '--set2', required=True, metavar='FILE',
		help='second input FastA or GenBank file '
		'(optionally gunzip compressed)')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--cpus', type=require_int_nonnegative,
		metavar='INT', default=0, help='number of CPUs [all]')
	opt.add_argument('-f', '--fraction', type=require_float_0_to_100,
		metavar='FLOAT', default=70.0, help='minimum fragment alignment '
		'length percentage [70.0]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-i', '--identity', type=require_float_0_to_100,
		metavar='FLOAT', default=30.0, help='minimum fragment alignment '
		'identity [30.0]')
	opt.add_argument('-l', '--length', type=require_int_nonnegative,
		metavar='INT', default=0, help='minimum fragment alignment length '
		'(incl gaps) [0]')
	opt.add_argument('-o', '--outpath', type=str, metavar='PATH',
		default=None, help='output directory [ANI--<date>_<time>]')
	opt.add_argument('-s', '--step-size', type=require_int_positive,
		metavar='INT', default=200, help='number of bases to shift during '
		'sequence fragmentation prior to alignments; to turn off steps to '
		'speed up, set to same length as the fragment size [200]')
	opt.add_argument('-v', '--version', action='version',
		version='%(prog)s v{}'.format(__version__))
	opt.add_argument('-w', '--fragment-size', type=require_int_positive,
		metavar='INT', default=1000, help='fragment lengths to slice '
		'nucleotides into prior to alignments [1000]')
	opt.add_argument('--bed', default=False, action='store_true',
		help='produce BED files listing coordinates (0-based) used in ANI '
		'calculation; requires Bedtools [off]')
	opt.add_argument('--decimal-places', type=require_int_nonnegative,
		metavar='INT', default=3, help='precision decimal points to round '
		'ANI values [3]')
	opt.add_argument('--fill', type=str, metavar='CHAR', default='',
		help='character to add to end of fragments that are less than the '
		'window size to force them to be the same length; turns on '
		'--keep-small-frags [None]')
	opt.add_argument('--keep-small-frags', default=False,
		action='store_true', help='include the final fragment (less than the '
		'window size) of each sequence record during alignment [off]')
	opt.add_argument('--min-ACGT', type=require_float_0_to_1, metavar='FLOAT',
		default=.97, help='minimum fraction of ACGT in each input set [0.97]')
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
		help='skip BLASTn and re-filter previous alignment data with '
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

def require_int_positive(x):
	try:
		if int(x) < 1 or '.' in str(x):
			sys.stderr.write('ERROR: {} must be a positive integer\n'.\
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

def genbank_to_fasta(infile, outdir):
	fasta_file = os.path.basename(infile).rsplit('.', 1)[0] + '.fa'
	outfile = os.path.join(outdir, fasta_file)
	records = []
	for rec in SeqIO.parse(infile, 'genbank'):
		if float(GC(rec.seq)) == 0:
			sys.stderr.write('ERROR: {} appears to lack nucleotides'.\
				format(rec.name))
			sys.exit(1)
		records.append(rec)
	with open(outfile, 'w') as ofh:
		SeqIO.write(records, ofh, 'fasta')
	return outfile

def get_seqlen_and_filter_ambiguous(infile, min_ACGT):
	seqlen, c = 0, Counter()
	for rec in SeqIO.parse(infile, 'fasta'):
		seqlen += len(rec.seq)
		c += Counter(list(rec.seq.upper()))
	ACGT_fraction = float(c['A'] + c['C'] + c['G'] + c['T']) / seqlen
	if ACGT_fraction < min_ACGT:
		sys.stderr.write('ERROR: low composition ({:.1f}%) of unambiguous '
			'nucleotides in {}\n'.format(ACGT_fraction * 100, 
			os.path.basename(infile)))
		sys.exit(1)
	if ACGT_fraction != 1.:
		sys.stderr.write('WARNING: {:.6f}% ambiguous nucleotides detected in '
			'{}\n'.format(100 - ACGT_fraction*100, os.path.basename(infile)))
	return seqlen

def fragment(seq, win, step, fill):
	iters = iter(seq)
	q = deque(islice(iters, win), maxlen=win)
	q.extend(fill for _ in range(win-len(q)))
	while True:
		yield q
		q.append(next(iters))
		q.extend(next(iters, fill) for _ in range(step-1))

def parse_coords(seq_id, aln_coords, outfile):
	hdr_name   = seq_id.split('__frg')[0]
	hdr_coords = [int(s) for s in seq_id.split('__pos')[-1].split('-')]
	aln_coords = [int(s) for s in aln_coords]
	if aln_coords[1] < aln_coords[0]:
		aln_coords[1], aln_coords[0] = aln_coords
	start = aln_coords[0] + hdr_coords[0]
	end   = aln_coords[1] + hdr_coords[0]
	with open(outfile, 'a') as ofh:
		ofh.write('{}\t{}\t{}\n'.format(hdr_name, start, end))

def count_fasta_records(infile):
	return sum([1 for ln in open(infile).readlines() if ln.startswith('>')])

def main():
	opts = parseArgs()
	require_dependency('blastn')
	require_dependency('makeblastdb')
	if opts.bed:
		require_dependency('bedtools')

	# I/O handling
	if opts.cpus < 1:
		cpus = str(cpu_count())
	else:
		cpus = str(opts.cpus)
	dec_pts = opts.decimal_places
	gsize, aln_frc = (0, 0), (0, 0, (0, 0))
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
		autocreated_dirname = 'ANI--' + strftime('%d%b%Y_%-I:%M:%S%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)
	if not os.path.exists(outpath):
		os.mkdir(outpath)

	if not opts.refilter:
		keep_small_frags = opts.keep_small_frags
		if len(opts.fill) > 0:
			keep_small_frags = True

		# Automatically handle gunzip compression and GenBank format input
		if set1.endswith('.gz'):
			set1 = decompress_file(set1, tmp)
		if set2.endswith('.gz'):
			set2 = decompress_file(set2, tmp)
		if set1.endswith(('.gbff', '.gbf', '.gbk', '.gb')):
			set1 = genbank_to_fasta(set1, tmp)
		if set2.endswith(('.gbff', '.gbf', '.gbk', '.gb')):
			set2 = genbank_to_fasta(set2, tmp)
		
		# Get input sequence lengths prior to fragmentation and
		# quit if either input set lacks primarily A,C,G,T nucleotides
		gsize = (get_seqlen_and_filter_ambiguous(set1, opts.min_ACGT),
			get_seqlen_and_filter_ambiguous(set2, opts.min_ACGT))

		# Fragment input sequences
		if1 = os.path.join(tmp, b1 + '.frags')
		if2 = os.path.join(tmp, b2 + '.frags')
		for f, s in [(if1, set1), (if2, set2)]:
			with open(f, 'w') as ofh:
				for rec in SeqIO.parse(s, 'fasta'):
					coord, i = 0, 1
					for frag in fragment(rec.seq, opts.fragment_size,
					opts.step_size, opts.fill):
						frag_seq = ''.join(list(frag))
						if keep_small_frags or \
						len(frag_seq) == opts.fragment_size:
							defln = ('{}__frg{}__pos{}-{}').format(rec.id, i,
								coord, coord + opts.fragment_size)
							ofh.write('>{}\n{}\n'.format(defln, frag_seq))
							i += 1
							coord += opts.step_size

		# Align fragments
		for s1, s2, sb1, sb2 in [(if1, if2, b1, b2), (if2, if1, b2, b1)]:
			aln = os.path.join(outpath, 'blast.' + sb1 + ',' + sb2 + '.tab')
			c1 = ['makeblastdb', '-in', s1, '-out', s1, '-dbtype', 'nucl']
			c2 = ['blastn', '-db', s1, '-query', s2, '-dust', 'no',
				'-max_hsps', '1', '-max_target_seqs', '1',
				'-num_threads', cpus, '-task', 'blastn',
				'-outfmt', '6 qseqid sseqid pident length mismatch gapopen'
				' qstart qend sstart send evalue bitscore qcovhsp gaps',
				'-out', aln]
			for cmd in (c1, c2):
				process = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
				_, err = process.communicate()
				if process.returncode != 0:
					sys.stderr.write(err)
					sys.stderr.write('ERROR: failed system call: {}\n'.\
						format(' '.join(cmd)))
					sys.exit(1)

	# Parse blast output
	ani, tot, filt = [[] for _ in range(3)], Counter(), [[] for _ in range(3)]
	# idx0 set1 qry, idx1 set2 qry, idx2 twoway
	d = {} #filtered aln data
	# d keys are (qseqid, sseqid), vals are (%ident, gapless_aln_len)
	b = os.path.join(outpath, 'blast.' + b2 + ',' + b1) #set1 as query
	with open(b + '.tab') as dat, \
	open(b + '.filt.tab', 'w') as passed:
		for line in dat:
			tot[0] += 1
			l = line.split('\t')
			if float(l[2])>=opts.identity and int(l[3])>=opts.length \
			and float(l[12])>=opts.fraction:
				ani[0].append(float(l[2]))
				aln_len = int(l[3]) - int(l[13])
				d[(str(l[0]), str(l[1]))] = float(l[2]), aln_len
				filt[0].append(aln_len)
				passed.write(line)
				if opts.bed:
					parse_coords(l[0], (l[6], l[7]), os.path.join(tmp,
						'coords.'+b1+'frags.'+b2+'ref,'+b1+'qry.bed'))
					parse_coords(l[1], (l[8], l[9]), os.path.join(tmp,
						'coords.'+b2+'frags.'+b2+'ref,'+b1+'qry.bed'))
	b = os.path.join(outpath, 'blast.' + b1 + ',' + b2) #set2 as query
	with open(b + '.tab') as dat, \
	open(b + '.filt.tab', 'w') as passed, \
	open(b + '.filt.two-way.tab', 'w') as passed2:
		for line in dat:
			tot[1] += 1
			l = line.split('\t')
			if float(l[2])>=opts.identity and int(l[3])>=opts.length \
			and float(l[12])>=opts.fraction:
				ani[1].append(float(l[2]))
				aln_len = int(l[3]) - int(l[13])
				filt[1].append(aln_len)
				passed.write(line)
				if opts.bed:
					parse_coords(l[0], (l[6], l[7]), os.path.join(tmp,
						'coords.'+b2+'frags.'+b1+'ref,'+b2+'qry.bed'))
					parse_coords(l[1], (l[8], l[9]), os.path.join(tmp,
						'coords.'+b1+'frags.'+b1+'ref,'+b2+'qry.bed'))
				pair = str(l[1]), str(l[0])
				if pair in d:
					ani[2].extend([float(l[2]), float(d[pair][0])])
					filt[2].append(aln_len)
					passed2.write(line)
					if opts.bed:
						parse_coords(l[0], (l[6], l[7]), os.path.join(tmp,
							'coords.'+b2+'frags.two-way.bed'))
						parse_coords(l[1], (l[8], l[9]), os.path.join(tmp,
							'coords.'+b1+'frags.two-way.bed'))
	if any(len(l) == 0 for l in ani + filt):
		sys.stderr.write('ERROR: no alignments between sets.\n')
		sys.exit(1)

	# Calculate overall alignment (absolute and relative) values
	if not opts.refilter:
		coverage = 1. * opts.fragment_size / opts.step_size
		set1_alnfrc = sum(filt[0]) / coverage / gsize[0]
		set2_alnfrc = sum(filt[1]) / coverage / gsize[1]
		bidirectional_alnfrc = (sum(filt[2]) / coverage / gsize[0],
			sum(filt[2]) / coverage / gsize[1])
		aln_frc = (set1_alnfrc, set2_alnfrc, bidirectional_alnfrc)
		tot = count_fasta_records(if1), count_fasta_records(if2)
	else:
		sys.stderr.write('WARNING: Total fragment counts, which are the '
			'denominators in the second column output reported in the ANI '
			'stats tab file, were estimated from the blast alignment data '
			'files rather than from tallied input fragments. These values '
			'are likely a underestimation.\n'
			'WARNING: Genome sizes and alignment fractions are reported as '
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

	# Write ANI data
	with open(os.path.join(outpath, 'ani.'+b1+','+b2+'.stats.tab'), 'w') \
	as ofh:
		ofh.write('Query_Sample(s)\tFiltered_Frags\tANI\tStDev\t'
			'Gapless_Aln\tInput_Size\tAln_Fraction\n')
		ofh.write('{}\t{}/{}\t{:.{p}f}%\t{:.{p}f}%\t{}\t{}\t{:.{p}f}\n'.\
			format(b1, len(filt[0]), tot[0], mean(ani[0]), std(ani[0]),
			sum(filt[0]), gsize[0], aln_frc[0], p=dec_pts))
		ofh.write('{}\t{}/{}\t{:.{p}f}%\t{:.{p}f}%\t{}\t{}\t{:.{p}f}\n'.\
			format(b2, len(filt[1]), tot[1], mean(ani[1]), std(ani[1]),
			sum(filt[1]), gsize[1], aln_frc[1], p=dec_pts))
		ofh.write('{},{}\t{}/{},{}/{}\t{:.{p}f}%\t{:.{p}f}%\t{},{}\t{},{}\t'
			'{:.{p}f},{:.{p}f}\n'.format(b1, b2, len(ani[2]) / 2, tot[0],
			len(ani[2]) / 2, tot[1], mean(ani[2]), std(ani[2]),
			sum(filt[2]), sum(filt[2]), gsize[0], gsize[1],
			aln_frc[2][0], aln_frc[2][1], p=dec_pts))

	# Optionally produce BED files of coordinates used for ANI calculations
	if opts.bed:
		beds = glob(os.path.join(tmp, 'coords.*.bed'))
		for f in beds:
			srt = sp.Popen(['bedtools', 'sort', '-i', f], stdout=sp.PIPE)
			mrg = sp.Popen(['bedtools', 'merge', '-i', 'stdin'],
				stdin=srt.stdout, stdout=sp.PIPE)
			with open(os.path.join(outpath, os.path.basename(f)), 'w') as ofh:
				ret = sp.Popen(['sort', '-V'], stdin=mrg.stdout, stdout=ofh)
	shutil.rmtree(tmp)

if __name__ == '__main__':
	main()
