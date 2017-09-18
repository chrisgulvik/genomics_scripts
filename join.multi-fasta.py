#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from itertools import chain
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(add_help=False,
		description='Reads more than one multi-FastA file and joins '
		'sequences together (in the order given) based on a shared '
		'identifier in each sequence record. Works on amino and nucleic '
		'acids, with and without gaps. All sample identifiers in the first '
		'file given must be present in all other files, and those '
		'absent in the first file will be ignored. For example: '
		'>s1_00012 CCC >s2_00475 C-C and '
		'>s1_00037 GGGG >s2_00029 GAAG >s3_00029 GTTG '
		'become >s1 CCCGGGG >s2 C-CGAAG',
		epilog='NOTE: Prokka generates FAA files with locus_tags being first '
		'in each defline, which is extracted from each SeqRecord object '
		'property as a name (--seq-property name). Sample identifiers can '
		'then be extracted from each locus_tag by taking the prefix when '
		'split from an underscore (--seq-delim \'_\'). Defaults are setup to '
		'handle this.')
	req = parser.add_argument_group('Required')
	req.add_argument('mfa', metavar='FILE', nargs='+',
		help='input multi-FastA file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-d', '--discard', metavar='\'STR\'', type=str, default='',
		help='remove sites containing at least one specified character(s) '
		'such as gaps or ambiguities; comma-separate for more than one '
		'character and flank with apostrophes if necessary to santize special '
		'characters such as \'N,n,X,x,-\' [none]')
	opt.add_argument('-n', '--min-alleles', metavar='INT', type=int, default=1,
		help='minimum number of allele variants (after optional --discard '
		'filtering option) required per site to report in output; '
		'e.g., 2 removes monomorphic sites and 3 removes '
		'invariant and biallelic sites [1]')
	opt.add_argument('-o', '--outfile', metavar='FILE', default=None,
		help='joined output multi-FastA file [stdout]')
	opt.add_argument('-x', '--max-alleles', metavar='INT', type=int, default=1000,
		help='maximum number of allele variants (after optional --discard '
		'filtering option) required per site to report in output; '
		'e.g., 2 removes all multiallelic sites [None]')
	opt.add_argument('--seq-delim', metavar='\'STR\'', type=str,
		default='\'_\'',
		help='extracted identifiers (from --seq-property) are split on this '
		'delimiter, and only the first item from this is used as a sample '
		'identifier for joining sequences. Flank with apostrophes if '
		'necessary to santize special characters. More than one character to '
		'singly match and split off is permitted. Set to \'\' if each '
		'identifier is not a locus_tag. [\'_\']')
	opt.add_argument('--seq-property',
		choices=['description', 'id', 'name'], default='name',
		help='SeqRecord object property to first extract identifiers from '
		'[name]')
	opt.add_argument('--sort-files', action='store_true',
		default=False,
		help='toggle on alphanumeric sorting of input sample files to '
		're-orient sequence order in the output [off]')
	return parser.parse_args()

def mfa_to_dic(seq_delim, seq_property, infile, keys_only=False):
	if len(seq_delim) > 0:
		if seq_property == 'name':
			d = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'),
				key_function=lambda rec: \
				rec.name.split(seq_delim)[0])
		elif seq_property == 'id':
			d = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'),
				key_function=lambda rec: \
				rec.id.split(seq_delim)[0])
		elif seq_property == 'description':
			d = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'),
				key_function=lambda rec: \
				rec.description.split(seq_delim)[0])
	else:
		if seq_property == 'name':
			d = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'),
				key_function=lambda rec: rec.name)
		elif seq_property == 'id':
			d = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'),
				key_function=lambda rec: rec.id)
		elif seq_property == 'description':
			d = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'),
				key_function=lambda rec: rec.description)
	if keys_only:
		d = sorted(list(d.keys()))
	return d

def main():
	opt = parseArgs()

	# Input file handling
	ifs = [os.path.abspath(os.path.expanduser(f)) for f in opt.mfa]
	if opt.sort_files:
		ifs.sort()

	# Get sample IDs from first file's deflines
	qp = opt.seq_property
	qd = opt.seq_delim.lstrip('\'').rstrip('\'')
	ids = mfa_to_dic(qd, qp, ifs[0], keys_only = True)

	# Capture all sequences for each sample ID from each input file
	d = {}
	for idx_file, file in enumerate(ifs):
		fa = mfa_to_dic(qd, qp, file)
		for sample in ids:
			d['{}_{}'.format(idx_file, sample)] = str(fa[sample].seq)

	# Prepare output
	o = []
	for sample in ids:
		seq = ''
		for i in range(len(ifs)):
			seq += d['{}_{}'.format(i, sample)]
		o.extend(['>' + str(sample), str(seq)])
	cnt_init = cnt_keep = len(o[1])


	# Optional site filtering
	if len(opt.discard) > 0 or opt.min_alleles != 1 or opt.max_alleles != 1000:
		unwanted = opt.discard.lstrip('\'').rstrip('\'').split(',')
		keep = []

		# Transpose sequences
		for site in zip(*(o[1::2])):
			alleles = set(site)
			filt_alleles = filter(lambda a: a not in unwanted, alleles)

			# Filter allelle contents
			if len(alleles) == len(filt_alleles):
				# Filter allelle quantities
				if len(alleles) >= opt.min_alleles and \
				len(alleles) <= opt.max_alleles:
					keep.append(''.join(site))
		cnt_keep = len(keep)

		# Interleave deflines with sequences
		o = list(chain.from_iterable(zip(o[0::2], zip(*keep))))

	# Report tallies
	for i, s in ((len(ids), 'sequence records found'),
		(cnt_init, 'initial sites per sequence'),
		(cnt_init - cnt_keep, 'sites discarded per sequence'),
		(cnt_keep, 'final sites per sequence')):
		sys.stderr.write('  {} {}\n'.format(i, s))

	# Write single output file
	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	for ln in o:
		ofh.write('{}\n'.format(''.join(ln)))

if __name__ == '__main__':
	main()
