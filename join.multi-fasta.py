#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
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
	req.add_argument('mfa0', metavar='FILE',
		help='input multi-FastA file')
	req.add_argument('mfa1', metavar='FILE', nargs='+',
		help='input multi-FastA file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', metavar='FILE', default=None,
		help='joined output multi-FastA file [stdout]')
	opt.add_argument('--seq-delim', metavar='\'STR\'', type=str,
		default='\'_\'',
		help='extracted identifiers (from --seq-property) are split on this '
		'delimiter, and only the first item from this is used as a sample '
		'identifier for joining sequences. Flanking apostrophes required. '
		'More than one character permitted. Set to \'\' if each identifier '
		'is not a locus_tag. [\'_\']')
	opt.add_argument('--seq-property',
		choices=['description', 'id', 'name'], default='name',
		help='SeqRecord object property to first extract identifiers from '
		'[name]')
	opt.add_argument('--sort-files', action='store_true',
		default=False,
		help='toggle on alphanumeric sorting of input sample files to '
		're-orient sequence order in the output [off]')
	return parser.parse_args()

def mfa_to_dic(seq_delim, seq_property, infile, keys_only = False):
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
	ifs = [os.path.abspath(os.path.expanduser(f)) for f in opt.mfa1]
	ifs.insert(0, os.path.abspath(os.path.expanduser(opt.mfa0)))
	if opt.sort_files:
		ifs.sort()

	# Get sample IDs from first file's deflines
	qp = opt.seq_property
	qd = opt.seq_delim[1:-1]
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

	# Write single output file
	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	for ln in o:
		ofh.write('{}\n'.format(ln))

if __name__ == '__main__':
	main()
