#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqUtils import GC

def parseArgs():
	parser = ArgumentParser(add_help=False,
		description='reports statistics of multi-FastA file per sequence '
		'record in tab-delimited format: Name, Length, %GC, %CT, # Gaps, # Ns')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input multi-FastA file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-f', '--filename', action='store_true', default=False,
		help='add input filename in front of first data column [off]')
	opt.add_argument('-g', '--nogaps', action='store_true', default=False,
		help='prior to statistics calculations, remove \'-\' gaps [off]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-n', '--noambiguous', action='store_true', default=False,
		help='prior to statistics calculations, remove N,n [off]')
	opt.add_argument('-o', '--outfile', metavar='FILE',
		help='tab-delimited output file [stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()

	# I/O handling
	infile = os.path.abspath(os.path.expanduser(opt.infile))
	mfasta = SeqIO.parse(infile, 'fasta')

	o = []
	for rec in mfasta:
		# Sequence and defline handling
		seq  = str(rec.seq).upper()
		if opt.nogaps:
			seq = seq.replace('-', '')
		if opt.noambiguous:
			seq = seq.replace('N', '')
		defln = rec.description
		if opt.filename:
			defln = '{}\t{}'.format(infile, defln)

		# Statistics
		seqlen = len(seq)
		ct = 100.0 * (seq.count('C') + seq.count('G')) / seqlen
		gc = GC(seq)
		gaps = seq.count('-')
		ambig = seq.count('N')
		o.append('{}\t{}\t{:.4f}\t{:.4f}\t{}\t{}'.format(
			defln, seqlen, gc, ct, gaps, ambig))

	# Write single output file
	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	for ln in o:
		ofh.write('{}\n'.format(''.join(ln)))

if __name__ == '__main__':
	main()
