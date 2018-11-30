#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Converts a FastA file into a BED '
		'file, which is a tab-delimited list of sequence record names and '
		'their boundaries')
	parser.add_argument('-i', '--infile', required=True,
		help='input FastA Format file')
	parser.add_argument('-o', '--outfile', required=False, default=None,
		help='output BED Format file [stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.infile))

	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	with open(ifh) as infasta:
		for rec in SeqIO.parse(infasta, 'fasta'):
			ofh.write('{}\t0\t{}\n'.format(rec.id, len(rec)))

if __name__ == '__main__':
	main()
