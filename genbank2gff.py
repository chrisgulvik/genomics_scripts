#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from BCBio import GFF  #pip install bcbio-gff
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Converts a GenBank file containing '
		'nucleotide sequences into a General Feature Format (GFF) file')
	parser.add_argument('-i', '--infile', required=True,
		help='input GenBank Format file <.gff||.gff3>')
	parser.add_argument('-o', '--outfile', required=False, default=None,
		help='output General Feature Format (.gff or .gff3) file [stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.infile))
	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	with open(ifh) as i:
		GFF.write(SeqIO.parse(i, 'genbank'), ofh)

if __name__ == '__main__':
	main()
