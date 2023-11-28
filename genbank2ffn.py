#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Converts a GenBank file containing '
		'nucleotide sequences into a FastA Format of genes (FFN) file')
	parser.add_argument('-i', '--infile', required=True,
		help='input GenBank Format file <.gbff||.gbk>')
	parser.add_argument('-o', '--outfile', required=False, default=None,
		help='output FastA Format of gene sequences (.ffn) file [stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	o = []
	with open(os.path.abspath(os.path.expanduser(opt.infile))) as ifh:
		for rec in SeqIO.parse(ifh, 'genbank'):
			for feat in rec.features:
				if feat.type == 'gene':
					seq = str(feat.location.extract(rec.seq))
					tag = feat.qualifiers['locus_tag'][0]
					o.append('>{}\n{}\n'.format(tag, seq))

	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	for rec in o:
		ofh.write(rec)

if __name__ == '__main__':
	main()
