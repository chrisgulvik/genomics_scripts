#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Converts a GenBank file containing '
		'nucleotide sequences into a Fasta Amino Acid format (FAA) file')
	parser.add_argument('-i', '--infile', required=True,
		help='input GenBank Format file <.gff||.gff3>')
	parser.add_argument('-o', '--outfile', required=False, default=None,
		help='output Fasta Amino Acid Format (.faa) file [stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.infile))
	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	with open(ifh) as i:
		for rec in SeqIO.parse(i, 'genbank') :
			for feat in rec.features:
				if feat.type == 'CDS' and \
				'translation' in feat.qualifiers:
					if len(feat.qualifiers['translation']) == 1:
						ofh.write('>{} {}\n{}\n'.format(
							feat.qualifiers['locus_tag'][0],
							rec.name,
							feat.qualifiers['translation'][0]))

if __name__ == '__main__':
	main()
