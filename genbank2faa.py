#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Creates a Fasta Amino Acid (FAA) '
		'formatted file of protein sequences from a GenBank flat file format',
		add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input GenBank file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', default=None,
		help='output Fasta Amino Acid (.faa) format file [stdout]')
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
