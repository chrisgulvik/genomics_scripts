#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqUtils import GC

def parseArgs():
	parser = ArgumentParser(description='Creates a FastA file from a '
	'GenBank file', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--input', required=True, metavar='FILE',
		help='input GenBank file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--output', metavar='FILE',
		help='output FastA file [basename.fa]')
	return parser.parse_args()


def main():
	opt = parseArgs()
	ifh = os.path.abspath(opt.input)
	if opt.output:
		ofh = os.path.abspath(opt.output)
	else:
		ofh = '.'.join(ifh.split('.')[:-1]) + '.fa'

	records = []
	with open(ifh) as i:
		for rec in SeqIO.parse(i, 'genbank'):
			# halt if nucleotides are absent rather than printing default 'N'
			if float(GC(rec.seq)) == 0:
				sys.exit('ERROR: {} appears to lack nucleotides'.format(ifh))
			records.append(rec)
	with open(ofh, 'w') as o:
		SeqIO.write(records, o, 'fasta')


if __name__ == '__main__':
	main()
