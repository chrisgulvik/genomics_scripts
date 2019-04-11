#!/usr/bin/env python


import gzip
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
		help='input GenBank file (optionally gunzip compressed)')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--output', metavar='FILE',
		help='output FastA file [basename.fa]')
	return parser.parse_args()


def main():
	opt = parseArgs()
	ifh = os.path.realpath(os.path.expanduser(opt.input))
	if opt.output:
		ofh = os.path.realpath(os.path.expanduser(opt.output))
	else:
		ofh = os.path.basename(ifh).rstrip('.gz').rsplit('.', 1)[0] + '.fa'

	records = []
	if ifh.endswith('.gz'):
		ifh = gzip.open(ifh)	
	for rec in SeqIO.parse(ifh, 'genbank'):
		# halt if nucleotides are absent rather than printing default 'N'
		if float(GC(rec.seq)) == 0:
			sys.stderr.write('ERROR: {} appears to lack nucleotides\n'.\
				format(ifh))
			sys.exit(1)
		records.append(rec)
	with open(ofh, 'w') as o:
		SeqIO.write(records, o, 'fasta')


if __name__ == '__main__':
	main()
