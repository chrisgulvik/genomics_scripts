#!/usr/bin/env python


import argparse
import os
import sys

def parseArgs():
	parser = argparse.ArgumentParser(description='extracts alignment '
		'coordinates in a Mauve or Parsnp backbone file and saves as '
		'a BED format', add_help=False, epilog='NOTE: unaligned sites '
		'can be identified with bedtools complement')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', metavar='FILE',
		required=True, help='input backbone file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--col', type=int, default=1, metavar='INT',
		help='start column for coordinate extraction [1]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', default=None, metavar='FILE',
		help='output BED file [stdout]')
	return parser.parse_args()

def main():
	args   = parseArgs()
	infile = os.path.abspath(os.path.expanduser(args.infile))
	col    = args.col
	if args.outfile is None:
		out = sys.stdout
	else:
		out = open(os.path.abspath(os.path.expanduser(args.outfile)), 'w')

	with open(infile) as ifh:
		hdr = next(ifh).rstrip('\n').split('\t')
		chrom = hdr[1-col].rstrip('_start').lstrip('>')
		for ln in ifh:
			l = ln.rstrip('\n').split('\t')
			out.write('{}\t{}\t{}\n'.format(
				chrom, int(l[1-col].lstrip('-')), int(l[col].lstrip('-'))))

if __name__ == '__main__':
	main()
