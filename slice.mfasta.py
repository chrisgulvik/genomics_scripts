#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Extracts sequence site(s) from a'
	' multi-record FastA file', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		type=str, help='input FastA sequence file')
	req.add_argument('-s', '--sites', nargs='+', metavar='INT|INT:INT',
		help='position(s) to extract from each sequence record;'
		' colon or hyphen separate for ranges')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', required=False, metavar='FILE',
		default=None, help='FastA-formatted output [stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.realpath(os.path.expanduser(opt.infile))
	sites = []
	for site in opt.sites:
		if ':' in site:
			sites.append(':'.join([int(n) + 1 for n in site.split(':')]))
		elif '-' in site:
			sites.append(':'.join([int(n) + 1 for n in site.split(':')]))
		else:
			sites.append(int(site) + 1)

	# Find sites in each seq record
	out = []
	for rec in SeqIO.parse(ifh, 'fasta'):
		seq = ''
		for site in sites:
			seq += str(rec[site].seq)

	# Output
	if opt.outfile is not None:
		ofh = os.path.abspath(os.path.expanduser(opt.outfile))
		SeqIO.write(out, ofh, 'fasta')
	else:
		SeqIO.write(out, sys.stdout, 'fasta')

if __name__ == '__main__':
	main()