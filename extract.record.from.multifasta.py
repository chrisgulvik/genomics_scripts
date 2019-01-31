#!/usr/bin/env python


import gzip
import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(
		description='extract record(s) from a multi-record FastA file that '
		'have a query string match to defline(s)', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input multi-FastA file, optionally gunzip compressed')
	req.add_argument('-q', '--query', required=True, type=str, metavar='STR',
		help='string to search deflines')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', default=None, metavar='FILE',
		help='FastA output [stdout]')
	opt.add_argument('-y', '--search-type', default='contains',
		choices=['contains', 'full_exact'], help='type of query match '
		'[default: contains]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	infile = os.path.abspath(os.path.expanduser(opt.infile))
	query  = opt.query

	if infile.endswith('.gz'):
		infile = gzip.open(infile)
	mfasta = SeqIO.parse(infile, 'fasta')
	query_match = []
	for record in mfasta:
		if opt.search_type == 'contains':
			if query in str(record.description):
				query_match.append(record)
		elif opt.search_type == 'full_exact':
			if query == str(record.description):
				query_match.append(record)

	if len(query_match) == 0:
		sys.stderr.write('ERROR: {} absent from deflines\n'.format(query))
		sys.exit(1)

	if opt.outfile:
		outfile = os.path.abspath(os.path.expanduser(opt.outfile))
		SeqIO.write(query_match, outfile, 'fasta')
	else:
		SeqIO.write(query_match, sys.stdout, 'fasta')

if __name__ == '__main__':
	main()