#!/usr/bin/env python


import os
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(
		description='extract record(s) from a multi-FastA file that match a query string')
	parser.add_argument('-i', '--infile', required=True,
		help='input multi-FastA file')
	parser.add_argument('-q', '--query', required=True,
		help='string to extract')
	parser.add_argument('-o', '--outfile', required=False, default=None,
		help='output file [`pwd`/<query>.fa]')
	return parser.parse_args()

def main():
	args = parseArgs()
	infile  = args.infile
	query   = args.query
	if args.outfile:
		outfile = args.outfile
	else:
		outfile = os.path.join(os.getcwd(), query + '.fa')
	if not os.path.exists(os.path.dirname(outfile)):
		os.mkdir(os.path.dirname(outfile))

	mfasta = SeqIO.parse(infile, 'fasta')
	query_match = []
	for record in mfasta:
		if query in record.description:
			query_match.append(record)

	if query_match:
		SeqIO.write(query_match, outfile, 'fasta')
	else:
		print '{} not found in {}'.format(query, infile)

if __name__ == '__main__':
	main()
