#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Slices out a single sequence from a '
	'FastA file', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		type=str, help='input FastA sequence file')
	req.add_argument('-d', '--defline', required=True, metavar='STR',
		type=str, help='header/defline name for target sequence record to '
		'slice out sequences from')
	req.add_argument('-b', '--begin', required=True, metavar='INT', type=int,
		help='initial position to extract from the specified sequence record')
	req.add_argument('-e', '--end', required=True, metavar='INT', type=int,
		help='final position to extract from the specified sequence record')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-m', '--defline-query-method', required=False,
		type=str, default='substr', choices=['full', 'substr'],
		help='search method for header/defline name [substr]')
	opt.add_argument('-o', '--outfile', required=False, metavar='FILE',
		default=None, help='FastA-formatted sliced output [stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.abspath(opt.infile)
	get = (opt.defline, opt.begin-1, opt.end)

	# Find seq record
	mfa = SeqIO.parse(ifh, 'fasta')
	fnd = []
	if opt.defline_query_method == 'full':
		for rec in mfa:
			if str(get[0]) == rec.description:
				fnd.append(rec)
	elif opt.defline_query_method == 'substr':		
		for rec in mfa:
			if str(get[0]) in rec.description:
				fnd.append(rec)
	else:
		sys.stderr.write('ERROR: unsupported {} search method; argparse '
			'should have prevented this\n'.format(opt.defline_query_method))
		sys.exit(1)

	# Halt if anything but one seq record match found
	if len(fnd) > 1:
		sys.stderr.write('ERROR: >1 defline match to {}\n'.format(get[0]))
		sys.exit(1)
	elif len(fnd) == 0:
		sys.stderr.write('ERROR: {} absent from deflines\n'.format(get[0]))
		sys.exit(1)

	# Slice seq record
	out = fnd[0][get[1]:get[2]]

	# Output
	if opt.outfile is not None:
		ofh = os.path.abspath(os.path.expanduser(opt.outfile))
		SeqIO.write(out, ofh, 'fasta')
	else:
		SeqIO.write(out, sys.stdout, 'fasta')

if __name__ == '__main__':
	main()