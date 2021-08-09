#!/usr/bin/env python


import gzip
import os
import re
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(
		description='extract record(s) from a GenBank file that contain one '
		'or more query and output in FastA format', add_help=False,
		epilog='NOTE: GenBank Feature Table Definition is described at '
		'http://www.insdc.org/files/feature_table.html')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input GenBank file, optionally gunzip compressed')
	req.add_argument('-q', '--query', required=True, metavar='STR',
		help='string to search each locus (record) name')
	req.add_argument('-r', '--range', required=True, metavar='INT:INT|INT-INT',
		help='range to extract within the locus')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', default=None, metavar='FILE',
		help='FastA output [default: stdout]')
	opt.add_argument('--search-type', default='contains',
		choices=['contains', 'full_exact'],
		help='type of query match [default: contains]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	infile = os.path.abspath(os.path.expanduser(opt.infile))
	if infile.endswith('.gz'):
		infile = gzip.open(infile)
	qry = opt.query
	positions = opt.range

	# Quick checks of positions
	inv_pattern = re.compile('[^0-9:-]')
	if any([inv_pattern.search(x) for x in positions]):
		sys.stderr.write('ERROR: sites must only contain [0-9:-]\n')
		sys.exit(1)
	if positions.startswith('0'):
		sys.stderr.write('ERROR: 0 position prefix found; site extraction is'
			' 1-based\n')
		sys.exit(1)
	if ':' in positions:
		start, stop = [int(x) for x in positions.split(':')]
	elif '-' in positions:
		start, stop = [int(x) for x in positions.split('-')]
	if start >= stop:
		sys.stderr.write('ERROR: {} range must have second value'
			' exceed first.\n'.format(positions))
		sys.exit(1)

	# Find the contig (locus or record) of interest
	query_match = []
	for rec in SeqIO.parse(infile, 'genbank'):
		if opt.search_type == 'contains':
			if qry in rec.name:
				query_match.append(rec)
		elif opt.search_type == 'full_exact':
			if qry == rec.name:
				query_match.append(rec)
	if len(query_match) == 0:
		sys.stderr.write('ERROR: {} absent in {}\n'.format(qry, infile))
		sys.exit(1)
	elif len(query_match) > 1:
		sys.stderr.write('ERROR: found >1 {} in {}\n'.format(qry, infile))
		sys.exit(1)

	# Slice GenBank record
	rec = query_match[0]
	if any(v > len(rec.seq) for v in (start, stop)):
		sys.stderr.write('ERROR: {} range larger than {} record'
			' ({} sites).\n'.format(positions, rec.name, len(rec.seq)))
		sys.exit(1)
	sliced_rec = rec[start-1:stop]

	# Warn about partial feature(s)
	feats = [feat for feat in rec.features if feat.type == 'CDS']
	desired = set(range(start-1, stop))
	for f in feats:
		span = set(range(f.location.start.position, f.location.end.position))
		intersection = span & desired
		min_intersection = min(intersection)
		if len(intersection) < len(desired):
			# NOTE: more complex comparison needed to handle multi-gene range  
			sys.stderr.write('INFO: partial incomplete sequence match to\n\n')
			sys.stderr.write('{}\n\n'.format(f))

	# Output genbank format
	if opt.outfile:
		outfile = os.path.abspath(os.path.expanduser(opt.outfile))
		SeqIO.write(sliced_rec, outfile, 'genbank')
	else:
		SeqIO.write(sliced_rec, sys.stdout, 'genbank')

if __name__ == '__main__':
	main()
