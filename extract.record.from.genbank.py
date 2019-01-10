#!/usr/bin/env python


import gzip
import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(
		description='extract record(s) from a GenBank file that contain a'
		' query string and output in FastA format', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True,
		help='input GenBank file')
	req.add_argument('-q', '--query', required=True,
		help='string to search deflines')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-f', '--query-feature', default='CDS', metavar='STR',
		help='genbank feature type to search in, e.g., CDS, gene, rRNA, '
		'source, tRNA, misc_feature')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', required=False, default=None,
		help='FastA output [stdout]')
	opt.add_argument('-u', '--query-qualifier', default='inference',
		metavar='STR', help='qualifier term within each genbank feature to '
		'search in, e.g., locus_tag, inference, codon_start, product, '
		'transl_table, translation')
	return parser.parse_args()

def main():
	opt = parseArgs()
	infile = os.path.abspath(os.path.expanduser(opt.infile))
	inbase = os.path.splitext(os.path.basename(infile))[0]
	query_term   = opt.query
	query_feat   = opt.query_feature
	query_qualif = opt.query_qualifier

	if infile.endswith('.gz'):
		infile = gzip.open(infile)

	query_match = []
	for rec in SeqIO.parse(infile, 'genbank'):
		for feature in rec.features:
			if query_feat == feature.type and \
			query_qualif in feature.qualifiers:
				list_to_search_against = feature.qualifiers[query_qualif]
				found = [s for s in list_to_search_against if query_term in s]
				if len(found) > 0:
					hit = ' '.join(found)
					query_match.append('>{} {} {}\n{}'.format(
						inbase, hit, rec.name, feature.extract(rec.seq)))

	if len(query_match) == 0:
		sys.stderr.write('ERROR: {} absent in {}\n'.format(
			query_term, infile))
		sys.exit(1)
	elif len(query_match) > 1:
		sys.stderr.write('WARNING: found >1 {} in {}\n'.format(
			query_term, infile))

	if opt.outfile:
		outfile = os.path.abspath(os.path.expanduser(opt.outfile))
		with open(opt.outfile, 'w') as ofh:
			for rec in query_match:
				ofh.write('{}\n'.format(rec))
	else:
		for rec in query_match:
			print rec

if __name__ == '__main__':
	main()