#!/usr/bin/env python


# Example: <script>.py -i input.gbk -s VFDB_setA_pro -f CDS -q inference -r product
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Parses a GenBank file for data '
		'given a search term and specified fields to query. An additional '
		'parsed field is also reported in tab-delimited format.', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input GenBank file')
	req.add_argument('-s', '--query-search', required=True, metavar='STR',
		type=str, help='search term to look for within the query feature '
		'and query qualifier')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-f', '--query-feature', default='CDS', metavar='STR',
		help='genbank feature type to search in, e.g., CDS, gene, rRNA, '
		'source, tRNA, misc_feature')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', default='parsed.tab', metavar='FILE',
		help='output tab-delimited file containing <locus_tag>\\t'
		'<query-qualifier>\\t<report-qualifier>  [./parsed.tab]')
	opt.add_argument('-q', '--query-qualifier', default='inference',
		metavar='STR', help='qualifier term within each genbank feature to '
		'search in, e.g., locus_tag, inference, codon_start, product, '
		'transl_table, translation')
	opt.add_argument('-r', '--report-qualifier', default='product',
		metavar='STR', help='additional qualifier term to parse data from '
		'and report when queries are found')
	return parser.parse_args()

def main():
	opts = parseArgs()
	gbk          = opts.infile
	query_feat   = opts.query_feature
	query_qualif = opts.query_qualifier
	query_str    = opts.query_search
	rep_qualif   = opts.report_qualifier
	outfile      = opts.outfile

	l = []
	for rec in SeqIO.parse(open(gbk, 'r'), 'genbank'):
		for feature in rec.features:
			if query_feat == feature.type and \
			query_qualif in feature.qualifiers:
				list_to_search_against = feature.qualifiers[query_qualif]
				found = [s for s in list_to_search_against if query_str in s]
				if len(found) > 0:
					hit = ' '.join(found)
					l.append('{}\t{}\t{}'.format(
						feature.qualifiers['locus_tag'][0], hit,
						''.join(feature.qualifiers[rep_qualif][0])))

	with open(opts.outfile, 'w') as ofh:
		for i in l:
			ofh.write('{}\n'.format(i))

if __name__ == '__main__':
	main()
