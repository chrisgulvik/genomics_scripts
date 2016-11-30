#!/usr/bin/env python


# Example usage: python <script>.py -i input.gbk -f rpoB -t gene -o rpoB.sample.fa

import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
	parser = ArgumentParser(description='Extracts a nucleotide sequence from '
		'a GenBank file', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input GenBank file with translations')
	req.add_argument('-f', '--find', required=True, metavar='STR',
		help='string to find')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-t', '--type', default='locus_tag', metavar='STR',
		help='genbank feature type to search in, e.g., CDS, gene, rRNA, '
		'source, misc_feature')
	opt.add_argument('-o', '--outfile', default='parsed.fas', metavar='FILE',
		help='output FastA file  [./parsed.fas]')
	return parser.parse_args()

def index_genbank(rec, feature_type, qualifier) :
	gb_idx = {}
	for (i, feat) in enumerate(rec.features):
		if feat.type == feature_type:
			if qualifier in feat.qualifiers:
				for s in feat.qualifiers[qualifier]:
					if s in gb_idx:
						sys.exit('ERROR: duplicate {}'.format(s))
					else:
						gb_idx[s] = i
	return gb_idx

def locus_tag2gene_seq(gbk, ft_type, query, outfile):
	r = 0
	for rec in SeqIO.parse(open(gbk, 'r'), 'genbank'):
		gbk_idx = index_genbank(rec, 'CDS', ft_type)
		if query in gbk_idx:
			gene_seq = str(rec.features[gbk_idx[query]].extract(rec.seq))
			product = str(rec.features[gbk_idx[query]].qualifiers['product'][0])
			record = SeqRecord(Seq(gene_seq, generic_dna), id=query,
				description=product, name='')
			r += 1
	if r == 1:
		SeqIO.write(record, outfile, 'fasta')
	elif r > 1:
		print 'ERROR: found >1 {}'.format(query)
	else:
		print 'ERROR: {} absent'.format(query)

def main():
	opts = parseArgs()
	locus_tag2gene_seq(opts.infile, str(opts.type), str(opts.find), opts.outfile)

if __name__ == '__main__':
	main()
