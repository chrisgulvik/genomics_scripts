#!/usr/bin/env python


import gzip
import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
	parser = ArgumentParser(description='Extracts a nucleotide sequence from '
		'a GenBank file', epilog='NOTE: if >1 record is found, a warning is '
		'printed to stderr, and all matches are sent as output', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input GenBank file, optionally gunzip compressed')
	req.add_argument('-f', '--find', required=True, metavar='STR',
		help='string to find')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-t', '--type', default='locus_tag', metavar='STR',
		help='genbank feature type to search in, e.g., CDS, gene, rRNA, '
		'source, misc_feature [locus_tag]')
	opt.add_argument('-o', '--outfile', required=False, default=None,
		help='FastA output [stdout]')
	return parser.parse_args()

def index_genbank(rec, feature_type, qualifier) :
	gb_idx = {}
	for (i, feat) in enumerate(rec.features):
		if feat.type == feature_type:
			if qualifier in feat.qualifiers:
				for s in feat.qualifiers[qualifier]:
					if s in gb_idx:
						sys.stderr.write('ERROR: >1 {} in input\n'.format(s))
						sys.exit(1)
					else:
						gb_idx[s] = i
	return gb_idx

def locus_tag2gene_seq(gbk, ft_type, query, outfile):
	if gbk.endswith('.gz'):
		gbk = gzip.open(gbk)
	records = []
	for rec in SeqIO.parse(gbk, 'genbank'):
		gbk_idx = index_genbank(rec, 'CDS', ft_type)
		if query in gbk_idx:
			gene_seq = str(rec.features[gbk_idx[query]].extract(rec.seq))
			product = str(rec.features[gbk_idx[query]].qualifiers['product'][0])
			record = SeqRecord(Seq(gene_seq, generic_dna), id=query,
				description=product, name='')
			records.append(record)

	if len(records) == 0:
		sys.stderr.write('ERROR: {} absent in {}\n'.format(query, gbk))
		sys.exit(1)
	elif len(records) > 1:
		sys.stderr.write('WARNING: found >1 {} in {}\n'.format(query, gbk))
	
	if outfile:
		outfile = os.path.abspath(os.path.expanduser(outfile))
		SeqIO.write(records, outfile, 'fasta')
	else:
		SeqIO.write(records, sys.stdout, 'fasta')

def main():
	opts = parseArgs()
	locus_tag2gene_seq(opts.infile, str(opts.type), str(opts.find), opts.outfile)

if __name__ == '__main__':
	main()
