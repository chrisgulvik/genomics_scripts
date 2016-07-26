#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from BCBio import GFF  #pip install bcbio-gff
from BCBio.GFF import GFFParser
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
	parser = ArgumentParser(description='Extracts gene sequences from a General Feature Format (GFF) file given a locus_tag or comma-separated list of locus_tags',
							epilog='NOTE: a BED format of extracted gene(s) is printed to standard out')
	parser.add_argument('-i', '--infile', required=True, help='input General Feature Format file <.gff||.gff3>')
	parser.add_argument('-l', '--locus_tag', required=True, default=None, help='locus tag(s) to extract; comma-separate if >1')
	parser.add_argument('-o', '--outpath', required=False, default=None, help='output directory [cwd]')
	args = parser.parse_args()
	return args

def extract_seq(start, end, strand, seq, locus_name):
	''' returns a SeqRecord obj given a nucleotide sequence and record
	info to extract specific sites including: start, stop, and strand '''
	extraction = seq[int(start):int(end)]
	if strand < 0:
		return SeqRecord(Seq(extraction, generic_dna).reverse_complement(),
						id=locus_name, description='')
	else:
		return SeqRecord(Seq(extraction, generic_dna),
						id=locus_name, description='')

def locus_tag2record_info(want, infile):
	''' returns a tuple containing record info 
	given a GFF file and locus_tag to locate '''
	parser = GFFParser()
	with open(infile, 'r') as gff:
		for record in parser.parse(gff, limit_info = dict(gff_type = ['gene', 'CDS', 'locus_tag', 'ID', 'product'])):
			for feature in record.features:
				if feature.type == 'gene' and 'locus_tag' in feature.qualifiers:
					locus_tag = feature.qualifiers.get('locus_tag', None)[0]
					if locus_tag == want:
						return (feature.location.start.position, feature.location.end.position, feature.strand, record.id, record.seq.tostring())

def main():
	opts = parseArgs()
	infile = opts.infile
	if ',' in opts.locus_tag:
		locus_want = opts.locus_tag.split(',')
	else:
		locus_want = opts.locus_tag
	if opts.outpath is not None:
		outpath = os.path.abspath(opts.outpath)
		if not os.path.exists(outpath):
			os.mkdir(outpath)
	else:
		outpath = os.getcwd()

	# Extract one or more locus tags
	if isinstance(locus_want, list):
		seq_rec = []
		for locus in locus_want:
			(start, end, strand, contig_id, contig_seq) = locus_tag2record_info(locus, infile)
			if strand < 0:
				strand_symbol = '-'
			else:
				strand_symbol = '+'
			print '{}\t{}\t{}\t{}\t0\t{}'.format(contig_id, start, end, locus, strand_symbol)
			sr = extract_seq(start, end, strand, contig_seq, locus)
			seq_rec.append(sr)
	else:
		(start, end, strand, contig_id, contig_seq) = locus_tag2record_info(locus_want, infile)
		if strand < 0:
			strand_symbol = '-'
		else:
			strand_symbol = '+'
		print '{}\t{}\t{}\t{}\t0\t{}'.format(contig_id, start, end, locus_want, strand_symbol)
		seq_rec = extract_seq(start, end, strand, contig_seq, locus_want)

	# Write one or more extracted genes to FastA
	SeqIO.write(seq_rec, os.path.join(outpath,'extracted.ffn'), 'fasta')

if __name__ == '__main__':
	main()
