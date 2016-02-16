#!/usr/bin/env python


import argparse
import os
import re
import sys
from Bio import SeqIO
from Bio.SeqUtils import GC


def parseArgs():
	parser = argparse.ArgumentParser(description='filters contigs (or scaffolds) based on length, coverage, GC skew, and complexity',
		epilog='NOTE: headers in IDBA contigs must have spaces removed or replaced with a non-whitespace character')
	parser.add_argument('-i', '--infile',
		required=True, help='input FastA file from SPAdes or Velvet')
	parser.add_argument('-o', '--outfile', help='output FastA file [basename.filtered.infile_ext]')
	parser.add_argument('-p', '--outpath', help='output path [cwd of input file]')
	parser.add_argument('-c', '--cov', type=int, default=5,
		help='minimum coverage (for SPAdes and Velvet) or minimum read count (for IDBA) [5]')
	parser.add_argument('-g', '--gcskew', default=True, action='store_false',
		help='switch to turn on saving >88 and <12%% GC contigs')
	parser.add_argument('-m', '--complex', default=True, action='store_false',
		help='switch to turn on saving low-complexity contigs')
	parser.add_argument('-l', '--len', type=int, default=500,
		help='minimum contig length (in bp) [500]')
	return parser.parse_args()

def filter_contig(record, min_len, min_cov, gc, complexity):
	''' removes short and low coverage contigs '''
	if len(record.seq) >= min_len:
		cov_pattern = re.compile('cov_([0-9.]+)_')  #SPAdes and Velvet
		cov_match = cov_pattern.search(record.name)
		if cov_match:
			if float((cov_match.group(0)).lstrip('cov_').rstrip('_')) >= min_cov:
				accepted_record = gc_filter(record, gc, complexity)
				if accepted_record:
					return accepted_record
		else:
			read_counts = re.compile('read_count_([0-9]+)') #IDBA (requires discarding or replacing spaces in deflines due to Biopython SeqIO's parsing)
			read_match = read_counts.search(record.name)
			if read_match:
				if int((read_match.group(0)).lstrip('read_count_')) >= min_cov:
					accepted_record = gc_filter(record, gc, complexity)
					if accepted_record:
						return accepted_record
			else:  #keep all contigs if deflines renamed and lack depth data
				accepted_record = gc_filter(record, gc, complexity)
				if accepted_record:
					return accepted_record

def gc_filter(record, gc, complexity):
	''' removes contigs with highly skewed GC '''
	if gc:
		gc_content = float(GC(record.seq))
		if 12 <= gc_content <= 88:
			accepted_record = complexity_filter(record, complexity)
			if accepted_record:
				return accepted_record
	else:
		accepted_record = complexity_filter(record, complexity)
		if accepted_record:
			return accepted_record

def complexity_filter(record, complexity):
	''' removes contigs with only 1 or 2 nucleotides represented '''
	if complexity:
		i = 0
		for nucleotide in ['A', 'T', 'C', 'G']:
			if nucleotide in record.seq:
				i += 1
		if i > 2:
			return record
	else:
		return record

def main():
	args = parseArgs()
	infile = args.infile
	outfile = args.outfile
	outpath = args.outpath
	min_cov = args.cov
	min_len = args.len
	gc = args.gcskew
	complexity = args.complex

	in_ext = os.path.splitext(os.path.basename(infile))[-1]
	if outfile is None:
		out_ext = '.filtered' + in_ext
		base = os.path.splitext(os.path.basename(infile))[0]
		outfile = base + out_ext
	if outpath is None:
		outpath = os.path.dirname(infile)

	with open(outfile, 'w') as filtered_fasta:
		with open(infile, 'rU') as input_fasta:
			for record in SeqIO.parse(input_fasta, 'fasta'):
				class_name = type(record)  #'Bio.SeqRecord.SeqRecord'
				r = filter_contig(record, min_len, min_cov, gc, complexity)
				if isinstance(r, class_name):
					SeqIO.write(r, filtered_fasta, 'fasta')

if __name__ == '__main__':
	main()
