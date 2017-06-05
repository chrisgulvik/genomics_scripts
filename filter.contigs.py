#!/usr/bin/env python


import argparse
import os
import re
import sys
from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqUtils import GC


def parseArgs():
	parser = argparse.ArgumentParser(description='filters contigs (or scaffolds) based on length, coverage, GC skew, and complexity',
		epilog='NOTE: Headers in IDBA contigs must have spaces removed or replaced with a non-whitespace character.', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', metavar='FILE',
		required=True, help='input FastA file from SPAdes or Velvet')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-b', '--baseheader', metavar='STR',
		help='contig header prefix (with \'_#\' as suffix) [basename infile]')
	opt.add_argument('-c', '--cov', type=int, default=5, metavar='INT',
		help='minimum coverage (for SPAdes and Velvet) or minimum read count '
		'(for IDBA) [5]; set to 0 if deflines lack coverage')
	opt.add_argument('-g', '--gcskew', default=True, action='store_false',
		help='switch to turn on saving >88 and <12%% GC contigs')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-l', '--len', type=int, default=500, metavar='INT',
		help='minimum contig length (in bp) [500]')
	opt.add_argument('-m', '--complex', default=True, action='store_false',
		help='switch to turn on saving low-complexity contigs')
	opt.add_argument('-o', '--outfile', metavar='FILE',
		help='output FastA file [./basename.filtered.infile_ext]')
	opt.add_argument('-ps', '--plasmid-spades', default=False, action='store_true',
		help='switch to split up contigs containing component_<int> headers '
		'into separate files; zero-based <int> in input deflines are '
		'converted to one-based in output; useful for plasmidSPAdes parsing')
	return parser.parse_args()

def filter_contig(record, min_len, min_cov, gc, complexity):
	''' removes short and low coverage contigs '''
	if min_cov == 0:
		if len(record.seq) >= min_len:
			accepted_record = gc_filter(record, gc, complexity)
			return accepted_record

	if len(record.seq) >= min_len:
		cov_pattern = re.compile('cov_([0-9.]+)')  #SPAdes v3.10.1 and Velvet v1.2.10 lack trailing underscore
		cov_match = cov_pattern.search(record.name)
		if cov_match:
			if float((cov_match.group(0)).lstrip('cov_')) >= min_cov:
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
	infile     = args.infile
	outfile    = args.outfile
	baseheader = args.baseheader
	min_cov    = args.cov
	min_len    = args.len
	gc         = args.gcskew
	complexity = args.complex

	if outfile is None:
		in_ext  = os.path.splitext(infile)[-1]
		out_ext = '.filtered' + in_ext
		base    = os.path.splitext(os.path.basename(infile))[0]
		outfile = base + out_ext
	if baseheader is None:
		baseheader = os.path.splitext(os.path.basename(infile))[0]

	if args.plasmid_spades:
		out_extension = os.path.splitext(outfile)[-1]
		out_base      = os.path.splitext(outfile)[0]
	else:
		out = outfile

	with open(infile, 'rU') as input_fasta:
		i = 1
		for record in SeqIO.parse(input_fasta, 'fasta'):
			class_name = type(record)  #'Bio.SeqRecord.SeqRecord'
			r = filter_contig(record, min_len, min_cov, gc, complexity)
			if isinstance(r, class_name):
				if args.plasmid_spades:
					if 'component_' in r.id:
						component_regex = re.compile('component_([0-9]+)')
						component_nr = str(int(component_regex.search(r.id).group(0).lstrip('component_')) + 1)
						out = out_base + '.' + component_nr + out_extension
						with open(out, 'a') as filtered_fasta:
							renamed_rec = SeqRecord.SeqRecord(id='{}_{}'.format(baseheader, i), seq=r.seq, description='')
							SeqIO.write(renamed_rec, filtered_fasta, 'fasta')
					else:
						sys.exit('ERROR: --plasmid-spades option invoked but component_ absent in {}'.format(r.id))
				else:
					with open(out, 'a') as filtered_fasta:
						renamed_rec = SeqRecord.SeqRecord(id='{}_{}'.format(baseheader, i), seq=r.seq, description='')
						SeqIO.write(renamed_rec, filtered_fasta, 'fasta')
				i += 1

if __name__ == '__main__':
	main()
