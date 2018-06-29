#!/usr/bin/env python


import argparse
import gzip
import os
import re
import sys
from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqUtils import GC


def parseArgs():
	parser = argparse.ArgumentParser(description='filters contigs (or '
		'scaffolds) based on length, coverage, GC skew, and complexity',
		epilog='NOTE: Headers in IDBA contigs must have spaces removed or '
		'replaced with a non-whitespace character.', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', metavar='FILE',
		required=True, help='input FastA file from IDBA, SKESA, SPAdes, or '
		'Velvet, optionally gunzip compressed')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-b', '--baseheader', default=None, metavar='STR',
		help='contig header prefix (with _<COUNTER> suffix) '
		'[basename infile .ext]')
	opt.add_argument('-c', '--cov', type=int, default=5, metavar='INT',
		help='minimum coverage (for SKESA, SPAdes, and Velvet) or minimum '
		'read count (for IDBA); set to 0 if deflines lack coverage [5]')
	opt.add_argument('-g', '--gcskew', default=True, action='store_false',
		help='keep >88 and <12%% GC contigs')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-l', '--len', type=int, default=500, metavar='INT',
		help='minimum contig length (in bp) [500]')
	opt.add_argument('-m', '--complex', default=True, action='store_false',
		help='keep low-complexity contigs with only 1 or 2 nucleotides '
		'represented')
	opt.add_argument('-o', '--outfile', default=None, metavar='FILE',
		help='output FastA file [./basename.filt.infile_ext]')
	opt.add_argument('--defline-fmt', choices={'rename', 'rename_keepOrig',
		'no_change'}, default='rename_keepOrig', help='rename: fully rename '
		'with -b prefix and _<COUNTER> suffix; '
		'rename_keepOrig: includes rename style but also appends a '
		'|OrigDefln=<InputDefln> suffix; no_change: retain exactly as input '
		'[rename_keepOrig]')
	opt.add_argument('--no-sort', default=False, action='store_true',
		help='skip sorting contigs by descending length')
	opt.add_argument('--plasmid-spades', default=False, action='store_true',
		help='split contigs containing component_<INT> headers into separate '
		'files according to their respective replicon')
	return parser.parse_args()

def filter_contig(record, min_len, min_cov, gc, complexity):
	''' removes short and low coverage contigs '''
	if min_cov == 0:
		if len(record.seq) >= min_len:
			accepted_record = gc_filter(record, gc, complexity)
			return accepted_record

	if len(record.seq) >= min_len:
		# SPAdes v3.12.0 and Velvet v1.2.10
		cov_velvet = re.compile('cov_([0-9]{1,}.[0-9]{1,})')
		cov_match = cov_velvet.search(record.name)
		if cov_match:
			if float((cov_match.group(0)).lstrip('cov_')) >= min_cov:
				accepted_record = gc_filter(record, gc, complexity)
				if accepted_record:
					return accepted_record
		else:
			# SKESA v2.2 appends _Circ after cov if suggested circular contig
			cov_skesa = re.compile('Contig_[0-9]{1,}_[0-9]{1,}.[0-9]{1,}')
			cov_match = cov_skesa.search(record.name)
			if cov_match:
				if float((cov_match.group(0)).split('_')[-1]) >= min_cov:
					accepted_record = gc_filter(record, gc, complexity)
					if accepted_record:
						return accepted_record
			else:
				# IDBA (requires discarding or replacing spaces in
				# deflines due to Biopython SeqIO's parsing)
				cov_idba = re.compile('read_count_([0-9]+)')
				read_match = cov_idba.search(record.name)
				if read_match:
					if int((read_match.group(0)).lstrip('read_count_')) >= \
					min_cov:
						accepted_record = gc_filter(record, gc, complexity)
						if accepted_record:
							return accepted_record
				else:
					# keep all contigs if deflines renamed and lack depth data
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

def rename_seqrec_defline(record, baseheader, i, rename_type):
	''' renames sequence record object according to specified type '''
	if rename_type == 'rename_keepOrig':
		return SeqRecord.SeqRecord(id='{}_{}|'
			'OrigDefln={}'.format(baseheader, i, record.id),
			seq=record.seq, description='')
	elif rename_type == 'rename':
		return SeqRecord.SeqRecord(id='{}_{}'.format(
			baseheader, i), seq=record.seq, description='')
	elif rename_type == 'no_change':
		return SeqRecord.SeqRecord(id=record.id, seq=record.seq,
			description='')
	else:
		sys.stderr.write('ERROR: {} is an unrecognized option\n'.format(
			rename_type))
		sys.exit(1)

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
		in_ext  = os.path.splitext(infile.strip('.gz'))[-1]
		out_ext = '.filt' + in_ext
		base    = os.path.splitext(os.path.basename(infile.strip('.gz')))[0]
		outfile = base + out_ext
	if baseheader is None:
		baseheader = os.path.splitext(os.path.basename(
			infile.strip('.gz')))[0]

	if args.plasmid_spades:
		out_extension = os.path.splitext(outfile)[-1]
		out_base      = os.path.splitext(outfile)[0]
	else:
		out = outfile

	if infile.endswith('.gz'):
		infile = gzip.open(infile)
	records = []
	i = 1
	for record in SeqIO.parse(infile, 'fasta'):
		class_name = type(record)  #'Bio.SeqRecord.SeqRecord'
		r = filter_contig(record, min_len, min_cov, gc, complexity)
		if isinstance(r, class_name):
			if args.plasmid_spades:
				if 'component_' in r.id:
					component_regex = re.compile('component_([0-9]+)')
					# change from 0- to 1-based counting
					component_nr = str(int(component_regex.search(r.id).\
						group(0).lstrip('component_')) + 1)
					out = out_base + '.' + component_nr + out_extension
					with open(out, 'a') as filtered_fasta:
						SeqIO.write(rename_seqrec_defline(r, baseheader, i,
							args.defline_fmt), filtered_fasta, 'fasta')
					i += 1
				else:
					sys.stderr.write('ERROR: --plasmid-spades option invoked '
						'but component_ absent in {} record\n'.format(r.id))
					sys.exit(1)
			else:
				records.append((r, len(r.seq)))

	if len(records) > 0:
		if not args.no_sort:
			records = sorted(records, key=lambda x: x[1], reverse=True)
		out_recs = []
		i = 1
		for r, _ in records:
			out_recs.append(rename_seqrec_defline(r, baseheader, i, 
				args.defline_fmt))
			i += 1
		SeqIO.write(out_recs, out, 'fasta')

if __name__ == '__main__':
	main()
