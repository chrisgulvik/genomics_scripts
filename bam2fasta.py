#!/usr/bin/env python


import argparse
import os
import sys

import pysam
from Bio import SeqIO, Seq, SeqRecord

def parseArgs():
	parser = argparse.ArgumentParser(description='extracts read sequences '
		'aligned to a reference from a BAM, CRAM, or SAM file and outputs '
		'the reads in FastA format', add_help=False, epilog='alternative to '
		'Picard\'s SamToFastq and HudsonAlpha\'s bam2fastq')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', metavar='FILE',
		required=True, help='input BAM or SAM file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', default=None, metavar='FILE',
		help='output FastA file [stdout]')
	return parser.parse_args()

def bam_to_seqrec(infile):
	''' converts SAM or BAM file to Biopython.SeqRecord object '''
	if infile.endswith('.bam')
		alns = pysam.AlignmentFile(infile, 'rb')
	elif infile.endswith('.sam')
		alns = pysam.AlignmentFile(infile, 'r')
	elif infile.endswith('.cram')
		alns = pysam.AlignmentFile(infile, 'rc')
	else:
		sys.stderr.write('ERROR: unable to detect whether input file is a BAM '
			'or SAM format based on filename suffix\n')
		sys.exit(1)
	for read in alns:
		seq = Seq.Seq(read.seq)
		if read.is_unmapped:
			continue
		if read.is_reverse: #read aligned to the negative strand of reference
			seq = seq.reverse_complement()
		rec = SeqRecord.SeqRecord(seq, read.qname, '', '')
		yield rec

def main(infile):
	args   = parseArgs()
	infile = os.path.abspath(os.path.expanduser(args.infile))
	if args.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(args.outfile)), 'w')
	else:
		ofh = sys.stdout

	SeqIO.write(bam_to_seqrec(infile), ofh, 'fasta')

if __name__ == '__main__':
	main()
