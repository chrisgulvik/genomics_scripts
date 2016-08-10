#!/usr/bin/env python

import argparse
import os
import sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
	parser = argparse.ArgumentParser(
		description='splits multi-FastA file into individual FastA files')
	parser.add_argument('-i', '--infile',
		required=True, help='input multi-FastA file')
	parser.add_argument('-o', '--outdir',
		help='output directory path [cwd of input file]')
	parser.add_argument('-g', '--nogaps', action='store_true', default=False,
		help='toggle on removal of \'~\' and \'-\' gaps [off]')
	return parser.parse_args()

def main():
	args = parseArgs()
	infile = args.infile
	outdir = args.outdir

	if outdir is None:
		outdir = os.path.dirname(infile)

	mfasta = SeqIO.parse(infile, 'fasta')

	for rec in mfasta:
		seq  = str(rec.seq).upper()
		if args.nogaps:
			seq = seq.replace('-', '').replace('~', '')
		name = (rec.id).split()[0]
		new_rec = SeqRecord(Seq(seq, generic_dna), id=name, description='')
		SeqIO.write(new_rec, '{}.fasta'.format(name), 'fasta')


if __name__ == '__main__':
	main()
