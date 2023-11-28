#!/usr/bin/env python


import gzip
import os
import sys
from argparse import ArgumentParser
from BCBio import GFF  #pip install bcbio-gff
from Bio import SeqIO
from Bio.SeqUtils import GC

def parseArgs():
	parser = ArgumentParser(description='Converts a GenBank file containing '
		'nucleotide sequences into a General Feature Format (GFF) file')
	parser.add_argument('-i', '--infile', required=True,
		help='input GenBank Format file <.gff||.gff3>')
	parser.add_argument('-o', '--outfile', required=False, default=None,
		help='output General Feature Format (.gff or .gff3) file [stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.realpath(os.path.expanduser(opt.infile))
	if opt.outfile is not None:
		ofh = os.path.realpath(os.path.expanduser(opt.outfile))
	else:
		ofh = sys.stdout

	records = []
	if ifh.endswith('.gz'):
		ifh = gzip.open(ifh)
	for rec in SeqIO.parse(ifh, 'genbank'):
		# halt if nucleotides are absent rather than printing default 'N'
		if float(GC(rec.seq)) == 0:
			sys.stderr.write('ERROR: {} appears to lack nucleotides\n'.\
				format(ifh))
			sys.exit(1)
		records.append(rec)
	with open(ofh, 'w') as o:
		GFF.write(records, o, include_fasta=True)

if __name__ == '__main__':
	main()
