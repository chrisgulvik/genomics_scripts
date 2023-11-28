#!/usr/bin/env python

# Given a genbank file, prints to stdout a tab-delimited summary
# of the proteome including the molecular weight in Daltons, where
# '~<value>' indicates the MW is estimated due to ambiguous amino acid(s)

# Usage: python script.py <input.gbk>

import os
import string
import sys
from argparse import ArgumentParser
from Bio import SeqIO, SeqUtils
from decimal import Decimal, ROUND_UP

def parseArgs():
	parser = ArgumentParser(description='Calculates the molecular weight of proteins in a GenBank file')
	parser.add_argument('-i', '--infile', required=True, help='input GenBank file with translations')
	args = parser.parse_args()
	return args

def calc_molec_weight(protein_seq):
	return Decimal(SeqUtils.molecular_weight(protein_seq, seq_type='protein')
		).quantize(Decimal('.001'), rounding=ROUND_UP)

def main():
	opts = parseArgs()
	infile = opts.infile
	gb = SeqIO.parse(infile, 'genbank')

	print 'Locus_Tag\tProduct\tProtein_MolecWt_[Da]\tProtein_Seq'  #Header
	for rec in gb:
		for s in rec.features:
			if (s.type == 'CDS' and 'translation') in s.qualifiers:
				locus   = s.qualifiers['locus_tag'][0]
				product = s.qualifiers['product'][0]
				protein = s.qualifiers['translation'][0]
				if 'X' in protein:
					num_X = protein.count('X')
					estimated_mw = Decimal(num_X) * Decimal('128.16')  #UNK=(C5)(H6)(N1)(O3)=128.16
					peptide_mw   = calc_molec_weight(protein.replace('X', ''))
					protein_mw   = '~' + str(peptide_mw + estimated_mw)
				else:
					protein_mw   = calc_molec_weight(protein)
				print '{}\t{}\t{}\t{}'.format(locus,product,protein_mw,protein)

if __name__ == '__main__':
	main()
