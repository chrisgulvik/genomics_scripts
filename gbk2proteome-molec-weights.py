#!/usr/bin/env python

# Given a genbank file, prints to stdout a tab-delimited summary
# of the proteome including the molecular weight in Daltons, where
# '~<value>' indicates the MW is estimated due to ambiguous amino acid(s)

# Usage: python script.py <input.gbk>

from Bio import SeqIO
from Bio import SeqUtils
from decimal import Decimal, ROUND_UP
import os
import string
import sys

def calc_molec_weight(protein_seq):
	return Decimal(SeqUtils.molecular_weight(protein_seq, seq_type='protein')
		).quantize(Decimal('.001'), rounding=ROUND_UP)

infile   = SeqIO.parse(sys.argv[1], 'genbank')
filename = os.path.splitext(os.path.splitext(os.path.basename(sys.argv[1]))[0])[0]

for rec in infile:
	for s in rec.features:
		if (s.type == 'CDS' and 'gene') in s.qualifiers:
			gene    = s.qualifiers['gene'][0]
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
			print '{}\t{}\t{}\t{}'.format(gene,locus,product,protein_mw)
