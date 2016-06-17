#!/usr/bin/env python

# Given a genbank file, prints to stdout a tab-delimited summary
# of the proteome including the molecular weight in Daltons
# Usage: python script.py <input.gbk>

from Bio import SeqIO
from Bio import SeqUtils
from decimal import Decimal, ROUND_DOWN, ROUND_UP
import os
import sys

infile = SeqIO.parse(sys.argv[1], 'genbank')
filename = os.path.splitext(os.path.splitext(os.path.basename(sys.argv[1]))[0])[0]

for rec in infile:
	for s in rec.features:
		if (s.type == 'CDS' and 'gene') in s.qualifiers:
			gene       = s.qualifiers['gene'][0]
			locus      = s.qualifiers['locus_tag'][0]
			product    = s.qualifiers['product'][0]
			protein_mw = Decimal(SeqUtils.molecular_weight(s.qualifiers['translation'][0],
							seq_type='protein')).quantize(Decimal('.001'), rounding=ROUND_UP)
			print '{}\t{}\t{}\t{}'.format(gene,locus,product,protein_mw)
