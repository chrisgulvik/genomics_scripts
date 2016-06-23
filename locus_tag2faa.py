#!/usr/bin/env python


# Usage: python script.py <input.gbk> <locus_tag>

from Bio import SeqIO
import sys

infile = SeqIO.parse(sys.argv[1], 'genbank')
want = str(sys.argv[2])

for record in infile:
	for s in record.features:
		if 'locus_tag' and 'translation' and 'protein_id' and 'product' in s.qualifiers:
			locus_record = s.qualifiers['locus_tag'][0]
			if locus_record == want:
				print '>' + s.qualifiers['locus_tag'][0] + '|' + s.qualifiers['protein_id'][0] + '|' + s.qualifiers['product'][0] + '\n' + s.qualifiers['translation'][0]
