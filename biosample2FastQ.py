#!/usr/bin/env python

import os
import sys

# Downloads raw FastQ reads given a BioSample accession (SAMN########)
# or SRA Run sccession (SRR#######)
# Depends on:
#     entrez direct scripts <http://www.ncbi.nlm.nih.gov/books/NBK179288/>
#     a SRA toolkit binary  <http://www.ncbi.nlm.nih.gov/books/NBK158899/>

os.system('esearch -db sra -query {} | efetch -format docsum | '
		  'grep \'<Runs>\' > /tmp/run_info'.format(sys.argv[1]))

download = []
with open('/tmp/run_info', 'r') as run_info:
	for l in run_info:
		if 'SRR' in l:
			SRR = l.lstrip('<Run acc="').split('"')[0]
			download.append(SRR)

print 'downloading {}...'.format(','.join(SRR))

for SRR in download:
	os.system('fastq-dump -O {} --dumpbase --split-files --readids -Q 33 '
			  '-defline-qual \'+\' --defline-seq \'@$ac_$sn[_$rn]/$ri\' '
			  '{}'.format(os.path.join(os.getcwd(), sys.argv[1]), SRR))
	print '\tdownloaded {}'.format(SRR)
