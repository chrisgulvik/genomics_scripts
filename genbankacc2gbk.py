#!/usr/bin/env python


import os
import re
import sys
import time
from Bio import Entrez

def multi_gbk(acc_list, gbk):
	fetched = Entrez.efetch(db='nucleotide', id=acc_list, rettype='gbwithparts', retmode='text').read().rstrip('\n')
	print '\tfound {} GenBank records'.format(str(len(fetched.split('//\n\n'))))
	if len(fetched) < 7500:
		sys.exit('ERROR: suspiciously short GenBank record')
	contigs = []
	for record in fetched.split('\n\n'):
		if acc_list[0] not in record:
			sys.exit('ERROR: improperly retrieved accession {}'.format(acc_list[0]))
		del acc_list[0]
		contig = record.rstrip('\n').split('\n')
		i = contig.index('ORIGIN      ')
		contig.insert(i, 'CONTIG')  #adding this line to multigenbank style enables biopython to parse in seqio as gb
		contigs.extend(contig)
	with open(gbk, 'w') as of:
		of.write('\n'.join(s for s in contigs))

def get_gbk(accession, gbk):
	fetched = Entrez.efetch(db='nucleotide', id=accession, rettype='gbwithparts', retmode='text').read()
	if accession not in fetched:
		sys.exit('ERROR: unable to retrieve accession number {}'.format(accession))
	if 'WGS         ' in fetched: #entry is the master record for a whole genome shotgun sequencing project and contains no sequence data
		for line in fetched.split('\n'):
			if line.startswith('WGS         '):
				wgs_accs = line.lstrip('WGS         ').split('-')
				wgs_acc_pref_first = re.sub('\d', '', wgs_accs[0])
				wgs_acc_pref_last  = re.sub('\d', '', wgs_accs[1])
				if wgs_acc_pref_first != wgs_acc_pref_last:
					sys.exit('ERROR: unable to parse WGS info')
				wgs_acc_suff_first = int(re.sub('\D', '', wgs_accs[0]))
				wgs_acc_suff_last  = int(re.sub('\D', '', wgs_accs[1]))
				acc_suffs = range(wgs_acc_suff_first, wgs_acc_suff_last, 1)
				wgs_acc_list = [wgs_acc_pref_first + str('{num:08d}'.format(num=s)) for s in acc_suffs]
				print '\t{} accession is an incomplete chromosome.\n\tRetrieving {} individual contigs...'.format(accession, str(len(wgs_acc_list)))
				time.sleep(3)  #be nice to NCBI before doing a second efetch request
				multi_gbk(wgs_acc_list, gbk)
	elif len(fetched) < 7500:
		sys.exit('{}\nERROR: suspiciously short record for {}'.format(fetched, accession))
	else:
		with open(gbk, 'w') as of:
			of.write(fetched)
		print '\tfound {}'.format(accession)

def main():
	if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
		sys.exit('\tUsage: genbankacc2gbk.py <acc> [additional GenBank accessions...]\n\n\tif more than one accession is provided all will be merged into a single output file\n')
	Entrez.email = 'yourname@univ.edu'
	if len(sys.argv) == 2:
		get_gbk(sys.argv[1], '{}.gbk'.format(sys.argv[1]))
	elif len(sys.argv) > 2:
		accs = [acc for acc in sys.argv[1:]]
		multi_gbk(accs, '{}.gbk'.format(','.join(accs)))

if __name__ == '__main__':
	main()
