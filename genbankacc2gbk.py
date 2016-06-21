#!/usr/bin/env python


import os
import re
import shutil
import sys
from Bio import Entrez

def multi_gbk(acc_list, gbk):
	# time.sleep(3)  #be nice to NCBI
	# fetched = Entrez.efetch(db='nucleotide', id=acc_list, rettype='gbwithparts', retmode='text').read()
	# for record in fetched.split('//\n'):
	# 	if wgs_acc_list[0] not in record:
	# 		sys.exit('ERROR: improperly retrieved accession {}'.format(wgs_acc_list[0]))
	# 	del wgs_acc_list[0]
	# with open(gbk, 'w') as of:
	# 	of.write(fetched)
	for acc in acc_list:
		url = 'http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={}&rettype=gbwithparts&retmode=gb'.format(acc)
		os.system('wget --no-glob -O {} {}'.format(acc + '.gbk.tmp', url))
		if acc not in open(acc + '.gbk.tmp', 'r').read():
			sys.exit('ERROR: unable to find {} in retrieved gbk.tmp file'.format(acc))
	with open(gbk, 'w') as outfile:
		for f in [s + '.gbk.tmp' for s in acc_list]:
			with open(f, 'r') as contig:
				shutil.copyfileobj(contig, outfile)
				os.remove(f)
	if len(fetched) < 7500:
		sys.exit('ERROR: suspiciously short GenBank record')

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
				wgs_acc_list = [wgs_acc_pref_first + str(s) for s in acc_suffs]
				print '\t{} accession is an incomplete chromosome.\n\tRetrieving {} individual contigs...'.format(accession, str(len(wgs_acc_list)))
				multi_gbk(wgs_acc_list, gbk)
	elif len(fetched) < 7500:
		sys.exit('{}\nERROR: suspiciously short record for {}'.format(fetched, accession))
	else:
		with open(gbk, 'w') as of:
			of.write(fetched)
		print '\tfound {}'.format(accession)

def main():
	Entrez.email = 'yourname@univ.edu'
	if len(sys.argv) == 2:
		get_gbk(sys.argv[1], '{}.gbk'.format(sys.argv[1]))
	elif len(sys.argv) > 2:
		accs = [acc for acc in sys.argv[1:]]
		multi_gbk(wgs_acc_list, '{}.gbk'.format(','.join(accs)))
	else:
		sys.exit('Usage: genbankacc2gbk.py <acc> [additional GenBank accessions...]')

if __name__ == '__main__':
	main()