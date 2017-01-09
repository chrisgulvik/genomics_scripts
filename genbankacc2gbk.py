#!/usr/bin/env python


import os
import re
import sys
import time
from argparse import ArgumentParser
from Bio import Entrez

def parseArgs():
	parser = ArgumentParser(description='Given accession(s), retrieve and '
		'save a GenBank file.', add_help=False, epilog='If more than one '
		'accession is provided, all will be merged into a single output file')
	req = parser.add_argument_group('Required')
	req.add_argument('acc', metavar='accession', nargs='+',
		help='accession number(s) to retrieve')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-l', '--min-length', type=int, metavar='INT',
		default=7500, help='minimum record length (base pairs) to save '
		'without error warnings [7500]')
	opt.add_argument('-o', '--outfile', type=str, metavar='FILE',
		default=None, help='output GenBank filename [<acc>.gbk] ')
	return parser.parse_args()

def multi_gbk(acc_list, minlen, gbk):
	fetched = Entrez.efetch(db='nucleotide', id=acc_list,
		rettype='gbwithparts', retmode='text').read().rstrip('\n')
	print '\tfound {} GenBank records'.format(
		str(len(fetched.split('//\n\n'))))
	if len(fetched) < minlen:
		sys.exit('ERROR: suspiciously short GenBank record')
	contigs = []
	for record in fetched.split('\n\n'):
		if acc_list[0] not in record:
			sys.exit('ERROR: improperly retrieved accession {}'.format(
				acc_list[0]))
		del acc_list[0]
		contig = record.rstrip('\n').split('\n')
		i = contig.index('ORIGIN      ')
		contig.insert(i, 'CONTIG')  #adding this line to multigenbank style enables biopython to parse in SeqIO as 'gb'
		contigs.extend(contig)
	with open(gbk, 'w') as of:
		of.write('\n'.join(s for s in contigs))

def get_gbk(accession, minlen, gbk):
	fetched = Entrez.efetch(db='nucleotide', id=accession,
		rettype='gbwithparts', retmode='text').read()
	if accession not in fetched:
		sys.exit('ERROR: unable to retrieve accession number {}'.format(
			accession))
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
				acc_suffs = range(wgs_acc_suff_first, wgs_acc_suff_last + 1, 1)
				wgs_acc_list = [wgs_acc_pref_first + 
					str('{num:08d}'.format(num=s)) for s in acc_suffs]
				print ('\t{} accession is an incomplete chromosome.\n'
				'\tRetrieving {} individual contigs...'.format(
					accession, str(len(wgs_acc_list))))
				time.sleep(3)  #be nice to NCBI before doing a second efetch request
				multi_gbk(wgs_acc_list, minlen, gbk)
	elif len(fetched) < minlen:
		sys.exit('{}\nERROR: suspiciously short record for {}'.format(
			fetched, accession))
	else:
		with open(gbk, 'w') as of:
			of.write(fetched)
		print '\tfound {}'.format(accession)

def main():
	opts = parseArgs()
	acc    = opts.acc
	minlen = opts.min_length
	Entrez.email = 'yourname@univ.edu'

	# handle outfile naming
	if opts.outfile is not None:
		out = os.path.abspath(opts.outfile)
	elif len(acc) == 1:
		out = os.path.join(os.getcwd(), '{}.gbk'.format(acc[0]))
	elif len(acc) > 1:
		out = os.path.join(os.getcwd(), '{}.gbk'.format(','.join(acc)))
	
	# fetch genbank record(s) and save as single file
	if len(acc) == 1:
		get_gbk(acc[0], minlen, out)
	elif len(acc) > 1:
		multi_gbk(acc, minlen, out)

if __name__ == '__main__':
	main()
