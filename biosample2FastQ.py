#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from shutil import rmtree
from tempfile import mkdtemp

def parseArgs():
	parser = ArgumentParser(description='Downloads FastQ (Phred 33) reads '
		'given a BioSample accession (SAMN########) or SRA Run accession '
		'(SRR#######). Depends on Entrez Direct scripts '
		'http://www.ncbi.nlm.nih.gov/books/NBK179288/ and SRA toolkit '
		'http://www.ncbi.nlm.nih.gov/books/NBK158899/', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-a', '--acc', required=True, metavar='STR',
		help='accession that points to FastQ reads to download')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outdir', metavar='DIR',
		default=None, help='output directory [./]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	acc = opt.acc
	tmp = mkdtemp()
	tfh = os.path.join(tmp, acc)
	if opt.outdir is not None:
		outdir = os.path.abspath(opt.outdir)
	else:
		outdir = os.getcwd()

	os.system('esearch -db sra -query {} | efetch -format docsum | '
			  'grep \'<Runs>\' > {}'.format(acc, tfh))

	download = []
	with open(tfh, 'r') as run_info:
		for ln in run_info:
			if '<Run acc="' in ln:
				SRR = ln.split('<Run acc="')[1].split('"')[0]
				download.append(SRR)
	rmtree(tmp)

	if len(download) > 1:
		print 'downloading {} ...'.format(', '.join(download))
	elif len(download) == 1:
		print 'downloading {} ...'.format(download[0])
	else:
		sys.exit('No matches for {}.'.format(acc))
	for SRR in download:
		os.system('fastq-dump -O {} --dumpbase --split-files --readids -Q 33 '
				  '-defline-qual \'+\' --defline-seq \'@$ac_$sn[_$rn]/$ri\' '
				  '{} 2> /dev/null '.format(outdir, SRR))
		print '\tdownloaded {}'.format(SRR)

if __name__ == '__main__':
	main()
