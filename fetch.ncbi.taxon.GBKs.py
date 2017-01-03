#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from shutil import rmtree
from subprocess import check_output
from tempfile import mkdtemp

def parseArgs():
	parser = ArgumentParser(description='Fetches GenBank files for all '
		'genomes within a specified taxon. Requires NCBI\'s E-utilities.',
		epilog='Note that although the proper assembly status term to '
		'describe completed genomes is \'Complete Genome\', the single word '
		'\'Complete\' is used here. When more than status is desired, the '
		'--status opt needs to be specified each time. For example, -s '
		'Complete -s Contig.', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--tax-id', required=True, metavar='INT',
		help='Taxonomy ID from NCBI. '
		'See https://www.ncbi.nlm.nih.gov/taxonomy')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outpath', metavar='PATH',
		default=None, help='output directory where GenBank files will be '
		'downloaded to [./]')
	opt.add_argument('-s', '--status', metavar='STR', action='append',
		choices=['Chromosome', 'Contig', 'Complete', 'Scaffold'],
		default=[], help='assembly status(es) of files to fetch '
		'[Chromosome, Complete, Contig, Scaffold]')
	return parser.parse_args()

def main():
	opt   = parseArgs()
	tmp   = mkdtemp()
	taxid = 'txid' + str(opt.tax_id)
	if not opt.status:  #Default is empty list to indicate fetch all
		status_want = ['Chromosome', 'Contig', 'Complete', 'Scaffold']
	else:
		status_want = opt.status
	if opt.outpath is not None:
		out = os.path.abspath(opt.outpath)
	else:
		out = os.getcwd()
	if not os.path.exists(out):
		os.mkdir(out)

	t = check_output('esearch -db taxonomy -query "{}[ORGN:exp]" | '
		'efetch -format docsum | xtract -pattern DocumentSummary -element '
		'ScientificName'.format(taxid), shell=True).rstrip('\n')
	print 'Fetching GenBank files for {} ...'.format(t)
	tmp_file = os.path.join(tmp, 'urls')
	os.system('esearch -db assembly -query "{}[ORGN:exp]" | '
		'efetch -format docsum | xtract -pattern DocumentSummary -element '
		'AssemblyStatus FtpPath_RefSeq 1> {}'.format(taxid, tmp_file))
	i = 0
	with open(tmp_file) as urls:
		for ln in urls:
			status, part_url = ln.rstrip('\n').split('\t')
			if any([x in status for x in status_want]):
				s = part_url.split('/')[-1]
				gbk_url = os.path.join(part_url, s + '_genomic.gbff.gz')
				gz_temp = os.path.join(tmp, s + '.gbk.gz')
				gbk_out = os.path.join(out, s + '.gbk')
				os.system('wget -q -O {} {}'.format(gz_temp, gbk_url))
				if os.path.getsize(gz_temp) > 0:
					os.system('gunzip -c {} > {}'.format(gz_temp, gbk_out))
					i += 1
				else:
					sys.exit('ERROR: unable to download {}'.format(s))	
	rmtree(tmp)
	print 'Downloaded {} GenBank files'.format(i)

if __name__ == '__main__':
	main()
