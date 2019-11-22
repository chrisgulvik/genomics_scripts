#!/usr/bin/env python


import multiprocessing as mp
import os
import sys
from argparse import ArgumentParser
from shutil import rmtree
from tempfile import mkdtemp

def parseArgs():
	parser = ArgumentParser(description='Downloads FastQ (Phred 33) reads'
		' given a BioSample accession (SAMN########) or SRA Run accession'
		' (SRR#######). Depends on Entrez Direct scripts'
		' http://www.ncbi.nlm.nih.gov/books/NBK179288/ and SRA toolkit'
		' http://www.ncbi.nlm.nih.gov/books/NBK158899/', add_help=False,
		epilog='NOTE: EDirect ver 7.50+ uses the API key value from the'
		' NCBI_API_KEY environment variable')
	req = parser.add_argument_group('Required')
	req.add_argument('-a', '--acc', required=True, metavar='STR',
		help='accession that points to FastQ reads to download')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outdir', metavar='DIR',
		default=None, help='output directory [./]')
	# opt.add_argument('--api-key', metavar='STR',
	# 	default=None, help='API key (64 characters) from NCBI; re-assigns NCBI_API_KEY'
	# 	' environment variable [None]')
	return parser.parse_args()

def check_dependency(dep):
	''' returns a tuple where first item is boolean and second item is the
	dependency path; symbolic links are not expanded to their real path '''
	found = False, None
	for path in os.environ['PATH'].split(os.pathsep):
		dependency_path = os.path.join(path, dep)
		if all([os.path.exists(dependency_path),
		os.path.isfile(dependency_path),
		os.access(dependency_path, os.X_OK)]):
			return True, dependency_path
	return found

def main():
	opt = parseArgs()
	acc = opt.acc
	tmp = mkdtemp()
	tfh = os.path.join(tmp, acc)
	if opt.outdir is not None:
		outdir = os.path.realpath(opt.outdir)
	else:
		outdir = os.getcwd()

	if not all([check_dependency('esearch')[0],
	check_dependency('efetch')[0]]):
		sys.stderr.write('ERROR: Entrez Direct scripts not found in $PATH.'
			' Install edirect.\n')
		sys.exit(1)
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
		print('downloading {} ...'.format(', '.join(download)))
	elif len(download) == 1:
		print('downloading {} ...'.format(download[0]))
	else:
		sys.exit('No matches for {}.'.format(acc))
	
	if check_dependency('fasterq-dump')[0]:
		threads_avail = mp.cpu_count()
		# https://github.com/ncbi/sra-tools/issues/242#issuecomment-549525056
		if threads_avail < 6:
			threads = threads_avail
		else:
			threads = 6
		for SRR in download:
			os.system('fasterq-dump {} -O {} --force --split-files'
					  ' --print-read-nr --threads {} --temp {} 1> {}'.format(
						SRR, outdir, threads, tmp, os.devnull))
			print('  downloaded {}'.format(SRR))
	elif check_dependency('fastq-dump')[0]:
		for SRR in download:
			os.system('fastq-dump -O {} --dumpbase --split-files --readids'
					  ' -Q 33 -defline-qual \'+\' '
					  ' --defline-seq \'@$ac_$sn[_$rn]/$ri\' {} 1> {}'.format(
						outdir, SRR, os.devnull))
			print('  downloaded {}'.format(SRR))
	else:
		sys.stderr.write('ERROR: fastq-dump or fasterq-dump is required but'
			'both are absent in $PATH. Install sratoolkit.\n')
		sys.exit(1)

if __name__ == '__main__':
	main()
