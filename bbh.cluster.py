#!/usr/bin/env python


import os
import shutil
import subprocess as sp
import sys
from argparse import ArgumentParser
from cmath import sqrt
from decimal import Decimal, ROUND_HALF_UP
from glob import glob
from multiprocessing import cpu_count
from tempfile import mkdtemp
from time import strftime

def parseArgs():
	parser = ArgumentParser(description='Lists and quantifies homologous'
		' sequence clusters from bidirectional best hits (BBH).',
		epilog='Note: Output is saved as bbh.clust.{tab,stats.txt}, so if'
		' iteratively testing parameters on the same dataset (e.g., inflation'
		' effects), rename it after each or it will be overwritten.')
	subparsers = parser.add_subparsers(title='subcommands',
		metavar='<input type>',
		description='One of the two is mandatory to specify input type.'
		' Executing either subcommand along with --help lists their optional'
		' and mandatory parameters.')
	a = subparsers.add_parser('DIR', help='Input directory containing pairs'
		' of BBH alignment files to cluster. Sample names are extracted from'
		' filenames and are used as query prefixes to locate sequences, i.e.,'
		' sequence identifiers must begin with the sample name.')
	a.add_argument('-c', '--data-column', metavar='INT', type=int, default=12,
		help='column number of data values to use for clustering; Default, bit'
		' score [12]')
	a.add_argument('-i', '--indir', metavar='PATH', required=True,
		help='input directory containing BBH files to cluster')
	a.add_argument('-k', '--keep', default=False, action='store_true',
		help='keep temporary files in outpath')
	a.add_argument('-p', '--pref', metavar='STR', type=str, default='bbh.',
		help='prefix of pairwise BBH filenames [bbh.]')
	a.add_argument('-s', '--suff', metavar='STR', type=str,
		default='.filt.tab', help='suffix of pairwise BBH filenames'
		' [.filt.tab]')
	a.add_argument('-d', '--delim', metavar='STR', type=str, default=',',
		help='delimiter between sample names in pairwise BBH filenames [,]')
	a.add_argument('-I', '--inflation', metavar='FLOAT', type=float,
		default=1.5, help='main inflation value for Markov clustering [1.5]')
	a.add_argument('-o', '--outpath', metavar='PATH', default=None,
		help='output directory [./BBH.clust--<date>_<time>]')
	a.add_argument('-t', '--threads', metavar='INT', type=int, default=0,
		help='number of threads [all]')
	a.add_argument('-x', '--xoptmcl', metavar='\'STR\'', type=str,
		default=None, help='extra commands to pass to mcl (e.g,.'
		' \'-scheme 7\' or\n\'-pct 95\') [none]')
	b = subparsers.add_parser('ABC', help='ABC format file summarizing'
		' all-vs-all BBHs')
	b.add_argument('-a', '--abc', metavar='FILE', required=True,
		help='input ABC format file summarizing all-vs-all BBHs to cluster')
	b.add_argument('-k', '--keep', default=False, action='store_true',
		help='keep temporary files in outpath')
	b.add_argument('-n', '--names', metavar='STR,STR,STR[...]', required=True,
		help='comma-delimited sample names present in ABC input file as well'
		' as prefixes to each sequence identifier (nodes and edges)')
	b.add_argument('-I', '--inflation', metavar='FLOAT', type=float,
		default=1.5, help='main inflation value for Markov clustering [1.5]')
	b.add_argument('-o', '--outpath', metavar='PATH', default=None,
		help='output directory [./BBH.clust--<date>_<time>]')
	b.add_argument('-t', '--threads', metavar='INT', type=int, default=0,
		help='number of threads [all]')
	b.add_argument('-x', '--xoptmcl', metavar='\'STR\'', type=str,
		default=None, help='extra commands to pass to mcl (e.g,.'
		' \'-scheme 7\' or\n\'-pct 95\') [none]')
	return parser.parse_args()

def require_dependency(dep):
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, dep)) and \
		not os.path.isdir(os.path.join(path, dep)):
			return True
	sys.stderr.write('ERROR: {} unavailable; not in $PATH\n'.format(dep))
	sys.exit(1)

def calc_sample_nr_from_pairwise_cnt(pairwise_file_count):
	# NOTE: equation is x^2 - x - n*2 = 0, where n is bbh file combos
	a, b, c = 1, -1, pairwise_file_count * -2
	discriminant = (b ** 2) - (4 * a * c)
	solution1 = (-b - abs(sqrt(discriminant))) / (2 * a)
	solution2 = (-b + abs(sqrt(discriminant))) / (2 * a)
	soln1 = int(Decimal(solution1).quantize(Decimal('1'),
		rounding=ROUND_HALF_UP))
	soln2 = int(Decimal(solution2).quantize(Decimal('1'),
		rounding=ROUND_HALF_UP))
	return [x for x in (soln1, soln2) if x > 0]

def main():
	opt = parseArgs()
	require_dependency('mcxload')
	require_dependency('mcl')

	# I/O handling
	tmp = mkdtemp()
	if opt.outpath is not None:
		outpath = os.path.realpath(os.path.expanduser(opt.outpath))
	else:
		autogen_dir = 'BBH.clust--' + strftime('%d%b%Y_%-I:%M%p').upper()
		outpath = os.path.join(os.getcwd(), autogen_dir)
	if not os.path.exists(outpath):
		os.mkdir(outpath)

	# Number of CPUs to use
	if opts.threads < 1:
		cpus = str(cpu_count())
	else:
		cpus = str(opts.threads)

	# Parse sample names and ABC file handling
	if opt.indir:
		pref = opt.pref.strip('\'').strip('"')
		suff = opt.suff.strip('\'').strip('"')
		delim = opt.delim.strip('\'').strip('"')
		indir = os.path.join(os.path.abspath(opt.indir))

		pair_files = glob(indir, pref + '*' + suff)
		if len(pair_files) < 3:
			sys.stderr.write('ERROR: too few (< 3) files found\n'.format(
				pair_files))

		# Estimate sample number(s) from quantity of bbh file quantity
		estimated_sample_cnts = calc_sample_nr_from_pairwise_cnt(
			len(pair_files))

		# Parse filenames to get unique sample names and generate ABC file
		sample_names = set()
		with open(os.path.join(tmp, 'bbh.groupsummary.abc'), 'w') as o:
			for pair_file in pair_files:
				b = os.path.basename(pair_file)
				s1, s2 = b.lstrip(pref).rstrip(suff).split(delim)
				sample_names.update([s1, s2])
				with open(pair_file) as pf:
					for line in pf:
						l = line.split('\t')
						o.write(l[0]+'\t'+l[1]+'\t'+l[opt.data_column-1]+'\n')
		if not any(x == len(sample_names) for x in estimated_sample_cnts):
			sys.stderr.write('ERROR: the observed sample name quantity ({})'
				' is unequal to the expected sample quantity from {} parsed'
				' bbh filenames\n'.format(len(sample_names),
					' or '.join(str(x) for x in estimated_sample_cnts),
					indir))
			sys.exit(1)
	elif opt.abc:
		with open(os.path.realpath(os.path.expanduser(opt.abc))) as f:
			if next(f).count('\t') != 3:
				sys.stderr.write('ERROR: expected 3-column tab-delim input\n')
				sys.exit()
		abc_file = os.path.realpath(os.path.expanduser(opt.abc))
		sample_names = opt.names.strip('\'').strip('"').split(',')

	# Execute clustering
	graph_file = os.path.join(tmp, 'graph.mci')
	seqids_file = os.path.join(tmp, 'seqids.tab')
	mcl_file = os.path.join(tmp, 'clust.mcl')
	cmd_graph = ['mcxload', '--stream-mirror', '-re', 'max', '-abc', abc_file,
				'-o', graph_file, '--write-binary', '-write-tab', seqids_file]
	cmd_clust = ['mcl', graph_file, '-I', str(opt.inflation), '-V', 'all',
				'-use-tab', seqids_file, '-o', mcl_file, '-te', cpus]
	if opt.xoptmcl:
		cmd_clust.extend(opt.xoptmcl.strip('\'').strip('"').split(' '))
	for cmd in [cmd_graph, cmd_clust]:
		with open(os.devnull, 'wb') as dump:
			return_code = sp.Popen(cmd, shell=False, stdout=dump, stderr=dump)
			if return_code.wait() != 0:
				sys.stderr.write('ERROR: failed system call\n{}\n'.format(
					' '.join(cmd)))
				sys.exit()

	# Parse mcl output
	with open(mcl_file) as mcl_data, \
	open(os.path.join(outpath, 'bbh.clust.tab'), 'w') as o: 
		i = 0
		for line in mcl_data:
			l = line.rstrip('\n').split('\t')
			match = ''
			for n in sample_names:
				match = match+','.join([m for m in l if m.startswith(n)])+'\t'
			o.write(match[:-1]+'\n' if match.endswith('\t') else match+'\n')
			i += 1
	with open(seqids_file) as f:
		seqs = sum(1 for l in f)
	with open(os.path.join(outpath, 'bbh.clust.stats.txt'), 'w') as stats:
		stats.write('Queried {} samples\nQueried {} total sequences\n' \
			'Clustered {} sequence groups\n'.format(len(sample_names), seqs,
				i))

	# Optionally keep intermediate files
	if opt.keep:
		if opt.indir:
			shutil.copy(os.path.join(tmp, 'bbh.groupsummary.abc'),
				os.path.join(outpath, 'bbh.groupsummary.abc'))
		for f in ['clust.mcl', 'graph.mci', 'seqids.tab']:
			shutil.copy(os.path.join(tmp, f), os.path.join(outpath, f))
	shutil.rmtree(tmp)

if __name__ == '__main__':
	main()
