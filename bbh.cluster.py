#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from glob import glob
from subprocess import Popen
from time import strftime

def parseArgs():
	parser = ArgumentParser(description='\nLists and quantifies homologous sequence clusters from bidirectional best\nhits (BBH).', formatter_class=RawTextHelpFormatter,
	epilog='Note: Output is saved as <outpath>/bbh.clust.{tab,stats.txt}, so if\niteratively testing parameters on the same dataset (e.g., inflation effects),\nrename it after each or it will be overwritten.')
	subparsers = parser.add_subparsers(title='subcommands', metavar='<input type>', description='One of these is mandatory to specify input type. Executing either subcommand\nalong with --help lists their optional and mandatory parameters.')
	a = subparsers.add_parser('DIR', help='Input directory containing BBH files to cluster.\nSample names are extracted from filenames and are used\nas query prefixes to locate sequences, i.e., sequence\nidentifiers must begin with the sample name.')
	a.add_argument('-i', '--indir', metavar='PATH', required=True, help='input directory containing BBH files to cluster')
	a.add_argument('-k', '--keep', default=False, action='store_true', help='keep temporary files in outpath')
	a.add_argument('-p', '--pref', metavar='\'STR\'', type=str, default='\'bbh.\'', help='file prefix to strip from each BBH file for extracting sample ID names [\'bbh.\']')
	a.add_argument('-s', '--suff', metavar='\'STR\'', type=str, default='\'.filt.tab\'', help='file suffix to strip from each BBH file for extracting sample ID names [\'.filt.tab\']')
	a.add_argument('-d', '--delim', metavar='\'STR\'', type=str, default='\',\'', help='delimiter between two sample ID names in each BBH file for extracting sample ID names [\',\']')
	a.add_argument('-I', '--inflation', metavar='FLOAT', type=float, default=1.5, help='main inflation value for Markov clustering [1.5]')
	a.add_argument('-o', '--outpath', metavar='PATH', default=None, help='output directory [./BBH.clust--<date>_<time>]')
	a.add_argument('-t', '--threads', metavar='INT', type=int, default='1', help='number of threads [1]')
	a.add_argument('-x', '--xoptmcl', metavar='\'STR\'', type=str, default=None, help='extra commands to pass to mcl (e.g,. \'-scheme 7\' or\n\'-pct 95\') [none]')
	b = subparsers.add_parser('ABC', help='ABC format file summarizing all-vs-all BBHs')
	b.add_argument('-a', '--abc', metavar='FILE', required=True, help='input ABC format file summarizing all-vs-all BBHs to cluster; passed as \'--abc <file>\' to mcxload')
	b.add_argument('-k', '--keep', default=False, action='store_true', help='keep temporary files in outpath')
	b.add_argument('-n', '--names', metavar='\'STR,STR,STR[...]\'', required=True, help='comma-delimited sample names present in ABC input file as well as prefixes to each sequence identifier (nodes and edges)')
	b.add_argument('-I', '--inflation', metavar='FLOAT', type=float, default=1.5, help='main inflation value for Markov clustering [1.5]')
	b.add_argument('-o', '--outpath', metavar='PATH', default=None, help='output directory [./BBH.clust--<date>_<time>]')
	b.add_argument('-t', '--threads', metavar='INT', type=int, default='1', help='number of threads [1]')
	b.add_argument('-x', '--xoptmcl', metavar='\'STR\'', type=str, default=None, help='extra commands to pass to mcl (e.g,. \'-scheme 7\' or\n\'-pct 95\') [none]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	if opt.outpath is not None:
		outpath = os.path.abspath(opt.outpath)
		if not os.path.exists(outpath):
			os.mkdir(outpath)
	else:
		autocreated_dirname = 'BBH.clust--' + strftime('%d%b%Y_%-I:%M%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)
		os.mkdir(outpath)
	if opt.indir:  
		# Summarize BBH pair files into single file
		abc_file     = os.path.join(outpath, 'bbh.groupsummary.abc')
		pair_files   = glob(os.path.join(os.path.abspath(opt.indir), opt.pref[1:-1]+'*'+opt.suff[1:-1]))
		sample_names = set()
		with open(abc_file, 'w') as o:
			for pair_file in pair_files:
				pref = opt.pref[1:-1]
				suff = opt.suff[1:-1]
				s1, s2 = os.path.basename(pair_file).lstrip(pref).rstrip(suff).split(opt.delim[1:-1])
				sample_names.update([s1, s2])
				with open(pair_file) as pf:
					for line in pf:
						l = line.split('\t')
						o.write(l[0]+'\t'+l[1]+'\t'+l[4]+'\n')
	elif opt.abc:
		with open(os.path.abspath(opt.abc)) as f:
			if next(f).count('\t') != 3:
				sys.stderr.write('ERROR: expected 3-column tab-delimited input\n')
				sys.exit()
		abc_file     = os.path.abspath(opt.abc)
		sample_names = opt.names[1:-1].split(',')

	# Execute clustering
	cmd_graph = ['mcxload', '--stream-mirror', '-re', 'max', '-abc', abc_file,
				'-o', os.path.join(outpath, 'graph.mci'), '--write-binary',
				'-write-tab', os.path.join(outpath, 'seqids.tab')]
	cmd_clust = ['mcl', os.path.join(outpath, 'graph.mci'), '-I', str(opt.inflation),
				'-use-tab', os.path.join(outpath, 'seqids.tab'), '-V', 'all',
				'-o', os.path.join(outpath, 'clust.mcl'), '-te', str(opt.threads)]
	if opt.xoptmcl:
		cmd_clust.extend(opt.xoptmcl[1:-1].split(' '))
	for cmd in [cmd_graph, cmd_clust]:
		with open(os.devnull, 'wb') as dump:
			return_code = Popen(cmd, shell=False, stdout=dump, stderr=dump)
			if return_code.wait() != 0:
				sys.stderr.write('ERROR: failed system call\n{}\n'.format(' '.join(cmd)))
				sys.exit()

	# Parse mcl output
	with open(os.path.join(outpath, 'clust.mcl')) as mcl_data, \
	open(os.path.join(outpath, 'bbh.clust.tab'), 'w') as o: 
		i = 0
		for line in mcl_data:
			l = line.rstrip('\n').split('\t')
			match = ''
			for n in sample_names:
				match = match + ','.join([m for m in l if m.startswith(n)]) + '\t'
			o.write(match[:-1]+'\n' if match.endswith('\t') else match+'\n')
			i += 1
	with open(os.path.join(outpath, 'seqids.tab'), 'r') as f:
		seqs = sum(1 for l in f)
	with open(os.path.join(outpath, 'bbh.clust.stats.txt'), 'w') as stats:
		stats.write('Queried {} samples\nQueried {} total sequences\n' \
			'Clustered {} sequence groups\n'.format(len(sample_names), seqs, i))

	# Optionally keep intermediate files
	if not opt.keep:
		if opt.indir: os.remove(abc_file)
		for f in ['clust.mcl', 'graph.mci', 'seqids.tab']:
			os.remove(os.path.join(outpath, f))

if __name__ == '__main__':
	main()
