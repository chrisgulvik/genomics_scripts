#!/usr/bin/env python


import os
from argparse import ArgumentParser
from glob import glob
from sys import exit
from time import strftime
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Extracts core orthologous sequences from clustered bidirectional best hits (BBH)', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', metavar='FILE', required=True, help='tab-delimited input listing clustered sequence identifiers')
	req.add_argument('-d', '--indir', metavar='DIR', required=True, help='directory of FastA files containing sequence files used to generate the infile')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--core', metavar='FLOAT', type=float, default=1, help='fraction of samples a sequence must be in (per cluster) to be considered core [1.0]')
	opt.add_argument('-e', '--ext', metavar='STR', type=str, default='.fa', help='file extension to append to each extracted FastA file [.fa]')
	opt.add_argument('-h', '--help', action='help', help='show this help message and exit')
	opt.add_argument('-o', '--outpath', metavar='PATH', default=None, help='output directory [./BBH.clust.extracts--<date>_<time>]')
	opt.add_argument('-p', '--pref', metavar='STR', type=str, default='', help='prefix to discard before extracting sequence names from input FastA files [None]')
	opt.add_argument('-s', '--suff', metavar='STR', type=str, default='.faa', help='suffix to discard when extracting sequence names from input FastA files [.faa]')
	return parser.parse_args()

def main():
	opts = parseArgs()
	# Grab sample ID names from files within input dir
	sample_ids = []
	sample_files = glob(os.path.join(opts.indir, opts.pref+'*'+opts.suff))
	for s in sample_files:
		sample_ids.append(os.path.basename(s).lstrip(opts.pref).rstrip(opts.suff))
	if len(sample_ids) != len(sample_files):
		exit('ERROR: sample IDs unable to be parsed. '
			'Revise input directory (-d) or \'-p\' and \'-s\' arguments.')
	if opts.outpath is not None:
		outpath = os.path.abspath(opts.outpath)
	else:
		autocreated_dirname = 'BBH.clust.extracts--' + strftime('%d%b%Y_%-I:%M%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)
	for p in [os.path.join(outpath, s) for s in sample_ids]:
		if not os.path.exists(p):
			os.makedirs(p)

	# Get core sequence identifiers
	clust, filt, par = 0, 0, 0
	core_ids = []
	with open(opts.infile) as f:
		for line in f:
			clust += 1
			l = line.rstrip('\n').split('\t')
			if len([s for s in l if s != ''])/float(len(l)) >= opts.core:
				if any(',' in s for s in l):
					par += 1
				else:
					filt += 1
					core_ids.append(l)
	trs_ids = sorted([sorted([row[i] for row in core_ids]) for i in range(len(core_ids[1]))])
	print '{} total clusters'.format(clust)

	# Extract core sequences
	for sid, sample_file, sample in zip(sample_ids, sample_files, trs_ids):
		mfasta_idx = SeqIO.index(sample_file, 'fasta')
		for locus in sample:
			with open(os.path.join(outpath, sid, locus+opts.ext), 'w') as o:
				o.write(mfasta_idx.get_raw(locus))
	print '{} cluster(s) of core homologs have at least one set of '\
			'paralogs'.format(par)
	print '{} ({}%) core orthologs found per sample'.format(
			filt, round(float(filt)/clust*100, 3))

if __name__ == '__main__':
	main()
