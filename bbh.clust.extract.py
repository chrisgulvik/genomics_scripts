#!/usr/bin/env python


import os
from argparse import ArgumentParser
from glob import glob
from time import strftime
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Extracts core orthologous sequences from clustered bidirectional best hits (BBH)', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', metavar='FILE', required=True, help='tab-delimited input listing clustered sequence identifiers')
	req.add_argument('-d', '--indir', metavar='DIR', required=True, help='directory of FastA files containing sequence files used to generate the infile')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-b', '--by-sample', action='store_true', default=False, help='output sequence files according to sample name rather than by orthologous cluster group')
	opt.add_argument('-c', '--core', metavar='FLOAT', type=float, default=1, help='fraction of samples a sequence must be in (per cluster) to be considered core [1.0]')
	opt.add_argument('-e', '--ext', metavar='STR', type=str, default='.fa', help='file extension to append to each extracted FastA file [.fa]')
	opt.add_argument('-h', '--help', action='help', help='show this help message and exit')
	opt.add_argument('-l', '--paralogs', action='store_true', default=False, help='output clusters with paralogs as well')
	opt.add_argument('-o', '--outpath', metavar='PATH', default=None, help='output directory [./BBH.clust.extracts--<date>_<time>]')
	opt.add_argument('-p', '--pref', metavar='STR', type=str, default='', help='prefix to discard before extracting sequence names from input FastA files [None]')
	opt.add_argument('-s', '--suff', metavar='STR', type=str, default='.faa', help='suffix to discard when extracting sequence names from input FastA files [.faa]')
	return parser.parse_args()

def main():
	opt = parseArgs()

	# Grab sample ID names from files within input dir
	sample_ids = []
	sample_files = sorted(glob(os.path.join(opt.indir, opt.pref+'*'+opt.suff)))
	for s in sample_files:
		sample_ids.append(os.path.basename(s).lstrip(opt.pref).rstrip(opt.suff))
	if opt.outpath is not None:
		outpath = os.path.abspath(opt.outpath)
	else:
		autocreated_dirname = 'BBH.clust.extracts--' + strftime('%d%b%Y_%-I:%M%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)

	if opt.by_sample:
		for p in [os.path.join(outpath, s) for s in sample_ids]:
			if not os.path.exists(p):
				os.makedirs(p)

		# Get core sequence identifiers
		clust, ortho, par = 0, 0, 0
		core_ids = []
		with open(opt.infile) as f:
			for line in f:
				clust += 1
				l = line.rstrip('\n').split('\t')
				if len([s for s in l if s != ''])/float(len(l)) >= opt.core:
					if any(',' in s for s in l):
						par += 1
						if opt.paralogs:
							core_ids.append(l)
					else:
						ortho += 1
						core_ids.append(l)
		trs_ids = sorted([sorted([row[i] for row in core_ids]) for i in range(len(core_ids[1]))])

		# Extract core sequences
		for sid, sample_file, sample in zip(sample_ids, sample_files, trs_ids):
			mfa_idx = SeqIO.index(sample_file, 'fasta')
			for locus in sample:
				if ',' in locus:
					tags = locus.split(',')
					for t in tags:
						with open(os.path.join(outpath, sid, t+opt.ext), 'w') as o:
							o.write(mfa_idx.get_raw(t))
				else:
					with open(os.path.join(outpath, sid, locus+opt.ext), 'w') as o:
						o.write(mfa_idx.get_raw(locus))

	else:
		# Read in all FastA files
		mfasta_idx = {}
		for faa in sample_files:
			faa_idx = SeqIO.index(faa, 'fasta')
			mfasta_idx.update(faa_idx)

		# Get sequences
		clust, ortho, par = 0, 0, 0
		with open(opt.infile) as f:
			for line in f:
				clust += 1
				l = line.rstrip('\n').split('\t')
				if len([s for s in l if s != ''])/float(len(l)) >= opt.core:
					if any(',' in s for s in l):
						par += 1
						if opt.paralogs:
							p = os.path.join(outpath, 'para-clust_' + str(par))
							if not os.path.exists(p):
								os.makedirs(p)
							for tag in l:
								if ',' in tag:
									tags = tag.split(',')
									for t in tags:
										SeqIO.write(mfasta_idx[t],
											os.path.join(p, t + opt.ext), 'fasta')
								else:
									SeqIO.write(mfasta_idx[tag],
										os.path.join(p, tag + opt.ext), 'fasta')
					else:
						ortho += 1
						p = os.path.join(outpath, 'ortho-clust_' + str(ortho))
						if not os.path.exists(p):
							os.makedirs(p)
						for tag in l:
							SeqIO.write(mfasta_idx[tag],
								os.path.join(p, tag + opt.ext), 'fasta')

	# Summarize counts
	print '{} total clusters'.format(clust)
	print '{} clusters of core homologs have at least one paralog'.format(par)
	print '{} ({}% of input clusters) core orthologs found per sample'.format(
			ortho, round(float(ortho)/clust*100, 3))

if __name__ == '__main__':
	main()
