#!/usr/bin/env python


import collections
import itertools
import os
import sys
from argparse import ArgumentParser
from math import sqrt

from Bio import SeqIO
from regex import findall
from scipy.stats.mstats import pearsonr

def parseArgs():
	parser = ArgumentParser(description='Given a pair of FastA nucleotide '
		'files, calculates observed and expected tetranucleotide frequencies, '
		'computes Z-scores for each tetranucleotide, and reports the Pearson '
		'correlation coefficient', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-1', '--set1', required=True, metavar='FILE',
		help='first input FastA sequence file')
	req.add_argument('-2', '--set2', required=True, metavar='FILE',
		help='second input FastA sequence file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-l', '--min-length', type=int, metavar='INT',
		default=1000, help='minimum sequence length (base pairs) to use '
		'in kmer frequency calculations [1000]')
	return parser.parse_args()

def nuc_cln(seq):
	if not len(set(seq) - set(['A', 'C', 'G', 'T'])):
		return True
	sys.stderr.write('WARNING: {} kmer sequence contains illegal '
					'characters; skipping it...\n'.format(seq))
	return False

def gen_kmers(k):
	for kmer in itertools.product(('A', 'C', 'G', 'T'), repeat=k):
		yield ''.join(kmer)  # faster than appending to string

def calc_zscores(di, tri, tet):
	# for each kmer, calc expected frequencies
	# according to the maximal order (k-1, k-2, ...) Markov model
	kmers_exp = {}
	for mer in [kmer for kmer in tet if nuc_cln(kmer)]:
		kmers_exp[mer] = 1. * tri[mer[:3]] * tri[mer[1:]] / di[mer[1:3]]

	# for each kmer, calc the stdev and Z-score
	kmers_std = {}
	kmers_z  = {}
	for mer, exp in kmers_exp.items():
		two = di[mer[1:3]]
		kmers_std[mer] = sqrt(exp * (two-tri[mer[:3]]) * (two-tri[mer[1:]]) / (two*two))
		try:
			kmers_z[mer] = (tet[mer]-exp) / kmers_std[mer]
		except ZeroDivisionError: # if an individual kmer's stdev=0
			kmers_z[mer] = 1 / (two*two)
	return kmers_z

def cnt_kmers(infile, minlen):
	di, tri, tet  = [collections.Counter() for _ in range(4 - 1)]
	for rec in SeqIO.parse(infile, 'fasta'):
		if len(rec.seq) < minlen:
			sys.stderr.write('INFO: skipping {} due to small size...\n'.format(rec.id))
			continue
		for seq in [str(rec.seq).upper(),
					str(rec.seq.reverse_complement()).upper()]:
			# NOTE: mono tallies unnecessary
			# for n in ('A', 'C', 'G', 'T'):
			# 	mono[n] += seq.count(n)
			# NOTE: count func doesnt tally overlaps so cant use for kmers
			# Prefer regex.findall overlapping func for speed (written in C); this
			# method avoids the check with nuc_cln(), so its quicker for shortmers
			for k, c in [(2, di), (3, tri)]:
				for kmer in list(gen_kmers(k)):
					c[kmer] += len(findall(kmer, seq, overlapped=True))
			# slice across the contig instead of iterating over each gen_kmers
			for i in range(len(seq) - 4 + 1):
				tet[seq[i:i + 4]] += 1
	return calc_zscores(di, tri, tet)

def main():
	opts = parseArgs()
	minlen = opts.min_length

	t1 = cnt_kmers(os.path.abspath(os.path.expanduser(opts.set1)), minlen)
	t2 = cnt_kmers(os.path.abspath(os.path.expanduser(opts.set2)), minlen)

	tetmers = sorted(list(gen_kmers(4)))
	zscores = [[t1.get(t,0) for t in tetmers], [t2.get(t,0) for t in tetmers]]
	r_val, p_val = pearsonr(zscores[0], zscores[1])
	print '{}\t{}\t{}'.format(r_val, opts.set1, opts.set2)

if __name__ == '__main__':
	main()
