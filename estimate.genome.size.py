#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from numpy import mean, std
from shutil import rmtree
from tempfile import mkdtemp

def parseArgs():
	parser = ArgumentParser(description='Estimates genome size and the '
			'percentage of the genome with duplication. Requires BBTools.',
			add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-1', '--R1', required=True, metavar='FILE',
		help='input R1 FastQ reads file')
	req.add_argument('-2', '--R2', required=True, metavar='FILE',
		help='input R2 FastQ reads file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-k', '--kmers', type=str, metavar='CSVs',
		default='21,33,55,77,99,127',
		help='comma-separated list of kmers to use [21,33,55,77,99,127]')
	opt.add_argument('-w', '--warn-std', metavar='FLOAT', type=float,
		default=5.0, help='print warning to stderr when the standard '
		'deviation exceeds a specified percentage [5.0]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	R1  = os.path.abspath(opt.R1)
	R2  = os.path.abspath(opt.R2)
	tmp = mkdtemp()

	kmers = opt.kmers.replace(' ', '').split(',')
	if len(kmers) < 3:
		sys.exit('ERROR: at least 3 kmers required to estimate genome size')

	dupes, sizes = [], []
	for k in kmers:
		peaks = os.path.join(tmp, k + '_peaks.txt')
		os.system('kmercountexact.sh ploidy=1 k={} peaks={} '
			'in1={} in2={} 2> {}'.format(k, peaks, R1, R2, os.devnull))
		i = 0
		with open(peaks) as dat:
			for ln in dat:
				if ln.startswith('#genome_size'):
					sizes.append(float(ln.split()[1]))
					i += 1
				elif ln.startswith('#percent_repeat'):
					dupes.append(float(ln.split()[1]))
					i += 1
		if i != 2:
			sys.exit('ERROR: unable to parse genome_size and percent_repeat '
				'from:  {}'.format(peaks))
	rmtree(tmp)
	avg_size = mean(sizes)/float(1000000)
	avg_dupe = mean(dupes)
	std_size = std(sizes)/float(1000000)
	std_dupe = std(dupes)
	print '{:.2f} Mbp mean estimated genome size'.format(avg_size)
	print '{:.3f} Mbp stdev estimated genome size'.format(std_size)
	print '{:.3}% mean estimated duplication in the genome'.format(avg_dupe)
	print '{:.3}% stdev estimated duplication in the genome'.format(std_dupe)
	if std_size/avg_size*100 > opt.warn_std:
		sys.stderr.write('WARNING: percentage error exceeded {}'.format(
			opt.warn_std))

if __name__ == '__main__':
	main()
