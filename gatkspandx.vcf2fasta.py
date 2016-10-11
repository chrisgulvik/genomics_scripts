#!/usr/bin/env python


import sys
from argparse import ArgumentParser
from collections import Counter


def parseArgs():
	parser = ArgumentParser(description='Parses a VCF file processed by GATK '
			'within SPANDx and prints a FastA to standard out. This is '
			'particularly useful when manipulating the out.vcf file. '
			'Multiallelic SNPs are handled, heterozygous SNPs are filtered, '
			'FILTER tags are ignored, and unmapped loci are reported as gaps. '
			'Sample names and their binary data must begin on the 10th '
			'column. ',
			add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--input', required=True, metavar='FILE',
		help='input VCF file from Parsnp')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	return parser.parse_args()

# Main
opts = parseArgs()

# Get sample IDs
with open(opts.input, 'r') as ifh:
	for ln in ifh:
		if any(ln.startswith('#CHROM') for x in ln):
			samples_d = {k:[] for k in ln.rstrip('\n').split('\t')[9:]}
			samples_l = ln.rstrip('\n').split('\t')[9:]
			break

# Parse the VCF file
ifh = open(opts.input, 'r')
for ln in ifh:
	if not ln.startswith('#'):
		data = ln.rstrip('\n').split('\t')
		ref = data[3]
		if ',' in data[4]:
			alt = data[4].split(',')
		else:
			alt = [data[4]]

		# filter genotype
		gt = 0
		for x in data[9:]:
			y, z = x.split(':')[0].split('/')
			if y != z:
				gt += 1

		if gt == 0:
			binary_data = [x.split(':')[0].split('/')[0] for x in data[9:]]
			for i, datum in enumerate(binary_data):
				if datum == '0':
					nuc = ref
				elif datum == '1':
					nuc = alt[0]
				elif datum == '2':
					nuc = alt[1]
				elif datum == '3':
					nuc = alt[2]
				elif datum == '.':
					nuc = '-'
				else:
					sys.stderr.write('ERROR: ALT value >3\n')
				sample_id = samples_l[i]
				samples_d[sample_id].append(nuc)
ifh.close()

# Write FastA to Std Out
snp_len = len(samples_d[sample_id])
for sample in samples_l:
	snps = ''.join(samples_d[sample])
	if len(snps) != snp_len:
		sys.stderr.write('ERROR: unequal sequence records printed\n')
	else:
		print '>{}'.format(sample)
		print snps
