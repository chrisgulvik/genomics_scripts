#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from itertools import groupby
from operator import itemgetter
from shutil import rmtree
from tempfile import mkdtemp
from Bio import SeqIO, Alphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation

def parseArgs():
	parser = ArgumentParser(description='Creates a GenBank file from a '
	'FastA file', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--input', required=True, metavar='FILE',
		help='input FastA file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--output', metavar='FILE', default=None,
		help='output GenBank file [stdout]')
	opt.add_argument('--label-ambiguous-features', action='store_true',
		default=False, help='label positions where N occurs as \'ambiguous\' '
		'[off]')
	return parser.parse_args()

def collapse_integers_to_ranges(sites):
	''' converts list of integers into list of pairs of begin,end sites '''
	ranges = []
	for _, g in groupby(enumerate(sites), key=lambda e: e[0]-e[1]):
		group = [val[1] for val in g]
		if len(group) > 1:
			ranges.append(group[::len(group) - 1])
		elif len(group) == 1:
			ranges.append([group[0], group[0]])
	return ranges

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.input))
	if opt.output is not None:
		ofh = os.path.abspath(os.path.expanduser(opt.output))
	else:
		ofh = sys.stdout

	if not opt.label_ambiguous_features:
		SeqIO.convert(ifh, 'fasta', ofh, 'genbank',
			alphabet=Alphabet.generic_dna)
		sys.exit(0)

	tmp = mkdtemp()
	temp_gbk = os.path.join(tmp, 'gb')
	SeqIO.convert(ifh, 'fasta', temp_gbk, 'genbank',
		alphabet=Alphabet.generic_dna)

	records = []
	with open(temp_gbk) as i:
		for rec in SeqIO.parse(i, 'genbank'):
			# Find ambiguous N nucleotides
			i, split_len = 0, 0
			sites = []
			seq = rec.seq.upper()
			while i >= 0:
				i = seq.find('N')
				if i >= 0:
					sites.append(i + split_len)
					split_len += i+1
					seq = seq[i+1:]

			# Label ambiguous sites in GenBank file
			ranges = collapse_integers_to_ranges(sites)
			for begin_end in ranges:
				sf = SeqFeature(FeatureLocation(begin_end[0], begin_end[1]+1),
					type='ambiguous')
				rec.features.append(sf)
			records.append(rec)

	SeqIO.write(records, ofh, 'genbank')
	sys.stderr.write('INFO: {} sequence records written\n'.format(
		len(records)))
	rmtree(tmp)

if __name__ == '__main__':
	main()
