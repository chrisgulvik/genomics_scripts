#!/usr/bin/env python


import gzip
import os
import re
import sys
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Alphabet import generic_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
	parser = ArgumentParser(description='Extracts sites from each record in a'
	' multi-record FastA file', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input multi-record FastA file, optionally gunzip compressed')
	req.add_argument('-s', '--sites', required=True, nargs='+', type=str,
		metavar='INT|RANGE', help='position(s) to extract; range requires a'
		' colon or hyphen; order given is the concatenated output order')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', default=None,
		help='FastA output [stdout]')
	opt.add_argument('--unequal-lengths', action='store_true', default=False,
		help='allow unequal record lengths [off]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.infile))
	positions = opt.sites

	# Quick checks of positions
	inv_pattern = re.compile('[^0-9:-]')
	if any([inv_pattern.search(x) for x in positions]):
		sys.stderr.write('ERROR: sites must only contain [0-9:-]\n')
		sys.exit(1)
	if any([x.startswith('0') for x in positions]):
		sys.stderr.write('ERROR: 0 position prefix found; site extraction is'
			' 1-based\n')
		sys.exit(1)

	# Extract all specified sites one file at a time
	if ifh.endswith('.gz'):
		ifh = gzip.open(ifh)
	records = []
	rec_length = len(next(SeqIO.parse(ifh, 'fasta')))
	for rec in SeqIO.parse(ifh, 'fasta'):
		if rec_length != len(rec):
			sys.stderr.write('WARNING: unequal record lengths\n')
			if not opt.unequal_lengths:
				sys.exit(1)
		rec_length = len(rec)
		seq = ''
		for pos in positions:
			if not any((x in (':', '-')) for x in pos):
				pos = int(pos) - 1
				if pos >= rec_length:
					sys.stderr.write('ERROR: {} value larger than {} record'
						' ({} sites).\n'.format(pos, rec.id, rec_length))
					sys.exit(1)
				seq += str(rec.seq[pos])
			elif ':' in pos:
				start, stop = [int(x) for x in pos.split(':')]
				if any(v > rec_length for v in (start, stop)):
					sys.stderr.write('ERROR: {} range larger than {} record'
						' ({} sites).\n'.format(pos, rec.id, rec_length))
					sys.exit(1)
				if start > stop:
					sys.stderr.write('ERROR: {} range must have second value'
						' exceed first.\n'.format(pos))
					sys.exit(1)
				seq += str(rec.seq[start-1:stop])
			elif '-' in pos:
				start, stop = [int(x) for x in pos.split('-')]
				if any(v > rec_length for v in (start, stop)):
					sys.stderr.write('ERROR: {} range larger than {} record'
						' ({} sites).\n'.format(pos, rec.id, rec_length))
					sys.exit(1)
				if start >= stop:
					sys.stderr.write('ERROR: {} range must have second value'
						' exceed first.\n'.format(pos))
					sys.exit(1)
				seq += str(rec.seq[start-1:stop])
			else:
				sys.stderr.write('ERROR: {} not detected a single position or'
					' a range of positions\n'.format(pos))
				sys.exit(1)
		defline = str(rec.description) + ' Sites(' + ' '.join(positions) + ')'
		record = SeqRecord(Seq(seq, generic_alphabet), id=defline,
			description='')
		records.append(record)

	# Output specified positions with FastA format
	if opt.outfile is not None:
		outfile = os.path.abspath(os.path.expanduser(opt.outfile))
		SeqIO.write(records, outfile, 'fasta')
	else:
		try:
			SeqIO.write(records, sys.stdout, 'fasta')
		except Exception, errno:
			if errno == 32: #IOError for py2, BrokenPipeError for py3
				pass

if __name__ == '__main__':
	main()
