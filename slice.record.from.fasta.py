#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(add_help=False,
		description='Finds FastA file sequence record(s) that contain a given '
		'string and extracts the specified range of sites.',
		epilog='NOTE: BLAST outputs 1-based closed intervals for positions, '
		'and the begin and end options here will fetch the same.')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', metavar='FILE',
		required=True, help='input FastA file')
	req.add_argument('-b', '--begin', metavar='INT', type=int,
		required=True, help='start site to extract')
	req.add_argument('-c', '--contains', metavar='STR', type=str,
		required=True, help='defline contains')
	req.add_argument('-e', '--end', metavar='INT', type=int,
		required=True, help='stop site to extract')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-a', '--all', action='store_false', default=True,
		help='enable all sequence records to be parsed if the -c STR occurs '
		'in >1 defline; default quits after the first record is parsed [off]')
	opt.add_argument('-o', '--outfile', metavar='FILE', default=None,
		help='output FastA file [stdout]')
	opt.add_argument('-p', '--prefix', metavar='STR', type=str, default='',
		help='prefix to add to output defline(s) before the input defline and '
		'extracted sites [None]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	fas = os.path.abspath(os.path.expanduser(opt.infile))
	if opt.begin < 1 or opt.end < 1:
		sys.stderr.write('\tERROR: begin and end positions are 1-based\n')
		sys.exit(1)

	o = []
	with open(fas) as ifh:
		for rec in SeqIO.parse(ifh, 'fasta'):
			if opt.contains in rec.id:
				if opt.begin < opt.end:
					o.append('>{}{} {}:{}'.format(
						opt.prefix, rec.id,opt.begin, opt.end))
					o.append(rec.seq[opt.begin-1:opt.end])
				elif opt.begin > opt.end:
					o.append('>{}{} {}:{}(rc)'.format(
						opt.prefix, rec.id, opt.end, opt.begin))
					o.append(rec.seq[opt.end-1:opt.begin].reverse_complement())
				elif opt.begin == opt.end:
					o.append('>{}{} {}:{}'.format(
						opt.prefix, rec.id, opt.begin, opt.end))
					o.append(rec.seq[opt.begin-1])
				if opt.all:
					break

	if len(o) == 0:
		sys.stderr.write('\t{} absent in deflines\n'.format(opt.contains))
	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	for ln in o:
		ofh.write('{}\n'.format(''.join(ln)))

if __name__ == '__main__':
	main()
