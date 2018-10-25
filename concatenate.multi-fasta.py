#!/usr/bin/env python


import re
import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(add_help=False,
		description='Combines all sequence records in a multi-FastA file '
		'into a single sequence record FastA file.')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input multi-FastA file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-c', '--character', type=str, default='N',
		metavar='STR', help='character to use for insertions between '
		'sequence records [N]')
	opt.add_argument('-d', '--defline', type=str, default=None,
		metavar='STR', help='output sequence record name; if including '
		'whitespace flank with apostrophes [basename infile]')
	opt.add_argument('-o', '--outfile', metavar='FILE', default=None,
		help='output FastA file [stdout]')
	opt.add_argument('-r', '--remove', metavar='\'STR\'', type=str,
		default=None,
		help='prior to combining the individual sequences records together, '
		'remove sites containing at least one specified character(s) '
		'such as gaps or ambiguities; comma-separate for more than one '
		'character and flank with apostrophes if necessary to santize special '
		'characters such as \'N,n,X,x,-\' [none]')
	opt.add_argument('-s', '--gap-size', type=int, default=1000, metavar='INT',
		help='number of characters to insert between sequence records [1000]')
	return parser.parse_args()

def make_linewrapped_fasta(text, defline, width=80):
	''' given a text string, a newline is added every width chars;
	faster than textwrap fill for FastA '''
	# import textwrap
	# wrapper = textwrap.TextWrapper()
	# wrapper.break_long_words = True
	# wrapper.width = width
	# wrapper.replace_whitespace = False
	# wrapper.expand_tabs = False
	# wrapper.break_on_hyphens = False
	# return '>{}\n{}\n'.format(defline, wrapper.fill(seq))	
	split_text = [text[i:i + width] for i in range(0, len(text), width)]
	return '>{}\n{}\n'.format(defline, '\n'.join(split_text))

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.infile))
	seqs = []
	for rec in SeqIO.parse(ifh, 'fasta'):
		seqs.append(str(rec.seq))

	# Optional site filtering
	if opt.remove is not None:
		regex = re.compile(r'{}'.format(
			opt.remove.strip('\'').replace(',', '|')))
		seqs = [re.sub(regex, '', s) for s in seqs]

	# Stitch sequences with character filler
	seq = (opt.character * opt.gap_size).join(seqs)
	if opt.defline is not None:
		defline = opt.defline.strip('\'')
	else:
		if '.' in os.path.basename(ifh):
			defline = '.'.join(os.path.basename(ifh).split('.')[:-1])
		else:
			defline = os.path.basename(ifh)
	out_fasta = make_linewrapped_fasta(seq, defline)

	# Write output file
	if opt.outfile is not None:
		ofh = open(os.path.abspath(os.path.expanduser(opt.outfile)), 'w')
	else:
		ofh = sys.stdout
	ofh.write(out_fasta)

if __name__ == '__main__':
	main()
