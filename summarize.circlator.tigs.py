#!/usr/bin/env python


import os
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqUtils import GC

def parseArgs():
	parser = ArgumentParser(description='Summarizes circlator '
		'output according to each sequence\'s length, GC content, and '
		'optionally its circularized status when the logfile is provided.',
		add_help=False, epilog='Circlator DOI: 10.1186/s13059-015-0849-0')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', metavar='FILE',
		required=True, help='input FastA file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-l', '--circlator-log', metavar='FILE',
		help='corresponding tab-delimited \'04.merge.circularise.log\' file '
		'from Circlator; required for reporting circularized status')
	opt.add_argument('-o', '--outfile', required=False, metavar='FILE',
		default=None, help='summary output [stdout]')
	opt.add_argument('-s', '--separate', required=False, action='store_true',
		default=False, help='switch to copy each sequence record (contig or '
		'unitig) into its own separate file in the cwd [OFF]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.infile))
	log = os.path.abspath(os.path.expanduser(opt.circlator_log))

	d = {}
	if log is not None:
		with open(log) as cl:
			for ln in cl:
				l = ln.rstrip('\n').split('\t')
				d[l[1]] = l[5]

	with open(ifh) as input_fasta:
		for rec in SeqIO.parse(input_fasta, 'fasta'):
			seq_id = str(rec.id)
			if opt.separate:
				# Create individual files for each seq rec (contig or unitig)
				SeqIO.write(rec, seq_id + '.fa', 'fasta')

			# Gather summary data
			if log is not None:
				if seq_id in d:
					v = d[seq_id]
					if v == '1':
						status = 'circularized'
					else:
						status = 'not_circularized'
				else:
					status = 'unitig_not_in_logfile'
			else:
				status = ''
			seq_len = str(len(rec.seq))
			GC_content = float(GC(rec.seq))
			data = '{}\t{}\t{:.2f}\t{}'.format(seq_id, seq_len, GC_content,
				status)

			# Output
			if opt.outfile is not None:
				ofh = os.path.abspath(os.path.expanduser(opt.outfile))
				with open(ofh, 'wa') as o:
					o.write(data + '\n')
			else:
				print(data)

if __name__ == '__main__':
	main()