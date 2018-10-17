#!/usr/bin/env python


from argparse import ArgumentParser


def parseArgs():
	parser = ArgumentParser(description='Filters BLAST output format 6 '
		'(tab-delimited) for best hit based on bitscore. '
		'Handles additional data columns following bitscore.',
		add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input file in NCBI\'s BLAST -outfmt 6 format')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-c', '--column', type=int, metavar='{1,2}', choices=[1, 2],
		default=1, help='report best hit per query label (1st column; \'1\') '
		'or target (2nd column; \'2\') [1]')
	opt.add_argument('-o', '--outfile', metavar='FILE',
		default=None, help='output file [stdout]')
	opt.add_argument('-s', '--bitscore', type=float, metavar='FLOAT',
		default=0.0, help='minimum alignment Bit score [0.0]')
	return parser.parse_args()

def main():
	opts = parseArgs()
	
	# Identify unique query labels and will not assume sorted file
	with open(opts.infile, 'r') as ifh:
		query_labels = set()
		for ln in ifh:
			query_labels.add(ln.split('\t')[opts.column-1])
	num_cols = len(ln.split('\t'))

	# Filter best hits
	best = {k: ['0']*num_cols for k in query_labels}
	with open(opts.infile, 'r') as ifh:
		for ln in ifh:
			data = ln.rstrip('\n').split('\t')
			if float(data[11]) > opts.bitscore and \
			float(data[11]) > float(best[data[opts.column-1]][11]):
				best[data[opts.column-1]] = data

	# Write output
	if opts.outfile is None:
		outfile = open(os.path.abspath(os.path.expanduser(opts.outfile)), 'w')
	else:
		outfile = opts.outfile
	for k,v in sorted(best.items()):
		if v != ['0']*num_cols:
			o.write('\t'.join(v) + '\n')

if __name__ == '__main__':
	main()
