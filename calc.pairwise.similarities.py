#!/usr/bin/env python


from argparse import ArgumentParser
from itertools import combinations
from glob import glob
from re import search
from shutil import rmtree
from tempfile import mkdtemp
import os

def parseArgs():
	parser = ArgumentParser(description='Computes pairwise percent '
		'similarities by performing Needleman-Wunsch global alignments. '
		'Requires the needle binary from The European Molecular Biology '
		'Open Software Suite (EMBOSS).', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--inpath', required=True, metavar='PATH',
		help='path containing input files')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-e', '--ext', type=str, metavar='STR', default='fasta',
		help='file extension to find all nucleotide FastA files within the '
		'specified inpath [fasta]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', metavar='FILE',
		default='Pairwise.Similarities.tab',
		help='output file [./Pairwise.Similarities.tab]')
	return parser.parse_args()

def main():
	opts = parseArgs()
	inpath  = os.path.abspath(opts.inpath)
	ext     = opts.ext
	tmp     = mkdtemp()
	outfile = os.path.abspath(opts.outfile)

	# Get sorted unique list of tuples
	files = glob(os.path.join(inpath, '*.{}'.format(ext)))
	uniq = set()
	for i in combinations(files, 2):
		if i[0] != i[1]:
			uniq.add(i)
	sort_uniq_filepairs = sorted(uniq, key=lambda tupl:tupl[0])

	with open(outfile, 'w') as ofh:
		for i, j in sort_uniq_filepairs:
			b1 = os.path.basename(i).rstrip('.' + ext)
			b2 = os.path.basename(j).rstrip('.' + ext)
			tmpfile = os.path.join(tmp, b1 + ',' + b2 + '.needle')
			
			# Run global alignments
			os.system('needle -asequence {} -snucleotide1 -bsequence '
			'{} -snucleotide2 -gapopen 10.0 -gapextend 0.5 -outfile '
			'{} -aformat3 markx2 -awidth3 60 -auto -error'.format(
				i, j, tmpfile))
			
			# Parse percent similarity values (Version: EMBOSS:6.6.0.0)
			dat = open(tmpfile).readlines()
			f = search(r'#\s+Similarity:\s+\d+\/\d+\s+\(\d+.\d+\%\)\n',
				''.join(dat))
			similarity = f.group(0).split('(')[-1].rstrip('%)\n')
			ofh.write('{}\t{}\t{}\n'.format(b1, b2, similarity))
	rmtree(tmp)

if __name__ == '__main__':
	main()