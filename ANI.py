#!/usr/bin/env python


import os
from argparse import ArgumentParser
from collections import deque
from itertools import islice
from numpy import mean, std
from shutil import rmtree
from subprocess import Popen
from sys import exit
from tempfile import mkdtemp
from time import strftime
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Computes the average nucleotide '
	'identity (ANI) between two nucleic acid sequence sets', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-1', '--set1', required=True, metavar='FILE',
		help='first input FastA sequence file')
	req.add_argument('-2', '--set2', required=True, metavar='FILE',
		help='second input FastA sequence file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--cpus', type=str, metavar='INT',
		default='1', help='number of CPUs [1]')
	opt.add_argument('-f', '--fraction', type=float, metavar='FLOAT',
		default=70.0, help='minimum alignment length percentage [70.0]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-i', '--identity', type=float, metavar='FLOAT',
		default=30.0, help='minimum percent identity [30.0]')
	opt.add_argument('-j', '--fill', type=str, metavar='CHAR', default='', 
		help='character to add to end of fragments that are less than the '
		'window size to force them to be the same length; requires the -k '
		'switch as well [None]')
	opt.add_argument('-k', '--keep-small-frags', default=False,
		action='store_true', help='keep remaining nucleotide sequences less '
		' than the window size')
	opt.add_argument('-l', '--length', type=int, metavar='INT',
		default=0, help='minimum alignment character length (sum of all '
		'aligned segments and all gaps) [0]')
	opt.add_argument('-o', '--outpath', metavar='PATH',
		default=None, help='output directory [./ANI--<date>_<time>]')
	opt.add_argument('-s', '--step-size', type=int, metavar='INT',
		default=200, help='nucleotide step size during sequence '
		'fragmentation prior to BLAST alignments; to turn off steps to '
		'speed up, set to same length as the fragment size [200]')
	opt.add_argument('-w', '--fragment-size', type=int, metavar='INT',
		default=1000, help='fragment lengths to slice nucleotide sequence '
		'sets into prior to BLAST alignments [1000]')
	opt.add_argument('--refilter', default=False, action='store_true',
		help='skip blast system commands and re-filter previously generated '
		'blast output from this script with different cutoff parameters, '
		'which overwrites loci and stats output')
	return parser.parse_args()

def fragment(seq, win, step, fill):
	iters = iter(seq)
	q = deque(islice(iters, win), maxlen=win)
	q.extend(fill for _ in range(win-len(q)))
	while True:
		yield q
		q.append(next(iters))
		q.extend(next(iters, fill) for _ in range(step-1))

def main():
	opts = parseArgs()
	set1 = os.path.abspath(opts.set1)
	set2 = os.path.abspath(opts.set2)
	b1   = os.path.basename(set1).rsplit('.')[0]
	b2   = os.path.basename(set2).rsplit('.')[0]
	if opts.outpath is not None:
		outpath = os.path.abspath(opts.outpath)
	else:
		autocreated_dirname = 'ANI--' + strftime('%d%b%Y_%-I:%M%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)
	if not os.path.exists(outpath):
		os.mkdir(outpath)

	if not opts.refilter:
		# Fragment input sequences
		tmp  = mkdtemp()
		ifh1 = os.path.join(tmp, b1 + '.frags')
		ifh2 = os.path.join(tmp, b2 + '.frags')
		for f, s in [(ifh1, set1), (ifh2, set2)]:
			with open(f, 'w') as ofh:
				for rec in SeqIO.parse(s, 'fasta'):
					i = 1
					for frag in fragment(rec.seq, opts.fragment_size,
					opts.step_size, opts.fill):
						frag_seq = ''.join(list(frag))
						if opts.keep_small_frags or \
						len(frag_seq) == opts.fragment_size:
							ofh.write('>{}\n{}\n'.format(
								str(i) + '__' + rec.id, frag_seq))
							i += 1

		# Execute bidirectional blast
		for s1, s2, sb1, sb2 in [(ifh1, ifh2, b1, b2), (ifh2, ifh1, b2, b1)]:
			c1 = ['makeblastdb', '-in', s1, '-out', s1, '-dbtype', 'nucl']
			c2 = ['blastn', '-db', s1, '-query', s2, '-dust', 'no',
				'-max_hsps', '1', '-max_target_seqs', '1',
				'-num_threads', opts.cpus, '-task', 'blastn',
				'-outfmt', '6 qseqid sseqid pident length bitscore qcovhsp',
				'-out', os.path.join(outpath, 'blast.'+sb1+','+sb2+'.tab')]
			for cmd in [c1, c2]:
				with open(os.devnull, 'wb') as dump:
					return_code = Popen(cmd, shell=False, stdout=dump)
					if return_code.wait() != 0:
						exit('ERROR: failed sys call\n{}'.format(' '.join(cmd)))
		rmtree(tmp)

	# Parse blast output
	ani1, d1 = [0], {}
	tot1, filt1, len1 = 0, 0, 0 #counts for set1 as query
	with open(os.path.join(outpath, 'blast.'+b2+','+b1+'.tab')) as dat:
		for line in dat:
			tot1 += 1
			l = line.split('\t')
			if float(l[2])>=opts.identity and int(l[3])>=opts.length \
			and float(l[5])>=opts.fraction:
				d1[(str(l[0]), str(l[1]))] = (l[2], l[3], l[4])
				filt1 += 1
				len1  += int(l[3])
				ani1.append(float(l[2]))
	tot2, filt2, len2, twowaylen = 0, 0, 0, 0 #counts for set2 as query
	ani2, twoway = [0], [0]
	with open(os.path.join(outpath, 'blast.'+b1+','+b2+'.tab')) as dat:
		for line in dat:
			tot2 += 1
			l = line.split('\t')
			if float(l[2])>=opts.identity and int(l[3])>=opts.length \
			and float(l[5])>=opts.fraction:
				filt2 += 1
				len2  += int(l[3])
				ani2.append(float(l[2]))
				if (str(l[1]), str(l[0])) in d1:
					twowaylen += int(l[3])
					twoway.extend([float(l[2]),
									float(d1[str(l[1]), str(l[0])][0])])
	# Write ANI data
	with open(os.path.join(outpath, 'ani.'+b1+','+b2+'.stats.tab'), 'w') \
	as ani:
		ani.write('Query_Sample(s)\tFiltered\tANI\tStandard_Deviation\t'
			'Alignment\n')
		ani.write('{},{}\t{}/{}\t{}%\t{}%\t{}\n'.format(b1, b2, len(twoway),
			tot1+tot2, round(mean(twoway), 3), round(std(twoway), 3),
			twowaylen))
		ani.write('{}\t{}/{}\t{}%\t{}%\t{}\n'.format(b1, filt1, tot1,
			round(mean(ani1), 3), round(std(ani1), 3), len1))
		ani.write('{}\t{}/{}\t{}%\t{}%\t{}\n'.format(b2, filt2, tot2,
			round(mean(ani2), 3), round(std(ani2), 3), len2))

if __name__ == '__main__':
	main()