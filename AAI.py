#!/usr/bin/env python


import os
from argparse import ArgumentParser
from glob import glob
from numpy import mean, std
from subprocess import Popen
from sys import exit
from time import strftime

def parseArgs():
	parser = ArgumentParser(description='Computes the average amino acid '
			'identity (AAI) between two protein sets', add_help=False)
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
	opt.add_argument('-l', '--length', type=int, metavar='INT',
		default=0, help='minimum alignment character length (sum of all '
		'aligned segments and all gaps) [0]')
	opt.add_argument('-o', '--outpath', metavar='PATH',
		default=None, help='output directory [./AAI--<date>_<time>]')
	opt.add_argument('-s', '--bitscore', type=float, metavar='FLOAT',
		default=0.0, help='minimum alignment Bit score [0.0]')
	opt.add_argument('--refilter', default=False, action='store_true',
		help='skip blast system commands and re-filter previously generated '
		'blast output from this script with different cutoff parameters, '
		'which overwrites loci and stats output')
	return parser.parse_args()

def main():
	opts = parseArgs()
	set1 = os.path.abspath(opts.set1)
	set2 = os.path.abspath(opts.set2)
	if opts.outpath is not None:
		outpath = os.path.abspath(opts.outpath)
	else:
		autocreated_dirname = 'AAI--' + strftime('%d%b%Y_%-I:%M%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)
	if not os.path.exists(outpath):
		os.mkdir(outpath)

	if not opts.refilter:
		# Execute bidirectional blast
		for s1, s2 in [(set1, set2), (set2, set1)]:  #s1 as subj; s2 as query
			s1b = os.path.basename(s1).rsplit('.')[0]
			s2b = os.path.basename(s2).rsplit('.')[0]
			c1 = ['makeblastdb', '-in', s1, '-input_type', 'fasta',
				'-out', os.path.join(outpath, s1b), '-dbtype', 'prot']
			c2 = ['blastp', '-db', os.path.join(outpath, s1b), '-max_hsps', '1',
				'-query', s2, '-max_target_seqs', '1', '-num_threads', opts.cpus,
				'-out', os.path.join(outpath, 'blast.'+s1b+','+s2b+'.tab'),
				'-outfmt', '6 qseqid sseqid pident length bitscore qcovhsp']
			for cmd in [c1, c2]:
				with open(os.devnull, 'wb') as dump:
					return_code = Popen(cmd, shell=False, stdout=dump)
					if return_code.wait() != 0:
						exit('ERROR: failed sys call\n{}'.format(' '.join(cmd)))
		junk = []
		for ext in ['*.[np]hr', '*.[np]in', '*.[np]sq']:
			junk.extend(list(glob(os.path.join(outpath, ext))))
		for j in junk:
			os.remove(j)

	# Parse blast output
	aai1, d1 = [], {}
	tot1, filt1 = 0, 0
	with open(os.path.join(outpath, 'blast.'+s2b+','+s1b+'.tab')) as dat:
		for line in dat:
			tot1 += 1
			l = line.split('\t')
			if float(l[2])>=opts.identity and int(l[3])>=opts.length \
			and float(l[4])>=opts.bitscore and float(l[5])>=opts.fraction:
				d1[(str(l[0]), str(l[1]))] = (l[2], l[3], l[4])
				filt1 += 1
				aai1.append(float(l[2]))
	tot2, filt2 = 0, 0
	aai2, twoway = [], []
	with open(os.path.join(outpath, 'blast.'+s1b+','+s2b+'.tab')) as dat, \
	open(os.path.join(outpath, 'aai.'+s2b+','+s1b+'.loci.tab'), 'w') as loci:
		for line in dat:
			tot2 += 1
			l = line.split('\t')
			if float(l[2])>=opts.identity and int(l[3])>=opts.length \
			and float(l[4])>=opts.bitscore and float(l[5])>=opts.fraction:
				filt2 += 1
				aai2.append(float(l[2]))
				if (str(l[1]), str(l[0])) in d1:
					loci.write(l[0]+'\t'+l[1]+'\n')
					twoway.extend([float(l[2]),
									float(d1[str(l[1]), str(l[0])][0])])
	# Write AAI data
	with open(os.path.join(outpath, 'aai.'+s2b+','+s1b+'.stats.tab'), 'w') as aai:
		aai.write('Sample(s)\tFiltered\tAAI\tStandard_Deviation\n')
		aai.write('{},{}\t{}/{}\t{}%\t{}%\n'.format(s2b, s1b, len(twoway),
			tot1+tot2, round(mean(twoway), 3), round(std(twoway), 3)))
		aai.write('{}\t{}/{}\t{}%\t{}%\n'.format(s2b, filt1, tot1,
			round(mean(aai1), 3), round(std(aai1), 3)))
		aai.write('{}\t{}/{}\t{}%\t{}%\n'.format(s1b, filt2, tot2,
			round(mean(aai2), 3), round(std(aai2), 3)))

if __name__ == '__main__':
	main()
