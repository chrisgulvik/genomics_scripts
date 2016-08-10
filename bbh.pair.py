#!/usr/bin/env python


import os
from argparse import ArgumentParser
from decimal import Decimal, ROUND_HALF_UP
from glob import glob
from subprocess import Popen
from sys import exit
from time import strftime

def parseArgs():
	parser = ArgumentParser(description='Lists and quantifies bidirectional best hits (BBH) of nucleotide or protein sets. Alignment filters are properly applied in both directions. Requires NCBI BLAST 2.2.29+ or newer', add_help=False,
	epilog='Output is saved as <outpath>/blast.<set1-ext>,<set2-ext>.tab, <outpath>/blast.<set2-ext>,<set1-ext>.tab, and <outpath>/bbh.<set1-ext>,<set2-ext>.{discard,filt,stats}.tab. If iteratively testing parameters on the same dataset (e.g., percent identity effects) with the --refilter option, rename bbh.*.tab files after each iteration or files will be overwritten.')
	req = parser.add_argument_group('Required')
	req.add_argument('-1', '--set1', metavar='FILE', required=True, help='first input FastA sequence file')
	req.add_argument('-2', '--set2', metavar='FILE', required=True, help='second input FastA sequence file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-b', '--bitscore', metavar='FLOAT', type=float, default=0, help='minimum alignment Bit score [0]')
	opt.add_argument('-c', '--cpus', metavar='INT', type=str, default='1', help='number of CPUs [1]')
	opt.add_argument('-e', '--ext', metavar='\'STR\'', type=str, default=None, help='file extension to strip from each set for generating paired output filenames; optionally flanked by single quotes [\'.<ext>\']')
	opt.add_argument('-f', '--fraction', metavar='FLOAT', type=float, default=30, help='minimum alignment length percentage [30]')
	opt.add_argument('-h', '--help', action='help', help='show this help message and exit')
	opt.add_argument('-i', '--identity', metavar='FLOAT', type=float, default=80, help='minimum percent identity [80]')
	opt.add_argument('-l', '--length', metavar='INT', type=int, default=0, help='minimum alignment character length (sum of all aligned segments and all gaps) [0]')
	opt.add_argument('-m', '--moleculetype', default='prot', choices=['nucl', 'prot'], help='input molecule type [prot]')
	opt.add_argument('-o', '--outpath', metavar='PATH', default=None, help='output directory [./BBH.pair--<date>_<time>]')
	opt.add_argument('--refilter', default=False, action='store_true', help='skip blast system commands and re-filter previously generated blast output from this script with different cutoff parameters')
	return parser.parse_args()

def main():
	opts = parseArgs()
	set1 = os.path.abspath(opts.set1)
	set2 = os.path.abspath(opts.set2)
	if opts.outpath is not None:
		outpath = os.path.abspath(opts.outpath)
		if not os.path.exists(outpath):
			os.mkdir(outpath)
	else:
		autocreated_dirname = 'BBH.pair--'+strftime('%d%b%Y_%-I:%M%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)
		os.mkdir(outpath)
	if opts.ext is not None:
		b1 = os.path.basename(set1).rstrip(opts.ext.lstrip('\'').rstrip('\''))
		b2 = os.path.basename(set2).rstrip(opts.ext.lstrip('\'').rstrip('\''))
	else:
		b1 = '.'.join(os.path.basename(set1).split('.')[:-1])
		b2 = '.'.join(os.path.basename(set2).split('.')[:-1])

	if not opts.refilter:
		# Execute bidirectional blast
		mol_type = opts.moleculetype 
		if mol_type == 'nucl':
			blast = 'blastn'
		elif mol_type == 'prot':
			blast = 'blastp'
		for s1, s2 in [(set1, set2), (set2, set1)]:  #s1 as subj; s2 as query
			s1b = os.path.basename(s1).rsplit('.')[0]
			s2b = os.path.basename(s2).rsplit('.')[0]
			cmd_makedb= ['makeblastdb', '-in', s1, '-input_type', 'fasta',
						'-out', os.path.join(outpath, s1b), '-dbtype', mol_type]
			cmd_blast = [blast, '-db', os.path.join(outpath, s1b), '-max_hsps', '1',
						'-query', s2, '-max_target_seqs', '1', '-num_threads', opts.cpus,
						'-out', os.path.join(outpath, 'blast.'+s1b+','+s2b+'.tab'),
						'-outfmt', '6 qseqid sseqid pident length bitscore qcovs']
			for cmd in [cmd_makedb, cmd_blast]:
				with open(os.devnull, 'wb') as dump:
					return_code = Popen(cmd, shell=False, stdout=dump)
					if return_code.wait() != 0:
						exit('ERROR: failed system call\n{}\nBLAST 2.2.29+ or newer required'.format(' '.join(cmd)))
		junk = []
		for ext in ['*.[np]hr', '*.[np]in', '*.[np]sq']:
			junk.extend(list(glob(os.path.join(outpath, ext))))
		for j in junk:
			os.remove(j)

	# Parse the blast output
	d = {}
	hit_s2 = 0
	with open(os.path.join(outpath, 'blast.'+b2+','+b1+'.tab')) as dat, \
	open(os.path.join(outpath, 'bbh.'+b2+','+b1+'.discard.tab'), 'w') as discfh:
		for line in dat:
			hit_s2 += 1
			l = line.split('\t')
			if float(l[2]) >= opts.identity and int(l[3]) >= opts.length and float(l[4]) >= opts.bitscore and float(l[5]) >= opts.fraction:
				d[(str(l[1]), str(l[0]))] = (float(l[2]), int(l[3]), float(l[4]), float(l[5]))
			else:
				discfh.write(line)
	filt_s2 = len(d)
	filt_s1 = hit_s1 = bbh = 0
	with open(os.path.join(outpath, 'blast.'+b1+','+b2+'.tab')) as dat, \
	open(os.path.join(outpath, 'bbh.'+b1+','+b2+'.discard.tab'), 'w') as discfh, \
	open(os.path.join(outpath, 'bbh.'+b1+','+b2+'.filt.tab'), 'w') as bbhfh:
		for line in dat:
			hit_s1 += 1
			l = line.split('\t')
			if float(l[2]) >= opts.identity and int(l[3]) >= opts.length and float(l[4]) >= opts.bitscore and float(l[5]) >= opts.fraction:
				filt_s1 += 1
				if (str(l[0]), str(l[1])) in d:
					bbhfh.write(line)
					bbh += 1
			else:
				discfh.write(line)

	# Report results
	pbbh = round(float(bbh)/((hit_s1+hit_s2)/2)*100, 3)
	with open(os.path.join(outpath, 'bbh.'+b1+','+b2+'.stats.tab'), 'w') as stats:
		stats.write('{},{}\tBBHs:\t{}\t({}%)\n'.format(b1, b2, bbh, pbbh))
		stats.write('{}\ttotal hits as subject:\t{}\n'.format(b1, hit_s1))
		stats.write('{}\tdiscarded hits as subject:\t{}\t({}%)\n'.format(b1,
			(hit_s2-filt_s2), round((hit_s2-filt_s2)/float(hit_s2)*100, 3) ))
		stats.write('{}\ttotal hits as subject:\t{}\n'.format(b2, hit_s2))
		stats.write('{}\tdiscarded hits as subject:\t{}\t({}%)\n'.format(b2,
			(hit_s1-filt_s1), round((hit_s1-filt_s1)/float(hit_s1)*100, 3) ))
	print '{} BBHs ({}% passed filters)'.format(bbh, pbbh)

if __name__ == '__main__':
	main()
