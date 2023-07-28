#!/usr/bin/env python3


import csv
import gzip
import multiprocessing as mp
import os
import shutil
import sys
from argparse import ArgumentParser
from glob import glob
from multiprocessing import cpu_count
from tempfile import mkdtemp
from io import StringIO

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def parseArgs():
	parser = ArgumentParser(description='Creates a sequence length '
		'distribution plot with a statistical summary', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--indir', required=True, metavar='DIR',
		help='input path containing read files')
	opt = parser.add_argument_group('Input Options')
	opt.add_argument('-e', '--extension', metavar='STR', type=str, nargs='+',
		default=['fastq', 'fastq.gz', 'fast5'],
		help='file extension(s) to assess [fastq, fastq.gz, fast5]')
	opt.add_argument('-p', '--prefix', type=str, metavar='STR', default='',
		help='filename prefix required in each sequence file [*]')
	fig = parser.add_argument_group('Figure Options')
	fig.add_argument('-o', '--outfile', metavar='FILE', type=str,
		default='read-lengths.pdf', help='output file [read-lengths.pdf]')
	fig.add_argument('--subtitle', metavar='\"STR\"', type=str,
		default=None, help='figure subtitle; surround with quotes [None]')
	fig.add_argument('--title', metavar='\"STR\"', type=str,
		default=None, help='figure title; surround with quotes [None]')
	fig.add_argument('--bins', type=int, metavar='INT', default=100,
		help='number of length groups for plotting [100]')
	fig.add_argument('--min-length', type=int, metavar='INT', default=1,
		help='minimum read length [1]')
	fig.add_argument('--max-length', type=int, metavar='INT', default=9999999,
		help='maximum read length [9999999]')
	mis = parser.add_argument_group('Misc Options')
	mis.add_argument('-c', '--cpus', type=require_int_nonnegative,
		metavar='INT', default='0', help='number of CPUs [all]')
	mis.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	mis.add_argument('-l', '--lengths-outfile', metavar='FILE', type=str,
		default=None, help='output file listing lengths [None]')
	return parser.parse_args()

def require_int_nonnegative(x):
	try:
		if int(x) < 0 or '.' in str(x):
			sys.stderr.write('ERROR: {} must be a non-negative integer\n'.\
				format(x))
			sys.exit(1)
	except ValueError:
		sys.stderr.write('ERROR: {} must be an integer\n'.format(x))
		sys.exit(1)
	return int(x)

def require_minimum_arg_value(arg_value, arg_name, min_value):
	if arg_value < min_value:
		sys.stderr.write('ERROR: minimum {} must be at least {}\n'.\
			format(arg_name, arg_value))
		sys.exit(1)

def calc_sequence_lengths(args):
	infile, tempdir, min_length, max_length = args
	# Auto-handle gunzip compressed input
	if infile.endswith('.gz'):
		tmp_outfile = os.path.join(tempdir,
			os.path.basename(infile).rstrip('.gz'))
		with gzip.open(infile) as ifh, open(tmp_outfile, 'w') as ofh:
			shutil.copyfileobj(ifh, ofh)
		infile = tmp_outfile
	outfile = os.path.join(tempdir, os.path.basename(infile) + '.len')

	# Count lengths within the sequence file
	i = 0
	if infile.endswith('.fastq'):
		with open(infile) as ifh, open(outfile, 'w') as ofh:
			for _, seq_string, _ in FastqGeneralIterator(ifh):
				sequence_length = len(seq_string)
				if min_length < sequence_length < max_length:
					i += 1
					ofh.write(str(sequence_length) + os.linesep)
	elif infile.endswith('.fast5'):
		import h5py
		# http://docs.h5py.org/en/stable/quick.html
		try:
			fast5_file = h5py.File(infile, 'r')
			names = []
			fast5_file.visit(names.append)
		except (IOError, RuntimeError) as e:
			sys.stderr.write('ERROR: {} appears corrupt\n{}\n'.\
				format(infile, e))
			sys.exit(1)
		keys = {
			'2D template':
			'Analyses/Basecall_1D_%03d/BaseCalled_template/Fastq',
			'2D complement':
			'Analyses/Basecall_2D_%03d/BaseCalled_complement/Fastq',
			'2D twodirections':
			'Analyses/Basecall_2D_%03d/BaseCalled_2D/Fastq',
			'1D template':
			'Analyses/Basecall_1D_000/BaseCalled_template/Fastq',
			'1D complement':
			'Analyses/Basecall_1D_%03d/BaseCalled_complement/Fastq'
		}
		with open(outfile, 'w') as ofh:
			# single read per Fast5
			for k, v in keys.iteritems():
				try:
					_ = fast5_file[v][()]
					# now that first read extraction works, do all of them
					for name in names:
						fastq = fast5_file[v][()]
						record = SeqIO.read(StringIO(fastq), 'fastq')
						sequence_length = len(rec.seq)
						if min_length < sequence_length < max_length:
							i += 1
							ofh.write(str(sequence_length) + os.linesep)
					sys.stderr.write('INFO: detected {} format for {}\n'.\
						format(k, infile))
					break
				except KeyError:
					pass
			# multi-Fast5 format has a different layout
			if i == 0:
				names = [str(s) for s in fast5_file.keys()]
				for k, v in keys.iteritems():
					try:
						_ = fast5_file[names[0]][v]
						# now that first read extraction works, do all of them
						# NOTE: find out how to get an iter instead of lookups
						for name in names:
							record = fast5_file[name][v].value
							record_chunks = record.rstrip('\n').split('\n')
							if len(record_chunks) == 4:
								seq_id, seq, sep, qual = record_chunks
								sequence_length = len(seq)
								if min_length < sequence_length < max_length:
									i += 1
									ofh.write(str(sequence_length) +
										os.linesep)
							else:
								sys.stderr.write('ERROR: unable to parse'
									' FastQ from {} in {}\n'.\
									format(name, infile))
						sys.stderr.write('INFO: detected {} format for {}\n'.\
							format(k, infile))
						break
					except KeyError:
						pass
		fast5_file.close()
	elif infile.endswith('.fasta') or infile.endswith('.fa'):
		with open(infile) as ifh, open(outfile, 'w') as ofh:
			for rec in SeqIO.parse(ifh, 'fasta'):
				sequence_length = len(rec.seq)
				if min_length < sequence_length < max_length:
					i += 1
					ofh.write(str(sequence_length) + os.linesep)
	if i < 1:
		sys.stderr.write('INFO: no sequences between {} and {:,.0f}'
			' found in {}\n'.format(min_length, max_length, infile))

def merge_count_files(count_files, merged_file):
	sequence_lengths = []
	for file in count_files:
		with open(file) as ifh:
			counts = ifh.read().splitlines()
		if len(counts) > 0:
			sequence_lengths.extend(counts)
		else:
			sys.stderr.write('INFO: no sequence lengths found in {}\n'.\
				format(file))

	sequence_lengths.sort(key=int, reverse=True)
	with open(merged_file, 'w') as ofh:
		for val in sequence_lengths:
			ofh.write(val + os.linesep)

	return sequence_lengths

def calculate_stats(merged_lengths):
	stats = {}
	stats['q1'] = int(np.percentile(merged_lengths, 25,
		interpolation='midpoint'))
	stats['median'] = int(np.median(merged_lengths))
	stats['q3'] = int(np.percentile(merged_lengths, 75,
		interpolation='midpoint'))
	stats['mean'] = int(np.mean(merged_lengths))
	stats['stdev'] = int(np.std(merged_lengths))
	stats['gt3kbp'] = int(np.sum(merged_lengths > 3000))
	stats['gt10kbp'] = int(np.sum(merged_lengths > 10000))
	stats['gt20kbp'] = int(np.sum(merged_lengths > 20000))
	stats['total_count'] = len(merged_lengths)
	stats['total_length'] = int(np.sum(merged_lengths))
	stats['largest'] = int(np.amax(merged_lengths))
	total_length = stats['total_length']
	limit = .5 * stats['total_length']
	for length in merged_lengths:
		total_length -= length
		if total_length <= limit:
			stats['n50'] = length
			break
	else:
		stats['n50'] = 0
	return stats

def plot_histogram(lengths, num_bins, stats, outfile, title, subtitle):
	sns.set_style('whitegrid')
	fig, ax = plt.subplots()
	ax.hist(lengths, bins=num_bins)

	# Place stats in upper right of axes coords
	textblock = '\n'.join((
		'Total Length = {:,.0f} bp'.format(stats['total_length']),
		'N50 = {:,.0f} bp'.format(stats['n50']),
		r'$\mu = {:,.0f}$ bp'.format(stats['mean']),
		r'$\sigma = {:,.0f}$ bp'.format(stats['stdev']),
		'Q1 = {:,.0f} bp'.format(stats['q1']),
		r'median = {:,.0f} bp'.format(stats['median']),
		'Q3 = {:,.0f} bp'.format(stats['q3']),
		'largest = {:,.0f} bp'.format(stats['largest']),
		'',
		'Total Reads = {:,.0f}'.format(stats['total_count']),
		'>3 kbp = {:,.0f} ({:.0f}%)'.format(stats['gt3kbp'],
			100. * stats['gt3kbp'] / stats['total_count']),
		'>10 kbp = {:,.0f} ({:.0f}%)'.format(stats['gt10kbp'],
			100. * stats['gt10kbp'] / stats['total_count']),
		'>20 kbp = {:,.0f} ({:.0f}%)'.format(stats['gt20kbp'],
			100. * stats['gt20kbp'] / stats['total_count'])
		))
	legend = dict(boxstyle='round', facecolor='grey', alpha=0.5)
	ax.text(0.6, 0.98, textblock, transform=ax.transAxes, fontsize=8,
		fontdict={'family':'monospace'}, verticalalignment='top', bbox=legend,
		multialignment='right')

	# Add labels
	plt.xlabel('Read Length [bp]', axes=ax)
	plt.ylabel('Quantity [#]', axes=ax)
	if title is not None:
		fig.suptitle(title, fontsize=10)
	if subtitle is not None:
		ax.set_title(subtitle, fontsize=8)

	# Save the plot
	plt.savefig(outfile, pad_inches=0.5, dpi=300)

def main():
	opt = parseArgs()
	min_length, max_length = opt.min_length, opt.max_length
	require_minimum_arg_value(min_length, 'length', 1)
	if min_length > max_length:
		sys.stderr.write('ERROR: --min-length cannot exceed --max-length\n')
		sys.exit(1)
	bins = opt.bins
	require_minimum_arg_value(bins, 'bins', 2)
	indir = os.path.realpath(os.path.expanduser(opt.indir))
	outfile = os.path.realpath(os.path.expanduser(opt.outfile))
	extensions, prefix = opt.extension, opt.prefix
	if opt.cpus < 1:
		cpus = cpu_count()
	else:
		cpus = opt.cpus

	# Locate sequence files
	sequence_files = []
	for ext in extensions:
		sequence_files.extend(glob(os.path.join(indir, '{}*.{}'.\
			format(prefix, ext))))
	if len(sequence_files) < 1:
		sys.stderr.write('ERROR: no sequence files found\n')
		sys.exit(1)
	else:
		sys.stderr.write('INFO: found {} sequence files...\n'.\
			format(len(sequence_files)))

	# Calculate lengths of each sequence within all sequence files
	tmp = mkdtemp()
	pool = mp.Pool(processes=cpus)
	args = [(file, tmp, min_length, max_length) for file in sequence_files]
	pool.map_async(calc_sequence_lengths, args).get(timeout=99999)

	# Summarize all length files
	length_files = glob(os.path.join(tmp, '*.len'))
	if len(length_files) < 1:
		sys.stderr.write('ERROR: no read length count files found\n')
		shutil.rmtree(tmp)
		sys.exit(1)
	merged_file = os.path.join(tmp, 'merged_file.txt')
	merged_lengths = merge_count_files(length_files, merged_file)
	if len(merged_lengths) < 3:
		sys.stderr.write('INFO: not enough sequence lengths to plot'
			' and calculate length statistics\n')
		shutil.rmtree(tmp)
		sys.exit(0)
	merged_lengths = np.array(merged_lengths, dtype=int)
	stats = calculate_stats(merged_lengths)
	plot_histogram(merged_lengths, bins, stats, outfile, opt.title,
		opt.subtitle)

	# Cleanup
	if opt.lengths_outfile is not None:
		shutil.copy(merged_file,
			os.path.realpath(os.path.expanduser(opt.lengths_outfile)))
	shutil.rmtree(tmp)

if __name__ == '__main__':
	main()
