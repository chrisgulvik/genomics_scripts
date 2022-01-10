#!/usr/bin/env python


import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess as sp
import sys
from multiprocessing import cpu_count
from tempfile import mkdtemp

from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import vcf


def parseArgs():
	parser = argparse.ArgumentParser(description='creates a scatter'
		' histogram plot of reads mapped to a reference', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', type=str, metavar='FILE', nargs='+',
		help='input FastQ file(s) optionally gunzip compressed; two inputs'
		' uses paired-end mapping')
	req.add_argument('-r', '--reference', metavar='FILE', required=True,
		type=str,
		help='FastA reference for reads to map to, optionally gunzip'
		' compressed')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--cpus', type=require_int_nonnegative,
		metavar='INT', default='0', help='number of CPUs [all]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-m', '--map-options', type=str,
		default='\'--no-unal --no-mixed --no-discordant -X 1000\'',
		metavar='\'STR\'',
		help='options passed to bowtie2 aligner; surround with single quotes'
		' [\'--no-unal --no-mixed --no-discordant -X 1000\']')
	opt.add_argument('-o', '--outfile', default=None, metavar='FILE',
		help='output plot filename [map-to-ref.pdf]')
	opt.add_argument('-x', '--xlsx', default=None, metavar='FILE',
		help='output XLSX filename [None]')
	return parser.parse_args()

def require_nargs_range(x, val_min, val_max):
	if val_min > len(x) > val_max:
		sys.stderr.write('ERROR: {} must not be lower than {} or higher than'
			' {}'.format(x, val_min, val_max))
		sys.exit(1)
	return str(x)

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

def require_dependency(dep):
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, dep)) and \
		not os.path.isdir(os.path.join(path, dep)):
			return True
	sys.stderr.write('ERROR: {} unavailable; not in $PATH\n'.format(dep))
	sys.exit(1)

def verify_file_exists_and_nonempty(infile):
	if os.path.exists(infile):
		if os.stat(infile).st_size == 0:
			sys.stderr.write('ERROR: {} file empty\n'.format(infile))
			sys.exit(1)
	else:
		sys.stderr.write('ERROR: {} file absent\n'.format(infile))
		sys.exit(1)

def system_call(cmd):
	process = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
	_, err = process.communicate()
	if process.returncode != 0:
		sys.stderr.write(err)
		sys.stderr.write('ERROR: failed system call: {}\n'.\
			format(' '.join(cmd)))
		sys.exit(1)

def cleanup_and_evaluate_mfasta(infile, outfile):
	'''
	avoid funky input defline headaches with alignment software by renaming
	and report contig lengths
	'''
	records, sequence_lengths = [], []
	i = 1
	for rec in SeqIO.parse(infile, 'fasta'):
		if float(GC(rec.seq)) == 0:
			sys.stderr.write('ERROR: {} appears to lack nucleotides\n'.\
				format(infile))
			sys.exit(1)
		sequence_lengths.append(len(rec.seq))
		# records.append(SeqRecord.SeqRecord(id='ctg_{}'.format(i),
		# 	seq=rec.seq, description=''))
		records.append(SeqRecord.SeqRecord(
			id=str(i), #plain integers for easy sorting VCF by chrom name
			seq=rec.seq, description=''))
		i += 1
	SeqIO.write(records, outfile, 'fasta')
	return sequence_lengths

def read_vcf_with_pyvcf(infile):
	'''
	takes VCF file as input, returns list of tuples containing:
	chrom, pos, depth, percent identity
	'''
	# strict allows you to parse files with spaces in the sample names
	vcf_reader = vcf.Reader(filename=infile, strict_whitespace=True)
	mapped_sites = [] #list of tuples where each tuple contains:
					# (chrom, pos, depth, %identity)
	for rec in vcf_reader:
		chrom, pos, alleles = rec.CHROM, rec.POS, rec.INFO['AD']
		total_depth = sum(alleles)
		if total_depth > 0:
			ref_depth = float(alleles[0])
			percent_identity = (ref_depth / total_depth) * 100
			print('{}\t{}\t{}\t{}'.format(chrom, pos, total_depth, percent_identity))
			mapped_sites.append(chrom, pos, total_depth, percent_identity)
	return mapped_sites

def read_vcf(infile, check_dupes=True):
	'''
	takes VCF file as input, returns pandas dataframe
	'''
	# load in VCF file
	skip_rows = 0
	with open(infile) as ifh:
		for ln in ifh:
			if ln.startswith('##'):
				skip_rows += 1
			else:
				break
	df = pd.read_csv(infile, sep='\t', skiprows=skip_rows, header=0,
		dtype={'#CHROM': int, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
		'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str},).\
		rename(columns={'#CHROM': 'CHROM'})
	if len(df) < 1:
		sys.stderr.write('ERROR: no alignment data present in VCF.'
			' Verify reads align to the input reference. Files are'
			' within {}\n'.format(os.path.dirname(infile)))
		sys.exit(1)
	# remove unused columns
	df.drop(columns=['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT'])

	if check_dupes:
		dupes = df[df.duplicated(subset=['CHROM','POS'], keep=False)]
		if len(dupes) > 0:
			sys.stderr.write('ERROR: duplicate positions found, '
				'perhaps due to SNPs and InDels. Only one position '
				'permitted, so either skip InDels or merge prior to trying '
				'to fill unmapped sites.\n')
			dupes = dupes[['CHROM', 'POS']]
			sys.stderr.write( dupes.to_string(index=False) + '\n')
			sys.exit(1)
	return df

def calc_stats_from_frequency_distribution_table(freqs, vals):
	values, freqs = np.array(vals), np.array(freqs)
	arg_sorted = np.argsort(values)
	values = values[arg_sorted]
	freqs = freqs[arg_sorted]
	count = freqs.sum()
	fx = values * freqs
	mean = fx.sum() / count
	variance = ((freqs * values**2).sum() / count) - mean**2
	variance = count / (count - 1) * variance
	std = np.sqrt(variance)
	minimum, maximum = np.min(values), np.max(values)
	cumcount = np.cumsum(freqs)
	Q1 = values[np.searchsorted(cumcount, 0.25*count)]
	Q2 = values[np.searchsorted(cumcount, 0.50*count)]
	Q3 = values[np.searchsorted(cumcount, 0.75*count)]
	idx = ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']
	df = pd.Series([count, mean, std, minimum, Q1, Q2, Q3, maximum], index=idx)
	return df

def convert_vcf_to_stats(vcf_df):
	'''
	calculates depth and percent identity from vcf dataframe
	'''
	# make Allelic Depths "AD" its own column
	sample_df = vcf_df.loc[:, vcf_df.columns.str.endswith('sample.sorted.sam')]
	# sample_df = df.filter(regex='sample.sorted.sam$', axis=1)
	if len(sample_df.columns) != 1:
		sys.stderr.write('ERROR: VCF file requires the mapped sample name to'
		' end with sample.sorted.sam\n')
		sys.exit(1)
	pl_ad = sample_df.iloc[:, 0].str.split(pat=':', expand=True)
	if len(pl_ad.columns) != 2:
		sys.stderr.write('ERROR: expected sample column in VCF to be split'
			' by : into two columns in pandas. The sample column should have'
			' ended with sample.sorted.sam (filename), and contained the PL'
			' and AD tags but only the AD second is used.\n')
		sys.exit(1)
	pl_ad.columns = ['PL', 'AD'] #AD,"Allelic depths (high-quality bases)"
	# calculate read depth per site
	allele_depths = pl_ad['AD'].str.split(pat=',', expand=True).\
		fillna(0).astype(int)
	total_depth_per_site = allele_depths.sum(axis=1)
	# calculate percent identity mapped to each reference site 
	ref_depths = allele_depths.iloc[:, 0].astype(float)
	percent_identity = (ref_depths / total_depth_per_site) * 100
	# form a new VCF df with just chrom and pos
	df = vcf_df[['CHROM','POS']]
	# add each data Series into dataframe
	df = pd.concat([df, total_depth_per_site.rename('depth'),
		percent_identity.rename('percent_identity')], axis=1)
	df = df[df.depth != 0]
	# calculate statistics
	alignment = {'length': len(df)}
	depth = {'mean': float(df.depth.mean()), 'std': float(df.depth.std())}
	s = calc_stats_from_frequency_distribution_table(
		df['depth'], df['percent_identity'])
	ani = {'mean': s['mean'], 'std': s['std'],
		'Q1': s['25%'], 'Q2': s['50%'], 'Q3': s['75%']}
	df['percent_identity'] = df['percent_identity'].round(0).astype(int)
	# print(df.dtypes)
	return df, alignment, depth, ani

def add_unmapped_sites_to_vcf(vcf_df, ref_seq_lengths, xlsx=None):
	'''
	identify positions in the input reference file that had no reads align
	and add those sites into the pandas dataframe with 0 data
	'''
	# dfs, chrom_names = [], []
	# for i, ref_seq_length in enumerate(ref_seq_lengths, start=1):
	# 	chrom_id = 'ctg_' + str(i) #cleanup_and_evaluate_mfasta() made names
	# 	mapped_positions_in_chrom = vcf_df.loc[vcf_df['CHROM'] == chrom_id]
	# 	if len(mapped_positions_in_chrom) == 0:
	# 		sys.stderr.write('INFO: no reads mapped to {}\n'.format(chrom_id))
	# 		continue
	# 	dfs.append(
	# 		mapped_positions_in_chrom.\
	# 		reindex(np.arange(1, ref_seq_length + 1)).\
	# 		fillna(0)).sort_values('POS')
	# 	chrom_names.append(chrom_id)
	# df = pd.concat(dfs, keys=chrom_names)
	# # require vcf dataframe to have a row for each reference site
	# total_ref_seq_length = sum(ref_seq_lengths)
	# positions_in_vcf = len(df)
	# if total_ref_seq_length != positions_in_vcf:
	# 	sys.stderr.write('ERROR: number of VCF rows unequal to FastA'
	# 		' reference input sequence lengths\n')
	# 	sys.exit(1)
	# return df

	# df = vcf_df.set_index(['CHROM', 'POS']).\
	# reindex((c, n) for c, v in enumerate(ref_seq_lengths, start=1) for n in np.arange(1, v+1)).reset_index()

	rotated_df = vcf_df.set_index(['CHROM', 'POS'])
	frames = []
	for chrom_name, chrom_length in enumerate(ref_seq_lengths, start=1):
		positions = np.arange(start=1, stop=chrom_length + 1)
		idx = pd.MultiIndex.from_product(([chrom_name], positions),
			names=['CHROM', 'POS'])
		# chrom_df = rotated_df.iloc[
		# 	rotated_df.index.get_level_values('CHROM') == chrom_name]
		# dupe_positions = chrom_df[chrom_df.duplicated(['POS'], keep=False)]
		# if len(dupe_positions) > 0:
		# 	sys.stderr.write('ERROR: duplicate positions in {} chrom found, '
		# 		'perhaps due to SNPs and InDels. Only one position '
		# 		'permitted, so either skip InDels or merge prior to trying '
		# 		'to fill unmapped sites.\n'.format(chrome_name))
		# 	sys.stderr.write( dupe_positions + '\n')
		# 	sys.exit(1)
		frame = rotated_df.reindex(idx, fill_value=0).\
			sort_values(by=['CHROM', 'POS']).reset_index()
		frames.append(frame)
	df = pd.concat(frames).sort_values(by=['CHROM', 'POS']).reset_index()
	positions_in_vcf = len(df)
	if positions_in_vcf != sum(ref_seq_lengths):
		sys.stderr.write('ERROR: number of VCF rows unequal to FastA'
			' reference input sequence lengths\n')
		sys.exit(1)
	if xlsx is not None:
		df.to_excel(os.path.realpath(os.path.expanduser(xlsx)))
	df = df[['depth', 'percent_identity']].to_numpy()
	return df

def filter_vcf(chrom_name=None, chrom_positions=None):
	return None

def window_average(array, window_size):
	a = np.nanmean(np.pad(array.astype(float),
		(0, 0 if array.size % window_size == 0 
			else window_size - array.size % window_size),
		mode='constant', constant_values=np.NaN).\
		reshape(-1, window_size), axis=1)
	return a

def scatter_hist(x, y, y2, chrom_lengths, ax, ax_histx, ax_histy):
	# no labels
	# x = 'Depth [bp of Reads / bp of Reference]'
	# y = 'Aligned Nucleotide Identity (ANI)'
	ax_histx.tick_params(axis='x', which='both', labelbottom=False)
	ax_histy.tick_params(axis='y', left=True, right=True, 
		which='both', labelright=True, labelrotation=90)

	# the scatter plot:
	ax.scatter(x, y2, marker='.', c='k', alpha=0.25)
	if chrom_lengths is not None:
		for length in chrom_lengths:
			ax.vlines(length, ymin=0, ymax=max(y), alpha=0.4)
	# ax.vlines([1200, 5000, 7500, 12500, 15000, 16000, 17000],
	# 	ymin=0, ymax=max(y), alpha=0.4)

	# # now determine nice limits by hand:
	# binwidth = 0.1
	# xy_max = max(np.max(np.abs(x)), np.max(np.abs(y2)))
	# lim = (int(xy_max/binwidth) + 1) * binwidth

	# bins = np.arange(1, lim + binwidth, binwidth)
	# print('bins'+str(bins))
	# bins=2000
	# ax_histx.hist(y, bins=bins)

	# ax_histx.bar(x, y, color='blue', align='center') # A bar chart
	# ax_histx.xlabel('Bins')
	# ax_histx.ylabel('Frequency')

	window_size = int(len(x) / 300)
	x = window_average(x, window_size)
	y = window_average(y, window_size)
	y2 = window_average(y2, window_size)

	ax_histx.bar(x, y, edgecolor='k')
	# ax_histx.bar(x, y, width=.01, edgecolor='k')

	# H, xedges, yedges = np.histogram2d(x, y, bins=100)
	# H = H.T  # Let each row list bins with common y range.
	# X, Y = np.meshgrid(xedges, yedges)
	# ax_histx.pcolormesh(X, Y, H)

	# window_size=200
	# x = window_average(x, window_size)
	# print(len(x))

	# ax_histy.bar(x, y2, edgecolor='k', orientation='horizontal')
	ax_histy.hist(y2, bins=x, orientation='horizontal')

	# COLORS https://matplotlib.org/api/_as_gen/matplotlib.pyplot.colors.html

def plot_scatterplot_and_histogram(x, y, y2, chrom_lengths, ani, 
	alignment, depth, num_bins, title, subtitle, outfile):
	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.55
	spacing = 0.015

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom + height + spacing, width, 0.3]
	rect_histy = [left + width + spacing, bottom, 0.2, height]
	rect_text = [left + width + spacing, bottom + height + spacing, 0.2, 0.3]

	# start with a square Figure
	fig = plt.figure(figsize=(10.5, 8))
	ax = fig.add_axes(rect_scatter)
	plt.xlabel('Reference Genome [bp]', axes=ax)
	plt.ylabel('Nucleotide Identity [%]', axes=ax)
	ax_histx = fig.add_axes(rect_histx, sharex=ax)
	plt.xlabel('', axes=ax)
	plt.ylabel('Reads Mapped (Depth of Coverage)', axes=ax)
	ax_histy = fig.add_axes(rect_histy)
	plt.xlabel('BB', axes=ax)
	plt.ylabel('Identity', axes=ax)
	text_area = fig.add_axes(rect_text, frameon=False)
	text_area.axis('off')

	# plt.ylabel('Depth', axes=ax_histx)
	# plt.xlabel('Quantity', axes=ax_histy)
	# plt.ylabel('Nucleotide Identity [%]', axes=ax_histy)

	textblock = '\n'.join((
		'ANI = {:.2f}%'.format(ani['mean']),
		r'$\sigma$ = {:.2f}%'.format(ani['std']),
		'Q1 = {:.2f}%'.format(ani['Q1']),
		'median = {:.2f}%'.format(ani['Q2']),
		'Q3 = {:.2f}%'.format(ani['Q3']),
		'Alignment = {:,} bp ({:.2f}%)'.format(alignment['length'], alignment['fraction']),
		r'Mean Depth = {:.2f}x $\pm$ {:.2f}x'.format(depth['mean'], depth['std'])
		))
	legend = dict(boxstyle='round', facecolor='grey', alpha=0.5)
	# fig.text(0.996, 0.98, textblock, transform=text_area.transAxes, fontsize=8,
	# 	fontdict={'family':'monospace'}, verticalalignment='top', bbox=legend,
	# 	multialignment='left')
	text_area.annotate(textblock, xy=(1, 1), xycoords='axes fraction',
		fontsize=8, 
		backgroundcolor='#DCDCDC', color='k', ha='right', va='top')

	# use the previously defined function
	scatter_hist(x, y, y2, chrom_lengths, ax, ax_histx, ax_histy)

	# if title is not None:
	# 	fig.suptitle(title, fontsize=10)
	# if subtitle is not None:
	# 	ax.set_title(subtitle, fontsize=8)
	if outfile.endswith('.png'):
		img_fmt = 'png'
	elif outfile.endswith('.svg'):
		img_fmt = 'svg'
	elif outfile.endswith('.eps'):
		img_fmt = 'eps'
	else:
		img_fmt = 'pdf'
	plt.savefig(outfile, pad_inches=0.5, dpi=300, format=img_fmt)

def main():
	opt = parseArgs()
	require_dependency('bowtie2-build')
	require_dependency('bowtie2')

	# I/O Handling
	require_nargs_range(opt.infile, 1, 2)	
	fastq_input = []
	for file in opt.infile:
		fastq_input.append(os.path.realpath(os.path.expanduser(file)))
	if len(fastq_input) == 2:
		mode = 'paired'
	elif len(fastq_input) == 1:
		mode = 'single'
	else:
		sys.stderr.write('ERROR: only 1 or 2 input files supported\n')
		sys.exit(1)
	if opt.outfile is not None:
		outfile = os.path.realpath(os.path.expanduser(opt.outfile))
		outdir = os.path.dirname(outfile)
		if not os.path.exists(outdir):
			os.mkdir(outdir)
	else:
		outfile = 'map-to-ref.pdf'
	input_ref = os.path.realpath(os.path.expanduser(opt.reference))
	tmp = mkdtemp()
	temp_ref = os.path.join(tmp, 'ref')
	# shutil.copy(input_ref, temp_ref)
	ref_seq_lengths = cleanup_and_evaluate_mfasta(input_ref, temp_ref)
	sys.stderr.write('INFO: reference has {} contigs and {} total length\n'.\
		format(len(ref_seq_lengths), sum(ref_seq_lengths)))
	if opt.cpus < 1:
		cpus = str(cpu_count())
	else:
		cpus = str(opt.cpus)

	# Map input reads to reference
	map_file = os.path.join(tmp, 'sample.sam')
	cmd = ['bowtie2-build', '--quiet', '--threads', cpus, temp_ref, temp_ref]
	system_call(cmd)
	for x in ['1', '2', '3', '4', 'rev.1', 'rev.2']:
		verify_file_exists_and_nonempty(temp_ref + '.' + x + '.bt2')
		# NOTE: if ref >4 Gbp the large 64 bits .bt2l extensions will fail
	cmd = ['bowtie2'] + opt.map_options.strip('\'').strip('\"').split() + \
		['-x', temp_ref, '-S', map_file]
	if mode == 'paired':
		cmd.extend(['-1', fastq_input[0], '-2', fastq_input[1]])
	elif mode == 'single':
		cmd.extend(['-U', fastq_input[0]])
	system_call(cmd)
	verify_file_exists_and_nonempty(map_file)

	# Sort the mapped file by chromosome position
	sort_file = os.path.join(tmp, 'sample.sorted.sam')
	cmd = ['samtools', 'sort', '--threads', cpus, '--output-fmt', 'SAM',
		'--reference', temp_ref, '-o', sort_file, map_file]
	system_call(cmd)
	verify_file_exists_and_nonempty(sort_file)
	# c = ['samtools', 'index', sort_file] #NOTE: cant index SAM, only BAM

	# Convert SAM to BED coordinates 0 to len(ref)
		# HOW TO get % nucl ident per site????

	# Convert mapped file to VCF
	# cmd1 = ['bcftools', 'view', '--threads', cpus, '-O', 'v', '-o', out_vcf, sort_file]
	# mpileup_file = os.path.join(tmp, 'sample.mpileup')
	# cmd2 = ['samtools', 'mpileup', '-aa', '--output-BP', '--output-MQ,' '--fasta-ref', temp_ref, '-o', mpileup_file, sort_file]
	# c1 = ['samtools', 'faidx', temp_ref]
	# system_call(c1)
	# verify_file_exists_and_nonempty(temp_ref + '.fai')
	vcf_file = os.path.join(tmp, 'sample.vcf')
	cmd = ['bcftools', 'mpileup', '--annotate', 'FORMAT/AD,INFO/AD',
		'--skip-indels', '--max-depth', '500', '--threads', cpus,
		'--fasta-ref', temp_ref, '-o', vcf_file,sort_file]
	system_call(cmd)
	verify_file_exists_and_nonempty(vcf_file)
#cd ~/test_plot.map-to-ref/2020_ND/3001927717_BacND2019/ANI-from-BAM
#bcftools mpileup --max-depth 500 -a "FORMAT/AD,INFO/AD" --fasta-ref ../../AmesAncestor/pXO2ISelemStretchWithcapABCDE.fasta -o vcf_file output.sam.sorted

	# Convert VCF file to Pandas DataFrame, then to Numpy N-dimensional array
	df_mapped_with_stats, alignment, depth, ani = convert_vcf_to_stats(
		read_vcf(vcf_file))
	shutil.rmtree(tmp)
	# convert df to class 'numpy.ndarray'
	arr_mapped_and_unmapped = add_unmapped_sites_to_vcf(
		df_mapped_with_stats, ref_seq_lengths, xlsx=opt.xlsx)
	# df_filt	= df[(df['depth'] >= opt.min_depth) &
	# 	(df['percent_identity'] > opt.min_identity)]
	# df_filt['Depth*Ident'] = df_filt['depth'] * df_filt['percent_identity']
	# cnt_identity = df_filt['Depth*Ident'].sum()
	# cnt_total_depth = df_filt['depth'].sum()
	# ani = 100. * cnt_identity / cnt_total_depth
	# percent_aligned = 100. * len(df_filt) / ref_length
	alignment['fraction'] = 100. * alignment['length'] / sum(ref_seq_lengths)

	# # df_filt = df_filt.round({'percent_identity': 0})  #problematic 0.5 val round
	# cnt_identity = np.around(df_filt['percent_identity'], 0).sum()
	# cnt_total_aln = len(df_filt)
	# ani = 100. * cnt_identity / cnt_total_aln
	# # df_filt = np.around(df)

	# Bin the coordinates as a list
	#### collapse df into (index_position, depth, identity)
	depths = arr_mapped_and_unmapped[:, 0] #dtype is int64
	per_id = arr_mapped_and_unmapped[:, 1] #dtype is int64

	# Plot Coords+Depth then Coords+%Ident
	coordinates = np.arange(1, sum(ref_seq_lengths) + 1)
	if len(ref_seq_lengths) > 1:
		chrom_lengths = np.add.accumulate(np.array(ref_seq_lengths))
	else:
		chrom_lengths = None
	plot_scatterplot_and_histogram(coordinates, depths, per_id, 
		chrom_lengths, ani, alignment, depth, 200, 'title', 'subtitle', outfile)

if __name__ == '__main__':
	main()
