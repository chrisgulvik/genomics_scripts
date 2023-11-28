#!/usr/bin/env python


import os
import string
import sys
from argparse import ArgumentParser
from decimal import Decimal, ROUND_HALF_UP
from subprocess import check_output, Popen
from time import strftime
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
	parser = ArgumentParser(description='Extracts proteins from a FastA or GenBank proteome. Providing a GenBank file extracts both protein and gene sequence hits, whereas a FastA input only allows for protein extractions.',
	epilog='Output directory structure: <outpath>/{indiv,log,sum}/<base>.<filenames>.<ext> Individually extracted loci are named <base>.<locus_tag>.<acc_in_db> with fsa extension for proteins and fas extension for genes. These are subsequently merged in the sum dir as <base>_found with faa extension for proteins and ffn extension for genes and summarized with a <base>.loci_abundance.tab file useful for plotting in R.')
	parser.add_argument('-i', '--infile', required=True, help='input proteome FastA file <.faa> or GenBank file <.gb||.gbk> with translations')
	parser.add_argument('-m', '--hmm', required=True, help='input HMM database file')
	parser.add_argument('-b', '--base', required=False, default=None, help='prefix string to append to output files [`basename infile .<ext>`]')
	parser.add_argument('-o', '--outpath', required=False, default=None, help='output directory  [./base--<date>_<time>]')
	parser.add_argument('-c', '--cpus', required=False, type=int, default='1', help='number of CPUs [1]')
	args = parser.parse_args()
	return args

def genbank2faa(in_queryfile, outpath):
	if in_queryfile.endswith('.gb'):
		proteome = os.path.basename(in_queryfile).replace('.gb', '.faa')
	elif in_queryfile.endswith('.gbk'):
		proteome = os.path.basename(in_queryfile).replace('.gbk', '.faa')
	with open(os.path.join(outpath, proteome), 'w') as p:
		for rec in SeqIO.parse(in_queryfile, 'genbank'):
			for feat in rec.features:
				if feat.type == 'CDS':
					p.write('>{} {}\n{}\n'.format(
						feat.qualifiers['locus_tag'][0],
						feat.qualifiers['product'][0],
						feat.qualifiers['translation'][0]))
	return os.path.join(outpath, proteome)

def index_genbank_features(gb_record, feature_type, qualifier) :
	gb_idx = dict()
	for (i, feat) in enumerate(gb_record.features):
		if feat.type == feature_type:
			if qualifier in feat.qualifiers:
				for s in feat.qualifiers[qualifier]:  #list of loci
					if s in gb_idx:  #make sure the locus_tag is unique in the record
						sys.exit('ERROR: malformed GenBank file; found duplicate of {}'.format(s))
					else:
						gb_idx[s] = i
	return gb_idx

def locus_tag2gene_seq(gbk, locus_tag, outfile):
	r = 0
	for gb_record in SeqIO.parse(open(gbk, 'r'), 'genbank'):
		locus_tag_index = index_genbank_features(gb_record, 'CDS', 'locus_tag')
		if locus_tag in locus_tag_index:
			gene_seq = str(gb_record.features[locus_tag_index[locus_tag]].extract(gb_record.seq))
			product = str(gb_record.features[locus_tag_index[locus_tag]].qualifiers['product'][0])
			record = SeqRecord(Seq(gene_seq, generic_dna), id=locus_tag, description=product, name='')
			r += 1
	if r == 1:
		SeqIO.write(record, outfile, 'fasta')  #gene seq
	elif r > 1:
		sys.exit('ERROR: malformed GenkBank file; found >1 occurance of {}'.format(locus_tag))
	else:
		sys.exit('ERROR: malformed GenkBank file; absent locus? {}'.format(locus_tag))

def main():
	opts = parseArgs()
	cpus = opts.cpus
	in_queryfile = str(os.path.abspath(opts.infile))  # Example: proteome.faa
	hmm = os.path.abspath(opts.hmm)                   # Example: prot-markers.hmm
	if opts.base:
		base = opts.base
	else:
		base = os.path.basename(os.path.abspath(opts.infile)).rsplit('.')[0]
	if opts.outpath is not None:
		outpath = os.path.abspath(opts.outpath)
		if not os.path.exists(outpath):
			os.mkdir(outpath)
	else:
		autocreated_dirname = base + '--' + strftime('%d%b%Y_%-I:%M%p').upper()
		outpath = os.path.join(os.getcwd(), autocreated_dirname)
		os.mkdir(outpath)

	# Setup dir structure
	for d in ['log', 'indiv', 'sum']:
		try:
			os.mkdir(os.path.join(outpath, d))
		except OSError as err:  #OSError: [Errno 17] File exists
			print 'WARNING: non-empty output directory {}/{} might result in overwritten files'.format(outpath, d)

	# Check input sequence file format
	if in_queryfile.endswith(('.gb', '.gbk')):
		gene = True
		proteome = genbank2faa(in_queryfile, os.path.join(outpath, 'sum'))
	elif in_queryfile.endswith(('.faa', '.fa', '.fasta', '.FastA', '.fas', '.fsa')):
		gene = False
		proteome = in_queryfile
	else:
		sys.exit('ERROR: unsupported {} file extension'.format(in_queryfile.split('.')[-1:]))

	# Verify HMMsearch is v3.0+
	hmm_ver = check_output('hmmsearch -h | head -n 4', shell=True).rstrip('\n')
	if '# HMMER 3.' not in hmm_ver:
		sys.exit('ERROR: HMMER v3.0+ required')

	# Run HMMsearch with proteome input against a HMM database
	search_results = os.path.join(outpath, 'log', base + '.targethits.txt')
	cmd_search = 'hmmsearch --tblout {} --cut_tc --notextw --cpu {} {} {}'.format(search_results, cpus, hmm, proteome)
	with open(os.path.join(outpath, 'log', base + '.targethits.log'), 'w') as f:
		return_code = Popen(cmd_search.split(), stdout=f, shell=False)
		if return_code.wait() != 0:
			sys.exit('ERROR: failed system call\n{}\n'.format(cmd_search))

	# Count number of HMM profiles queried
	cmd_getloci = "hmmstat {} | grep -Ev '(^#|^$)' | awk '{{print $2}}' | tr '\\n' '\\t'".format(hmm)
	query_names = check_output(cmd_getloci, shell=True)
	query_names_list = query_names.rstrip('\t').split('\t')

	# Parse HMM search results
	input_records = SeqIO.index(proteome, 'fasta')
	query_names_found = {}  #key=query_name;val=(target_name,accession,description_of_target)
	with open(search_results, 'r') as results:
		for line in results:
			if not line.startswith("#"):
				data = line.split()
				if data[2] not in query_names_found:
					query_names_found[data[2]] = (data[0],data[3],data[18])
				else:
					old_val = query_names_found.pop(data[2])
					query_names_found[data[2]] = [old_val, (data[0],data[3],data[18])]
				SeqIO.write(input_records[data[0]], '{}.fsa'.format(
					os.path.join(outpath, 'indiv', base + '.' + data[0].split()[0] + '.' + data[3].split()[0])),
					'fasta') #protein seq
				if gene:
					locus_tag2gene_seq(in_queryfile, input_records[data[0]].id,
						os.path.join(outpath, 'indiv', base + '.' + data[0].split()[0] + '.' + data[3].split()[0] + '.fas'))
		input_records.close()

	# Create merged/summary extracted sequence file(s)
	if gene:
		os.remove(os.path.join(outpath, 'sum', base + '.faa'))  #rm tmp proteome from gbk input
		os.system('ls {} | xargs cat > {}'.format(os.path.join(outpath, 'indiv', '*.fas'), os.path.join(outpath, 'sum', base + '.found.ffn')))  #genes
	os.system('ls {} | xargs cat > {}'.format(os.path.join(outpath, 'indiv', '*.fsa'), os.path.join(outpath, 'sum', base + '.found.faa')))  #proteins

	# Generate input for R
	counts = []
	for qn in query_names_list:
		if qn in query_names_found:
			if type(query_names_found.get(qn)) is tuple:
				counts.append('1')
				hit = query_names_found.get(qn)
			elif type(query_names_found.get(qn)) is list:
				counts.append(str(len(query_names_found.get(qn))))
				hits = query_names_found.get(qn)
		else:
			counts.append('0')
	with open(os.path.join(outpath, 'sum', '{}.loci_abundance.tab'.format(base)), 'w') as outabund:
		outabund.write('#Sample\t{}\n'.format('\t'.join(query_names_list)))
		outabund.write('{}\t{}\n'.format(base, '\t'.join(counts)))

	# Count number of HMM profiles queried
	cmd_stat = "hmmstat {} | grep -Ev '(^#|^$)' | wc -l".format(hmm)
	num_profiles = check_output(cmd_stat, shell=True).rstrip('\n')

	# Report absent loci, single hits, and multiples
	num_absent  = counts.count('0')
	num_single  = counts.count('1')
	num_gt_once = len(counts) - num_absent - num_single
	with open(os.path.join(outpath, 'sum', '{}.stats.tab'.format(base)), 'w') as outsum:
		outsum.write('{:>7}% or {}/{} loci not found\n'.format(
			str(Decimal(num_absent/float(num_profiles)*100).quantize(Decimal('.01'), rounding=ROUND_HALF_UP)),
			str(num_absent), num_profiles))
		outsum.write('{:>7}% or {}/{} loci found once\n'.format(
			str(Decimal(num_single/float(num_profiles)*100).quantize(Decimal('.01'), rounding=ROUND_HALF_UP)),
			str(num_single), num_profiles))
		outsum.write('{:>7}% or {}/{} loci found more than once'.format(
			str(Decimal(num_gt_once/float(num_profiles)*100).quantize(Decimal('.01'), rounding=ROUND_HALF_UP)),
			str(num_gt_once), num_profiles))

if __name__ == '__main__':
	main()
