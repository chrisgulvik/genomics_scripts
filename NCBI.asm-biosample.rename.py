#!/usr/bin/env python3


import csv
import json
import os
import re
import shutil
import sys
from argparse import ArgumentParser
from glob import glob

def parseArgs():
	parser = ArgumentParser(description='Renames assembly files',
		add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--indir', required=True, metavar='FILE',
		help='input path containing files to rename')
	req.add_argument('-j', '--json', required=True, metavar='FILE',
		help='input JSON database file of BioSample metadata')
	req.add_argument('-k', '--keys', required=True, metavar='FILE',
		help='TSV file; first column are BioSample Accessions'
		' and second column are corresponding strings (e.g.,'
		' Assembly Accessions) in input filenames to match')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outdir', metavar='FILE', default=None,
		help='path for renamed files [cwd]')
	rno = parser.add_argument_group('Renaming Options')
	rno.add_argument('--empty', metavar='STR', type=str, default='missing',
		help='character or string for empty values [missing]')
	rno.add_argument('--metadata', choices=['geo_loc_name',	'isolate', 'SRA',
		'strain'], nargs='+', default=['geo_loc_name'],
		help='biosample field(s) from JSON input to add to corresponding'
		' filenames [geo_loc_name]')
	rno.add_argument('--save-method', choices=['copy', 'move', 'symlink'],
		default='symlink',
		help='method for saving renamed files [symlink]')
	rno.add_argument('--substring-rename-method', default='append-suffix',
		choices=['append-prefix', 'append-suffix', 'replace'],
		help='filename renaming strategy [append-suffix]')
	fso = parser.add_argument_group('File Matching Options')
	fso.add_argument('--extension', metavar='STR', type=str,
		default=['.gbff.gz'], nargs='+',
		help='extension(s) of files to rename [.gbff.gz]')
	fso.add_argument('--no-recursive', action='store_true', default=False,
		help='turn off recursive search for files within indir [off]')
	qry = parser.add_argument_group('Word Query (2nd col --keys <TSV>) Options')
	qry.add_argument('--query-absent', default='fail',
		choices=['fail', 'skip'],
		help='when word not found in any files [fail]')
	# qry.add_argument('--query-type', default='filenames',
	# 	choices=['filenames', 'content'],
	# 	help='where to search for query words [filenames]')
	qry.add_argument('--substring-match-location', default='anywhere',
		choices=['anywhere', 'full', 'prefix', 'suffix'],
		help='where TSV values exist in each filename [anywhere]')
	return parser.parse_args()

def get_metadata(queries, data, field, search_type):
	'''finds each query string from a list of queries in a data dictionary
	with specified search field and search type'''
	all_found_records = {}
	if search_type == 'anywhere':
		for qry in queries:
			found = {k:v for k, v in data.items() if qry in v[field]}
			if len(found) == 0:
				sys.stderr.write('ERROR: {} {} absent\n'.format(qry, field))
				sys.exit(1)
			all_found_records.update(found)
	elif search_type == 'full':
		for qry in queries:
			found = {k: v for k, v in data.items() if v[field] == qry}
			if len(found) == 0:
				sys.stderr.write('ERROR: {} {} absent\n'.format(qry, field))
				sys.exit(1)
			all_found_records.update(found)
	elif search_type == 'prefix':
		for qry in queries:
			found = {k:v for k, v in data.items() if v[field].startswith(qry)}
			if len(found) == 0:
				sys.stderr.write('ERROR: {} {} absent\n'.format(qry, field))
				sys.exit(1)
			all_found_records.update(found)
	elif search_type == 'suffix':
		for qry in queries:
			found = {k:v for k, v in data.items() if v[field].endswith(qry)}
			if len(found) == 0:
				sys.stderr.write('ERROR: {} {} absent\n'.format(qry, field))
				sys.exit(1)
			all_found_records.update(found)
	elif search_type == 'key':
		for qry in queries:
			found = data.get(qry, None)
			if found is None:
				sys.stderr.write('ERROR: {} BioSample absent\n'.format(qry))
				sys.exit(1)
			all_found_records.update({qry: found})
	sys.stderr.write('INFO: {} entries remain after {} filter\n'.format(
		len(all_found_records), field))
	return all_found_records

def find_files(query, extensions, files, search_type, query_absent):
	'''finds files within a list of files with the filename containing a query
	 string with specified file extensions and search type'''
	found_files = []
	for file in files:
		f = os.path.basename(file)
		if search_type == 'anywhere':
			if query in f:
				found_files.append(file)
		elif search_type == 'full':
			for ext in extensions:
				b = f.rstrip(ext)
				if query == b:
					found_files.append(file)
		elif search_type == 'prefix':
			if f.startswith(query):
				found_files.append(file)
		elif search_type == 'suffix':
			for ext in extensions:
				if f.endswith(query + ext):
					found_files.append(file)
	if len(found_files) == 0 and query_absent == 'fail':
		sys.stderr.write('ERROR: no filenames contain {} {}\n'.format(query,
			search_type))
		sys.exit(1)
	return found_files

def main():
	opt = parseArgs()
	recursive = not opt.no_recursive
	if sys.version_info < (3, 5):
		sys.stderr.write('WARNING: Python 3.5+ is required to do the'
			'recursive file searching with the stdlib simply in glob\n')
		recursive = False
	indir = os.path.realpath(os.path.expanduser(opt.indir))
	extensions = opt.extension
	match_location = opt.substring_match_location
	if opt.outdir is None:
		outdir = os.getcwd()
	else:
		outdir = os.path.realpath(os.path.expanduser(opt.outdir))
		if not os.path.exists(outdir):
			os.makedirs(outdir)

	# Identify files with matching extensions
	files = []
	for ext in extensions:
		if not opt.no_recursive:
			found = glob(os.path.join(indir, '**', '*' + ext), recursive=True)
		else:
			found = glob(os.path.join(indir, '*' + ext))
		if len(found) > 0:
			files.extend(found)
	if len(files) == 0:
		sys.stderr.write('ERROR: no files in {} with {} extension\n'.format(
			indir, ','.join(extensions)))
		sys.exit(1)

	# Load TSV input keys
	rename_keys = {}
	with open(os.path.realpath(os.path.expanduser(opt.keys))) as ifh:
		reader = csv.DictReader(ifh, fieldnames=['BioSample', 'Query_Word'],
			delimiter='\t')
		biosample_regex = re.compile('^SAM(D|N|E([AG]?))[0-9]+$')
		for row in reader:
			if not bool(re.match(biosample_regex, row['BioSample'])):
				sys.stderr.write('ERROR: first column must be a BioSample'
					' accession and {} doesn\'t appear to be one. EBI'
					' explains, "BioSample accessions always begin with SAM.'
					' The next letter is either E or N or D depending if the'
					' sample information was originally submitted to EBI or'
					' NCBI or DDBJ respectively. After that, there may be an'
					' A or a G to denote an Assay sample or a Group of'
					' samples. Finally there is a numeric component that may'
					' or may not be zero-padded."\n'.format(row['BioSample']))
				sys.exit(1)
			rename_keys[str(row['BioSample'])] = row['Query_Word']
	sys.stderr.write('INFO: input renaming keys has {} entries\n'.format(
		len(rename_keys)))

	# Load JSON input with biosample accessions as keys
	with open(os.path.realpath(os.path.expanduser(opt.json))) as ifh:
		json_d = json.load(ifh)
	cnt_biosamples = len(json_d)
	sys.stderr.write('INFO: input has {} biosample entries\n'.format(
		cnt_biosamples))

	# Compare rename biosample key input with JSON input
	absent_biosamples = [s for s in rename_keys.keys() if not s in json_d]
	cnt_absent = len(absent_biosamples)
	if cnt_absent == len(rename_keys):
		sys.stderr.write('ERROR: all keys provided are absent in JSON'
			' database\n')
		sys.exit(1)
	elif cnt_absent > 0:
		sys.stderr.write('WARNING: {} keys provided are absent in JSON'
			' database and will not be renamed:\n  {}\n'.format(cnt_absent,
			' '.join(sorted(absent_biosamples))))
		for key in absent_biosamples:
			del rename_keys[key]

	# Rename files
	cnt_renamed = 0
	for biosample, qry_word in rename_keys.items():
		bs_data = json_d[biosample]
		new_name = ''
		for wanted_field in opt.metadata:
			if wanted_field in bs_data:
				# new_name += '_' + bs_data[wanted_field].split(':')[0]
				new_name += '_' + bs_data[wanted_field]
			else:
				new_name += '_' + opt.empty
		found_files = find_files(qry_word, extensions, files, match_location,
			opt.query_absent)
		for f in found_files:
			b = os.path.basename(f)
			if opt.substring_rename_method == 'append-prefix':
				filename = os.path.join(outdir, new_name + b)
			elif opt.substring_rename_method == 'append-suffix':
				idx = [i for i, s in enumerate(extensions) if b.endswith(s)]
				for idx, ext in enumerate(extensions):
					if b.endswith(ext):
						filename = (b.rstrip(extensions[idx]) + new_name
						+ extensions[idx])
						break
			elif opt.substring_rename_method == 'replace':
				filename = b.replace(qry_word, new_name)
			dest = os.path.join(outdir, filename)
			if opt.save_method == 'copy':
				shutil.copyfile(f, dest)
			elif opt.save_method == 'move':
				shutil.move(f, dest)
			elif opt.save_method == 'symlink':
				try:
					os.symlink(os.path.relpath(f, os.path.dirname(dest)), dest)
				except FileExistsError:
					pass
			cnt_renamed += 1
	sys.stderr.write('INFO: {} files renamed\n'.format(cnt_renamed))

if __name__ == '__main__':
	main()
