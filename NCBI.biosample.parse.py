#!/usr/bin/env python3


import csv
import json
import os
import sys
from argparse import ArgumentParser

def parseArgs():
	parser = ArgumentParser(description='Parses data from a JSON database of'
		' BioSample information from NCBI.biosample.summary.py output',
		add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input JSON database file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-e', '--empty', metavar='STR', type=str, default='',
		help='character or string to denote empty value [None]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-j', '--out-json', metavar='FILE', default=None,
		help='JSON output file [None]')
	opt.add_argument('-o', '--output', metavar='FILE', default=None,
		help='output tab-delimited file [stdout]')
	opt.add_argument('--max-records', metavar='INT', type=int, default=None,
		help='maximum number of records to output [None]')
	opt.add_argument('--min-records', metavar='INT', type=int, default=1,
		help='minimum number of records to output [1]')
	qry = parser.add_argument_group('Query Fields to Output')
	qry.add_argument('--find-biosample', metavar='STR', type=str, nargs='*',
		default=None, help='biosample identifier(s) [all]')
	qry.add_argument('--find-organism', metavar='STR', type=str, nargs='*',
		default=None, help='organism name prefix(es) [all]')
	qry.add_argument('--find-taxid', metavar='INT', type=int, nargs='*',
		default=None, help='taxonomic identifier(s) assigned by NCBI [all]')
	qry.add_argument('--find-taxname', metavar='STR', type=str, nargs='*',
		default=None, help='taxonomic name prefix(es) [all]')
	xtr = parser.add_argument_group('Additional Fields to Output')
	xtr.add_argument('--collected-by', action='store_true', default=False,
		help='name of persons or institute who collected the sample')
	xtr.add_argument('--collection-date', action='store_true', default=False,
		help='collection date')
	xtr.add_argument('--description', action='store_true', default=False,
		help='sample\'s description')
	xtr.add_argument('--geo', action='store_true', default=False,
		help='geographic origin')
	xtr.add_argument('--host', action='store_true', default=False,
		help='organism from which the sample was obtained')
	xtr.add_argument('--isolate', action='store_true', default=False,
		help='isolate number')
	xtr.add_argument('--sample-type', action='store_true', default=False,
		help='sample type')
	xtr.add_argument('--sra', action='store_true', default=False,
		help='Sequence Read Archive (SRA) accession')
	xtr.add_argument('--strain', action='store_true', default=False,
		help='strain number')
	return parser.parse_args()

def find_entries(queries, data, field, search_type):
	all_found_records = {}
	if search_type == 'full':
		for qry in queries:
			try:
				found = {k: v for k, v in data.items() if v[field] == qry}
			except KeyError:
				sys.stderr.write('ERROR: {} absent in db\n'.format(field))
				sys.exit(1)
			if len(found) == 0:
				sys.stderr.write('ERROR: {} {} absent\n'.format(qry, field))
				sys.exit(1)
			all_found_records.update(found)
	elif search_type == 'prefix':
		for qry in queries:
			try:
				found = {k: v for k, v in data.items() if v[field].startswith(qry)}
			except KeyError:
				sys.stderr.write('ERROR: {} absent in db\n'.format(field))
				sys.exit(1)
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

def main():
	opt = parseArgs()
	
	# Load JSON input
	with open(os.path.abspath(os.path.expanduser(opt.infile))) as ifh:
		json_d = json.load(ifh)
	cnt_biosamples = len(json_d)
	sys.stderr.write('INFO: input has {} biosample entries\n'.format(
		cnt_biosamples))
	if cnt_biosamples < 1:
		sys.stderr.write('ERROR: at least one biosample entry required\n'.\
			format(cnt_biosamples))
		sys.exit(1)

	# Record filtering
	find_l = [(opt.find_biosample, 'BioSample', 'key'),
		(opt.find_organism, 'Organism', 'prefix'),
		(opt.find_taxid, 'TaxID', 'full'),
		(opt.find_taxname, 'Taxonomy_Name', 'prefix')]
	for queries, field, search_type in find_l:
		if queries is not None:
			json_d = find_entries(queries, json_d, field, search_type)

	# Save additional fields specified
	want = ['Organism', 'Taxonomy_Name']
	if opt.collected_by:
		want.append('collected_by')
	if opt.collection_date:
		want.append('collection_date')
	if opt.description:
		want.append('description')
	if opt.geo:
		want.append('geo_loc_name')
	if opt.host:
		want.append('host')
	if opt.isolate:
		want.append('isolate')
	if opt.sample_type:
		want.append('sample_type')
	if opt.sra:
		want.append('SRA')
	if opt.strain:
		want.append('strain')
	d = {}
	for biosample, val in json_d.items():
		d[biosample] = {k: v for k, v in val.items() if k in want}

	# Quit if too few or too many records found
	num_records = len(d)
	if num_records < opt.min_records:
		sys.stderr.write('ERROR: {} records matching query less than {}'
			' --min-records specified\n'.format(num_records,
			opt.min_records))
		sys.exit(1)
	if opt.max_records is not None:
		if opt.max_records < num_records:
			sys.stderr.write('ERROR: {} records matching query greater than'
				' {} --max-records specified\n'.format(num_records,
				opt.max_records))
			sys.exit(1)

	# Optionally write data as JSON
	if opt.out_json is not None:
		out_json = os.path.realpath(os.path.expanduser(opt.out_json))
		with open(out_json, 'w') as o:
			json.dump(d, o)

	# Output TSV data
	if opt.output is not None:
		ofh = open(os.path.realpath(os.path.expanduser(opt.output)), 'w')
	else:
		ofh = sys.stdout
	w = csv.DictWriter(ofh, fieldnames=['BioSample'] + want,
		delimiter='\t', extrasaction='ignore')
	w.writeheader()
	for k in sorted(d):
		data = {field: d[k].get(field, opt.empty) for field in want}
		w.writerow({'BioSample': k, **data})

if __name__ == '__main__':
	main()
