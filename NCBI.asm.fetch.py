#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
import sqlite3

def parseArgs():
	parser = ArgumentParser(description='Fetches assembly files from NCBI\'s '
		'web server with rsync. Unless a query is provided to find a subset, '
		'all database entries will be fetched.', add_help=False,
		epilog='NOTE: files are named exactly as listed on NCBI')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--input', required=True, metavar='FILE',
		help='input SQLite3 database file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--file-conflicts', choices=['force', 'keep',
		'update'], default='update',
		help='how file conflicts are handled; force always overwrites local, '
		'keep always retains local, and update only overwrites local when '
		'NCBI has the same file with a newer timestamp [update]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-l', '--list', metavar='FILE', default=False,
		help='output list of files (line-by-line) to fetch but skip '
		'downloading them')
	opt.add_argument('-m', '--format', choices=['fna', 'gbff', 'gff'],
		default='gbff', help='file format to fetch [gbff]')
	opt.add_argument('-o', '--outpath', metavar='PATH',
		help='output path where retrieved files are stored [cwd]')
	opt.add_argument('-q', '--quiet', action='store_true', default=False,
		help='hide messages when entries are updated or added [off]')
	opt.add_argument('-t', '--table', default='asm_summary', type=str,
		help='SQL table name to fetch entries [asm_summary]')
	opt.add_argument('-s', '--query-search', default=None,
		help='term to search within the database\'s query-feature [None]')
	opt.add_argument('-u', '--query-feature', choices=['assembly_accession',
		'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid',
		'species_taxid', 'organism_name', 'infraspecific_name', 'isolate',
		'version_status', 'assembly_level', 'release_type', 'genome_rep',
		'seq_rel_date', 'asm_name', 'submitter', 'gbrs_paired_asm',
		'paired_asm_comp', 'ftp_path', 'excluded_from_refseq',
		'relation_to_type_material'], metavar='<db col title>', 
		default='assembly_accession',
		help='feature (column) in database to query (e.g., bioproject, '
		'biosample, taxid, organism_name, asm_name) [assembly_accession]')
	opt.add_argument('-y', '--query-type', choices=['exact', 'substring'],
		default='exact', help='fetch all assemblies where query term is '
		'found within the entry or matches the full entry [exact]')
	return parser.parse_args()

def sql_open(sql_out):
	try:
		con = sqlite3.connect(sql_out)
		cur = con.cursor()
		return con, cur
	except Error as e:
		sys.stderr.write('ERROR: {}\n'.format(e))
		sys.exit(1)

def sql_close(con, cur):
	cur.close()
	con.close()

def make_list_of_webfiles(sql_data, file_format):
	if len(sql_data) == 0:
		sys.stderr.write('ERROR: no entries match the query\n')
		sys.exit(1)
	get = []
	for row in sql_data:
		if not row[0].startswith('ftp:'):
			sys.stderr.write('ERROR: expect FTP path but instead found '
				'{}\n'.format(row[0]))
			sys.exit(1)
		b = os.path.basename(row[0])
		get.append('{}/{}_genomic.{}.gz'.format(row[0], b, file_format))
	return get

def main():
	opt = parseArgs()

	# I/O handling
	ifh = os.path.abspath(os.path.expanduser(opt.input))
	if opt.outpath is not None:
		out = os.path.abspath(os.path.expanduser(opt.outpath))
	else:
		out = os.getcwd()
	if os.path.exists(out) and opt.file_conflicts == 'update' \
	and not opt.quiet:
		sys.stderr.write('WARNING: conflicting files will be overwritten if '
			'a newer assembly file is found...\n\n')
	elif not os.path.exists(out):
		os.mkdir(out)

	# Gather paths of assembly files to fetch
	tbl = opt.table
	con, cur = sql_open(ifh)
	reload(sys)
	sys.setdefaultencoding('UTF-8')
	if opt.query_search is None:
		cur.execute('SELECT ftp_path FROM {}'.format(tbl))
		dat = cur.fetchall()
		get = make_list_of_webfiles(dat, opt.format)
	else:
		if opt.query_type == 'exact':
			cur.execute('SELECT ftp_path FROM {} WHERE ({} = \'{}\' '
				'AND version_status = "latest")'.format(
				tbl, opt.query_feature, opt.query_search))
		elif opt.query_type == 'substring':
			cur.execute('SELECT ftp_path FROM {} WHERE ({} LIKE \'%{}%\' '
				'AND version_status = "latest")'.format(
				tbl, opt.query_feature, opt.query_search))
		dat = cur.fetchall()
		get = make_list_of_webfiles(dat, opt.format)
	sql_close(con, cur)

	# Fetch assembly files (or just output webfile list)
	if not opt.quiet:
		sys.stderr.write('INFO: fetching {} assemblies...\n\n'.format(
			len(get)))
	if not opt.list:
		sys_cmd = 'rsync --copy-links --times --human-readable '
		if opt.quiet:
			sys_cmd += '--quiet '
		else:
			sys_cmd += '--verbose '
		if opt.file_conflicts == 'update':
			sys_cmd += '--update '
		elif opt.file_conflicts == 'force':
			sys_cmd += '--ignore-times '
		elif opt.file_conflicts == 'keep':
			sys_cmd += '--ignore-existing '
		for url in get:
			url = url.replace('ftp:', 'rsync:')
			os.system(sys_cmd + '{} {} 1>&2'.format(url, out))
			f = os.path.join(out, os.path.basename(url))
			if (not os.path.exists(f) or os.path.getsize(f) == 0) \
			and not opt.quiet:
				sys.stderr.write('ERROR: unable to fetch {}\n'.format(url))
	else:
		with open(os.path.abspath(os.path.expanduser(opt.list)), 'w') as o:
			for ln in get:
				o.write('{}\n'.format(''.join(ln)))

if __name__ == '__main__':
	main()
