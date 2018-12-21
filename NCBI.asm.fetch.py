#!/usr/bin/env python


import multiprocessing as mp
import os
import subprocess as sp
import sqlite3
import sys
from argparse import ArgumentParser
from shutil import rmtree
from tempfile import mkdtemp

def parseArgs():
	parser = ArgumentParser(description='Fetches assembly files from NCBI\'s '
		'web server with rsync or wget. Unless a query is provided to find a '
		'subset, all database entries will be fetched.', add_help=False,
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
	opt.add_argument('-m', '--format', choices=['fna', 'gbff', 'gff'],
		default='gbff', help='assembly file format to fetch [gbff]')
	opt.add_argument('-o', '--outpath', metavar='PATH',
		help='output path where retrieved files are stored [cwd]')
	opt.add_argument('-p', '--protocol', choices=['rsync', 'wget'],
		default='rsync', help='download protocol [rsync]')
	opt.add_argument('-q', '--quiet', action='store_true', default=False,
		help='hide messages when entries are updated or added [off]')
	opt.add_argument('-t', '--table', metavar='STR', default='asm_summary',
		 type=str, help='SQL table name to fetch entries [asm_summary]')
	opt.add_argument('-s', '--query-search', metavar='STR', default=None,
		help='term to search within the database\'s query-feature [None]')
	opt.add_argument('-u', '--query-feature', choices=['assembly_accession',
		'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid',
		'species_taxid', 'organism_name', 'infraspecific_name', 'isolate',
		'version_status', 'assembly_level', 'release_type', 'genome_rep',
		'seq_rel_date', 'asm_name', 'submitter', 'gbrs_paired_asm',
		'paired_asm_comp', 'ftp_path', 'excluded_from_refseq',
		'relation_to_type_material'], metavar='STR', 
		default='assembly_accession',
		help='feature (column) in database to query (e.g., bioproject, '
		'biosample, taxid, organism_name, asm_name) [assembly_accession]')
	opt.add_argument('-y', '--query-type', choices=['exact', 'substring'],
		default='exact', help='fetch all assemblies where query term is '
		'found within the entry or matches the full entry [exact]')
	opt.add_argument('--connections', metavar='INT', default=1, type=int,
		help='number of concurrent download connections [1]')
	opt.add_argument('--connect-timeout', metavar='INT', default=3, type=int,
		help='seconds to wait for a connection to NCBI before quitting [3]')
	opt.add_argument('--info', metavar='FILE', default=None,
		help='output tab-delimited list of sample metadata fetched [None]')
	opt.add_argument('--list', metavar='FILE', default=None,
		help='output list of URLs (line-by-line); useful if you want to do '
		'your own rsync --files-from or wget --input-file [None]')
	opt.add_argument('--no-download', action='store_true', default=False,
		help='skips downloading assembly files [off]')
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
		if not row[19].startswith('ftp:'):
			sys.stderr.write('ERROR: expect FTP path but instead found '
				'{}\n'.format(row[19]))
			sys.exit(1)
		b = os.path.basename(row[19])
		get.append('{}/{}_genomic.{}.gz'.format(row[19], b, file_format))
	return get

def generate_metadata(sql_data, file_format):
	''' returns a list, where each item is a tab-delimited string. The first
	column is the full filepath to fetch the file from NCBI's FTP server,
	the remaining columns are identical to NCBI's assembly summary format
	which hold 22 columns, e.g., Assembly Accn, BioSample, Organism, etc. '''
	nfo = []
	for row in sql_data:
		b = os.path.basename(row[19])
		# need a char for empty cells to read while-loop array with IFS tab
		metadata = [str(s) if str(s).strip() else '.' for s in row]
		nfo.append('{}/{}_genomic.{}.gz\t{}'.format(
			row[19], b, file_format, '\t'.join(metadata)))
	return nfo

def write_file_list(outfile, url_list, bool_quiet):
	with open(outfile, 'w') as o:
		for ln in url_list:
			o.write('{}\n'.format(''.join(ln)))
	if not bool_quiet:
		sys.stderr.write('INFO: saved file list for {} assemblies\n'.\
			format(len(url_list)))

def fetch_asm(sys_cmd_string):
	sys_cmd_list = sys_cmd_string.split()
	with open(os.devnull, 'wb') as dump:
		process = sp.Popen(sys_cmd_list, stdout=dump, stderr=sp.PIPE)
		_, err = process.communicate()
		sys.stderr.write(err)
		if process.returncode == 0:
			return True
		else:	
			sys.stderr.write('ERROR: failed system call: {}\n'.\
				format(sys_cmd_string))
			return False

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
			'a newer assembly file is found...\n')
	elif not os.path.exists(out):
		os.mkdir(out)
	cpu_processes = opt.connections

	# Gather paths of assembly files
	tbl = opt.table
	con, cur = sql_open(ifh)
	reload(sys)
	sys.setdefaultencoding('UTF-8')
	if opt.query_search is None:
		cur.execute('SELECT * FROM {}'.format(tbl))
		dat = cur.fetchall()
		get = make_list_of_webfiles(dat, opt.format)
	else:
		if opt.query_type == 'exact':
			cur.execute('SELECT * FROM {} WHERE ({} = \'{}\' '
				'AND version_status = "latest")'.format(
				tbl, opt.query_feature, opt.query_search))
		elif opt.query_type == 'substring':
			cur.execute('SELECT * FROM {} WHERE ({} LIKE \'%{}%\' '
				'AND version_status = "latest")'.format(
				tbl, opt.query_feature, opt.query_search))
		dat = cur.fetchall()
		get = make_list_of_webfiles(dat, opt.format)
	sql_close(con, cur)

	# Optionally save assembly metadata to an output file
	if opt.info is not None:
		with open(os.path.abspath(os.path.expanduser(opt.info)), 'w') as o:
			for ln in generate_metadata(dat, opt.format):
				o.write('{}\n'.format(''.join(ln)))
		if not opt.quiet:
			sys.stderr.write('INFO: saved metadata for {} assemblies\n'.\
				format(len(get)))

	# Optionally save assembly URLs to an output file
	if opt.list is not None:
		write_file_list(os.path.abspath(os.path.expanduser(opt.list)), get,
			opt.quiet)
	
	# Optionally skip downloading assembly files from NCBI
	if opt.no_download:
		sys.exit(0)

	# Make generic download command before adding url to fetch
	## rsync <options> fetch-url output-dir
	## wget -P <output-dir> <options> fetch-url
	if opt.protocol == 'rsync':
		cmd = 'rsync --copy-links --times --human-readable '
	elif opt.protocol == 'wget':
		cmd = 'wget --directory-prefix={} '.format(out)
	if opt.connect_timeout:
		if opt.protocol == 'rsync':
			cmd += '--contimeout {} '.format(opt.connect_timeout)
		elif opt.protocol == 'wget':
			cmd += '--connect-timeout={} '.format(opt.connect_timeout)
	conflicts = {'rsync': 
					{'force':  '--ignore-times ',
					 'keep':   '--ignore-existing ',
					 'update': '--update '},
				 'wget':
					{'force':  '',
					 'keep':   '--no-clobber ',
					 'update': '--timestamping '}}
	cmd += conflicts[opt.protocol][opt.file_conflicts]
	if opt.quiet:
		cmd += '--quiet '
	else:
		sys.stderr.write('INFO: fetching {} assemblies into {}...\n'.format(
			len(get), out))

	# Now add file (and outdir) to each generic fetch command and download
	if opt.protocol == 'rsync':
		cmds = [cmd + x.replace('ftp:', 'rsync:') + ' ' + out for x in get]
	elif opt.protocol == 'wget':
		cmds = [cmd + x.replace('ftp:', 'https:') for x in get]
	
	# Fetch the first assembly file
	first_try = cmds.pop(0)
	if not fetch_asm(first_try):
		sys.stderr.write('ERROR: retrieval issue with {}. Consider an '
			'alternative connection protocol.\n'.format(opt.protocol))
		sys.exit(1)

	# Connection and download successful, so now multiprocess the rest
	pool = mp.Pool(processes=cpu_processes)
	try:	
		for output in pool.imap_unordered(fetch_asm, cmds):
			if not output:
				sys.stderr.write('ERROR: retrieval issue with {}. Consider '
					'an alternative connection protocol.\n'.format(
						opt.protocol))
				pool.terminate()
				pool.join()
				sys.exit(1)
		pool.close()
		pool.join()
	except KeyboardInterrupt:
		pool.terminate()
		pool.join()
		sys.exit(1)

if __name__ == '__main__':
	main()
