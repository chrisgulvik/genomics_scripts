#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from csv import reader
from shutil import rmtree
from tempfile import mkdtemp
import sqlite3

def parseArgs():
	parser = ArgumentParser(description='Creates (or updates) a SQLite3 '
		'database of genome assembly information from NCBI\'s FTP site '
		'(or a given summary file)', add_help=False,
		epilog='The output database file is especially useful for retrieving '
		'genome assembly files with the aid of NCBI.asm.fetch.py')
	req = parser.add_argument_group('Required')
	req.add_argument('-o', '--output', required=True, metavar='FILE',
		help='output SQLite3 database file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--compress', action='store_true', default=False,
		help='compress SQLite3 output database; especially useful when the '
		'infile database entries overwrite fetched or provided summary data '
		'or contains many different additional entries [off]')
	opt.add_argument('-d', '--db', choices=['genbank', 'refseq'],
		default='refseq', help='database to summarize from NCBI [refseq]')
	opt.add_argument('-f', '--force', action='store_true', default=False,
		help='when input SQL entry has a conflicting assembly_accession '
		'entry with the summary data fetched or provided, overwrite/replace '
		'entry in the summary data with the infile database entry [off]')
	opt.add_argument('-g', '--group', choices=['all', 'archaea', 'bacteria',
		'fungi', 'invertebrate', 'metagenomes', 'other', 'plant', 'protozoa',
		'vertebrate_mammalian', 'vertebrate_other', 'viral'],
		metavar='{archaea, bacteria, viral}',
		default='bacteria', help='group of assembly data to fetch [bacteria]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-i', '--infile', metavar='FILE',
		help='input SQLite3 database file which is appended to the summary '
		'data output; table name must be \'<--db>_asm_summary\'; see '
		'--force for handling entry conflicts')
	opt.add_argument('-q', '--quiet', action='store_true', default=False,
		help='hide messages when entries are updated or added [off]')
	opt.add_argument('-r', '--rows', metavar='INT', type=int, default=200000,
		help='number of data summary rows to simultaneously add to SQLite '
		'database at once; useful for speed optimization; default takes 5 '
		'sec for refseq whereas naive insertion of 1 takes 7 min [200000]')
	opt.add_argument('-s', '--summary-data', metavar='FILE',
		help='tab-delimited local file of assembly summary data (22-column '
		'format as NCBI has specified in ftp://ftp.ncbi.nlm.nih.gov/genomes/'
		'ASSEMBLY_REPORTS/README_assembly_summary.txt) rather than directly '
		'fetching from NCBI\'s FTP site; especially useful if FTP is '
		'blocked or creating custom database')
	opt.add_argument('-u', '--url', metavar='STR', type=str, default=None,
		help='web address to file of assembly summary data (22-column); '
		'for bypassing NCBI access (e.g., EMBL-EBI) without locally storing '
		'the file')
	return parser.parse_args()

def sql_open(sql_out):
	try:
		con = sqlite3.connect(sql_out)
		cur = con.cursor()
		return con, cur
	except Error as e:
		sys.stderr.write('ERROR: {}\n'.format(e))
		sys.exit(1)

def sql_compress(cur):
	cur.execute('VACUUM')

def sql_close(con, cur):
	cur.close()
	con.close()

def sql_create_table(cur, tbl):
	cur.execute('''
		CREATE TABLE IF NOT EXISTS {} (
		assembly_accession TEXT UNIQUE NOT NULL,
		bioproject TEXT,
		biosample TEXT,
		wgs_master TEXT,
		refseq_category TEXT,
		taxid INT,
		species_taxid INT,
		organism_name TEXT,
		infraspecific_name TEXT,
		isolate TEXT,
		version_status TEXT,
		assembly_level TEXT,
		release_type TEXT,
		genome_rep TEXT,
		seq_rel_date TEXT,
		asm_name TEXT,
		submitter TEXT,
		gbrs_paired_asm TEXT,
		paired_asm_comp TEXT,
		ftp_path TEXT NOT NULL,
		excluded_from_refseq TEXT,
		relation_to_type_material TEXT
	)'''.format(tbl))

def sql_data_entry(con, cur, tbl, dat, bulk=False, overwrite=False):
	if bulk and not overwrite:
		cur.executemany('INSERT INTO {} VALUES({})'.format(
			tbl, ('?,' * 22)[:-1]), dat)
	elif bulk and overwrite:
		cur.executemany('REPLACE INTO {} VALUES({})'.format(
			tbl, ('?,' * 22)[:-1]), dat)
	elif not bulk and overwrite:
		cur.execute('REPLACE INTO {} VALUES({})'.format(
			tbl, ('?,' * 22)[:-1]), dat)
	else:
		cur.execute('INSERT INTO {} VALUES({})'.format(
			tbl, ('?,' * 22)[:-1]), dat)
	con.commit()

def sql_query(cur, tbl, query, quiet):
	cur.execute('SELECT * FROM {} WHERE assembly_accession=?'.format(
		tbl), (query,))
	fnd = cur.fetchall()
	if len(fnd) >= 1:
		if not quiet:
			sys.stderr.write('INFO: {} already found in db\n'.format(query))
		if len(fnd) == 1:
			return fnd[0]
		elif len(fnd) > 1:
			if not quiet:
				sys.stderr.write('WARNING: expect 1 but found {} entries for '
					'{}. Using first value {}...\n'.format(
						len(fnd), query, fnd[0]))
			return fnd[0]
	else:
		return False

def chunk_data(data, rows_per_chunk=5000):
	chunk = []
	for i, row in enumerate(data):
		if i % rows_per_chunk == 0 and i > 0:
			yield chunk
			del chunk[:]
		chunk.append([str(s).decode('UTF-8') for s in row])
	yield chunk

def main():
	opt = parseArgs()
	ofh = os.path.abspath(os.path.expanduser(opt.output))
	if opt.infile:
		ifh = os.path.abspath(os.path.expanduser(opt.infile))
		if ifh == ofh:
			sys.stderr.write('ERROR: in and out files cannot be the same\n')
			sys.exit(1)

	# Prepare or fetch summary file of assembly data
	if opt.summary_data:
		summary_file = os.path.abspath(os.path.expanduser(opt.summary_data))
		with open(summary_file) as f:
			ln = next(f)
			num_col = len(ln.split('\t'))
		if num_col != 22:
			sys.stderr.write('ERROR: expect 22 data columns from input '
				'summary file; contains {} columns\n'.format(num_col))
			sys.exit(1)
	else:
		tmp = mkdtemp()
		tfh = os.path.join(tmp, 'asm-summary.tab')
		if opt.url is not None:
			url = opt.url
		elif opt.group == 'all':
			url = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/'
				'assembly_summary_{}.txt'.format(opt.db))
		else:
			url = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/{}/{}/'
				'assembly_summary.txt'.format(opt.db, opt.group))
		os.system('curl --max-time 300 --silent --fail --show-error '
			'{} -o {}'.format(url, tfh))
		if os.path.exists(tfh) and os.path.getsize(tfh) > 0:
			summary_file = tfh
		else:
			sys.stderr.write('ERROR: unable to fetch assembly summary file '
				'from NCBI\'s FTP site\n')
			rmtree(tmp)
			sys.exit(1)

	# Output sqlite3 db from fetched or provided assembly data
	tbl = 'asm_summary'
	con, cur = sql_open(ofh)
	sql_create_table(cur, tbl)
	with open(summary_file) as f:
		r = reader((row for row in f if not row.startswith('#')),
			delimiter='\t')
		for chunk in chunk_data(r, rows_per_chunk=opt.rows):
			sql_data_entry(con, cur, tbl, chunk, bulk=True)

	# Handle comparisons of input to output sqlite3 dbs
	if opt.infile:
		sql_close(con, cur)
		con, cur = sql_open(ifh)
		cur.execute('SELECT * FROM {}'.format(tbl))
		input_data = cur.fetchall()
		sql_close(con, cur)
		con, cur = sql_open(ofh)
		reload(sys)
		sys.setdefaultencoding('UTF-8')
		if opt.force:
			for chunk in chunk_data(input_data, rows_per_chunk=opt.rows):
				sql_data_entry(con, cur, tbl, chunk, bulk=True, overwrite=True)
		else:
			for entry in input_data:
				a = sql_query(cur, tbl, str(entry[0]), opt.quiet)
				if not a:
					sql_data_entry(con, cur, tbl,
						[str(s).decode('UTF-8') for s in entry])

	# Finish and cleanup
	if opt.compress:
		sql_compress(cur)
	sql_close(con, cur)
	if not opt.summary_data:
		rmtree(tmp)

if __name__ == '__main__':
	main()
