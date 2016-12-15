#!/usr/bin/env python


import os
from argparse import ArgumentParser
from shutil import copyfile
from string import replace
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
	parser = ArgumentParser(add_help=False, description='Converts VFDB '
		'protein FastA file into a Prokka-friendly FastA format '
		'to create a custom annotation database which can then be provided '
		'to Prokka as `--proteins <outfile>`. If the database has been '
		'parsed for a particular genus, naming the output file as the '
		'genus.faa and specifying the Prokka genus database path (-d) enables '
		'the `--usegenus` option in Prokka.')
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, help='input multi-FastA'
		' protein file from http://www.mgc.ac.cn/VFs/download.htm '
		'doi:10.1093/nar/gkv1239', metavar='FILE')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-d', '--setup-db', action=None, metavar='PATH',
		help='copy output file to specified database location, e.g., '
		'$HOME/.linuxbrew/opt/prokka/db/genus '
		'and format the database so that Prokka can now use it; requires '
		'prokka to be exported in $PATH')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', default=None, metavar='FILE',
		help='output file [./<infile>.faa]')
	return parser.parse_args()

def main():
	args = parseArgs()
	infile  = args.infile
	cfg_db  = args.setup_db
	if args.outfile:
		outfile = args.outfile
	else:
		outfile = os.path.join(os.getcwd(), os.path.basename(infile) + '.faa')
	if not os.path.exists(os.path.dirname(outfile)):
		os.mkdir(os.path.dirname(outfile))

	mfasta = SeqIO.parse(infile, 'fasta')
	revised_recs = []
	for rec in mfasta:
		rd = rec.description
		s  = rd.split(' ')
		accession   = s[0]
		gene_name   = s[1].lstrip('(').rstrip(')')
		description = replace(' '.join(s[2:]), '] [', ' from ')
		prot_seq    = str(rec.seq)
		defline     = '~~~{}~~~{}'.format(gene_name, description)
		revised_rec = SeqRecord(Seq(prot_seq, generic_protein), id=accession,
			description=defline)
		revised_recs.append(revised_rec)
	SeqIO.write(revised_recs, outfile, 'fasta')

	if cfg_db is not None:
		db = os.path.join(cfg_db, os.path.basename(outfile).split('.')[:-2][0])
		copyfile(outfile, db)
		os.system('prokka --setupdb')

if __name__ == '__main__':
	main()
