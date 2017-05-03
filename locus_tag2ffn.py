#!/usr/bin/env python


import os
import sys
from argparse import ArgumentParser
from BCBio import GFF  #pip install bcbio-gff
from BCBio.GFF import GFFParser
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
	parser = ArgumentParser(description='Extracts gene sequences from a '
		'General Feature Format (GFF) file given a locus_tag or '
		'comma-separated list of locus_tags. BED file of extraction(s) is '
		'printed to stdout.',
		epilog='a BED format of extracted gene(s) is printed to standard out')
	parser.add_argument('-i', '--infile', required=True,
		help='input General Feature Format file <.gff||.gff3>')
	parser.add_argument('-l', '--locus_tag', required=True, default=None,
		help='locus tag(s) to extract; comma-separate if >1')
	parser.add_argument('-o', '--outfile', required=False, default=None,
		help='[<locus_tag(s)>.ffn]')
	return parser.parse_args()

def extract_seq(beg, end, sns, seq, locus_name):
	''' returns a SeqRecord obj given a nucleotide sequence and record
	info to extract specific sites including: begin, stop, and strand sense '''
	extraction = seq[int(beg):int(end)]
	if sns < 0:
		return SeqRecord(Seq(extraction, generic_dna).reverse_complement(),
						id=locus_name, description='')
	else:
		return SeqRecord(Seq(extraction, generic_dna),
						id=locus_name, description='')

def locus_tag2rec_nfo(want, infile):
	''' returns a tuple containing record info 
	given a GFF file and locus_tag to locate '''
	parser = GFFParser()
	with open(infile, 'r') as gff:
		rec_nfo = (None, ) * 5
		for rec in parser.parse(gff, limit_info = dict(
			gff_type = ['gene', 'CDS', 'locus_tag', 'product'])):
			for feat in rec.features:
				if feat.type == 'gene' and 'locus_tag' in feat.qualifiers:
					locus_tag = feat.qualifiers.get('locus_tag', None)
					if want in locus_tag:
						rec_nfo = (feat.location.start.position,
							feat.location.end.position, feat.strand,
							rec.id, rec.seq.tostring())
	return rec_nfo

def main():
	opt = parseArgs()
	ifh = opt.infile
	if ',' in opt.locus_tag:
		loc_wnt = opt.locus_tag.split(',')
	else:
		loc_wnt = opt.locus_tag

	if opt.outfile is not None:
		ofh = os.path.abspath(os.path.expanduser(opt.outfile))
	else:
		ofh = opt.locus_tag + '.ffn'

	# Extract one or more locus tags
	if isinstance(loc_wnt, list):
		seq_rec = []
		for loc in loc_wnt:
			rec_nfo = locus_tag2rec_nfo(loc, ifh)
			if all(rec_nfo):
				(beg, end, sns, ctg_id, ctg_seq) = rec_nfo
				if sns < 0:
					sym = '-'
				else:
					sym = '+'
				print '{}\t{}\t{}\t{}\t0\t{}'.format(ctg_id, beg, end, loc, sym)
				sr = extract_seq(beg, end, sns, ctg_seq, loc)
				seq_rec.append(sr)
			else:
				sys.exit('ERROR: absent {}'.format(loc))
	else:
		rec_nfo = locus_tag2rec_nfo(loc_wnt, ifh)
		if all(rec_nfo):
			(beg, end, sns, ctg_id, ctg_seq) = rec_nfo
			if sns < 0:
				sym = '-'
			else:
				sym = '+'
			print '{}\t{}\t{}\t{}\t0\t{}'.format(ctg_id, beg, end, loc_wnt, sym)
			seq_rec = extract_seq(beg, end, sns, ctg_seq, loc_wnt)
		else:
			sys.exit('ERROR: absent {}'.format(loc_wnt))

	# Write one or more extracted genes to FastA
	SeqIO.write(seq_rec, ofh, 'fasta')

if __name__ == '__main__':
	main()
