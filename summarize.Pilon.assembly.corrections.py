#!/usr/bin/env python


import sys
from Bio import SeqIO

def usage():
	sys.exit('Quantifies assembly errors that are corrected with Illumina '
		'reads by Pilon. DOI: 10.1371/journal.pone.0112963.\n\n'
		'Usage: {} <pilon>.changes <pilon.fasta>.gbk\n'
		'where the GenBank file has identical LOCUS ID tags as the '
		'FastA given to Pilon. Note: `sed -i \'s/_pilon//g\' <pilon.fasta>` '
		'is likely required prior to Prokka annotation.'.format(sys.argv[0]))

if len(sys.argv) != 3:
	usage()
elif any(h == sys.argv[1] for h in ['-h', '-help', '--help']):
	usage()

idx = {}
for rec in SeqIO.parse(open(sys.argv[2]), 'genbank'):
	contig_idx = set()
	for (i, feat) in enumerate(rec.features):
		if feat.type == 'CDS':
			loc = str(feat.location)
			if loc.endswith('(-)') or loc.endswith('(+)'):
				b,e = loc.lstrip('[').rstrip('](-)').rstrip('](+)').split(':')
				contig_idx.update(range(int(b), int(e) + 1))
			else:
				sys.exit('ERROR')
	idx[rec.id] = contig_idx

with open(sys.argv[1]) as ifh:
	# So many empties!
	i, cnt_sub, cnt_sub_cds = 0, 0, 0
	cnt_ins, cnt_ins2, cnt_del, cnt_del2 = 0, 0, 0, 0
	cnt_ins_cds, cnt_del_cds, cnt_ins2_cds, cnt_del2_cds = 0, 0, 0, 0
	sub_from, sub_to, ins, ins2, deL, del2 = '', '', '', '', '', ''
	for ln in ifh:
		# Verify FastA headers given to Pilon corresponds to GenBank provided
		chrom = str(ln.split()[0].split(':')[0])
		pos   = int(ln.split()[0].split(':')[1].split('-')[0])
		if chrom not in idx:
			sys.exit('ERROR: {} not in GenBank entry as a LOCUS identifier.'
				' The FastA deflines given to Pilon must match with GenBank'
				' LOCUS ID entries.'.format(chrom))

		# Count and classify each error type
		l = ln.split()[2:4]
		if '.' == l[0]:  # insertion
			if len(l[1]) == 1:
				cnt_ins += 1
				ins     += l[1]
				if pos in idx[chrom]:
					cnt_ins_cds += 1
			else:
				cnt_ins2 += 1
				ins2     += l[1]
				if pos in idx[chrom]:
					cnt_ins2_cds += 1
		elif '.' == l[1]:  # deletion
			if len(l[0]) == 1:
				cnt_del += 1
				deL     += l[0]
				if pos in idx[chrom]:
					cnt_del_cds += 1
			else:
				cnt_del2 += 1
				del2     += l[0]
				if pos in idx[chrom]:
					cnt_del2_cds += 1
		else:  # substitution
			cnt_sub  += 1
			sub_from += l[0]
			sub_to   += l[1]
			if pos in idx[chrom]:
				cnt_sub_cds += 1
		i += 1

# Report error counts
for s, t in [('substitu', cnt_sub),
(' 1 bp inser', cnt_ins), ('>1 bp inser', cnt_ins2),
(' 1 bp dele', cnt_del), ('>1 bp dele', cnt_del2)]:
	print 'Total {}tions: {}'.format(s, t)

# Report errors being within or outside of coding sequences
for s, t, u in [('sub', cnt_sub, cnt_sub_cds),
('1 bp ins', cnt_ins, cnt_ins_cds), ('>1 bp ins', cnt_ins2, cnt_ins2_cds),
('1 bp del', cnt_del, cnt_del_cds), ('>1 bp del', cnt_del2, cnt_del2_cds)]:
	for m, n in [(t - u, 'non-'), (u, '')]:
		print '   {} in {}coding sequences: {}'.format(s, n, m)

# Report nucleotide frequency data
for s, t, u in [('sub', 'from', sub_from), ('sub', 'to', sub_to),
('1 bp', 'ins', ins), ('>1 bp', 'ins', ins2),
('1 bp', 'del', deL), ('>1 bp', 'del', del2)]:
	for n in ['A', 'T', 'C', 'G']:
		print '      {} {} {}: {}'.format(s, t, n, u.count(n))
