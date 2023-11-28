#!/usr/bin/env python


import sys
from Bio import SeqIO


def usage():
	sys.exit('Quantifies assembly errors that are corrected with Illumina '
		'reads by Pilon. DOI: 10.1371/journal.pone.0112963.\n\n'
		'Usage: {} <pilon>.changes <pilon.fasta>.gbk 1>indel.tab 2>sub.tab\n'
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
	k = ['ins_CDS', 'ins_noCDS', 'ins2_CDS', 'ins2_noCDS',
	'del_CDS', 'del_noCDS', 'del2_CDS', 'del2_noCDS']
	d = dict.fromkeys(k, '')
	nucs = ['A', 'T', 'C', 'G']
	e = {'CDS': dict((n, dict((n, 0) for n in nucs)) for n in nucs),
	'noCDS': dict((n, dict((n, 0) for n in nucs)) for n in nucs)}
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
				if pos in idx[chrom]:
					d['ins_CDS'] += l[1]
				else:
					d['ins_noCDS'] += l[1]
			else:
				if pos in idx[chrom]:
					d['ins2_CDS'] += l[1]
				else:
					d['ins2_noCDS'] += l[1]
		elif '.' == l[1]:  # deletion
			if len(l[0]) == 1:
				if pos in idx[chrom]:
					d['del_CDS'] += l[0]
				else:
					d['del_noCDS'] += l[0]
			else:
				if pos in idx[chrom]:
					d['del2_CDS'] += l[0]
				else:
					d['del2_noCDS'] += l[0]
		else:  # substitution
			if pos in idx[chrom]:
				cds = 'CDS'
			else:
				cds = 'noCDS'
			x = ''  #e.g., 'TG' is T->G
			for y in (l[0], l[1]):
				for z in nucs:
					if y == z:
						x += z
			e[cds][x[0]][x[1]] += 1

# Report nucleotide frequency data for indels to stdout
head_indel = ('.\tCDS\tCDS\tCDS\tCDS\tnon-CDS\tnon-CDS\tnon-CDS\tnon-CDS\n'
'Nucleotide\tInsertion (1 bp)\tDeletion (1 bp)\t'
'Insertion (> 1 bp)\tDeletion (>1 bp)\tInsertion (1 bp)\t'
'Deletion (1 bp)\tInsertion (> 1 bp)\tDeletion (>1 bp)\n')
data_indel = ''
for j in ['A', 'T', 'C', 'G']:
	data_indel += j
	for k in ('CDS', 'noCDS'):
		for l in ('_', '2_'):
			for m in ('ins', 'del'):
				data_indel += '\t' + str((d[m + l + k]).count(j))
	data_indel += '\n'
print head_indel + data_indel

# Report nucleotide frequency data for substitutions to stderr
head_subst = ('Nucleotide (From)\tNucleotide (To)\tSubstitution in '
'Coding Sequence\tSubstitution in Non-Coding Sequence\n')
data_subst = ''
for s in nucs:
	for t in nucs:
		if s != t:
			c = str(e['CDS'][s][t]) + '\t' + str(e['noCDS'][s][t])
			data_subst += s + '\t' + t + '\t' + c + '\n'
sys.stderr.write(head_subst + data_subst)
