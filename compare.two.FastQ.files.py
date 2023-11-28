#!/usr/bin/env python3


import gzip
import os
import shutil
import sys
from argparse import ArgumentParser
from tempfile import mkdtemp

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import numpy as np
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def parse_args():
    parser = ArgumentParser(
        description='Compares two FastQ files and reports differences between'
        ' the two. Initially developed to audit mutated FastQ files and'
        ' report the number of shared and unshared sequence headers as well'
        ' as the SNP frequency of artificially mutated sequences.',
        epilog='NOTE: order between the two files is expected to be sorted '
        'based on the header names with no missing sequences between the two',
        add_help=False,
    )
    req = parser.add_argument_group('Required')
    req.add_argument(
        '-a',
        '--after',
        required=True,
        metavar='FILE',
        help='modified or processed FastQ sequence file'
        ' , optionally gunzip compressed',
    )
    req.add_argument(
        '-b',
        '--before',
        required=True,
        metavar='FILE',
        help='initial FastQ sequence file , optionally gunzip compressed',
    )
    opt = parser.add_argument_group('Optional')
    opt.add_argument(
        '-o',
        '--outdir',
        metavar='PATH',
        type=str,
        default=None,
        help='outpath to save two files (after and before) with differing'
        ' sequence records [None]',
    )
    opt.add_argument(
        '-h', '--help', action='help', help='show this help message and exit'
    )
    return parser.parse_args()


def decompress_file(infile, outdir):
    # Create a new output file name
    uncompressed_file = os.path.basename(infile).rstrip('.gz')
    outfile = os.path.join(outdir, uncompressed_file)

    # Decompress the gunzipped file
    with gzip.open(infile, 'rb') as ifh, open(outfile, 'wb') as ofh:
        shutil.copyfileobj(ifh, ofh)

    # Send the full path of the decompressed file back
    return outfile


def compare_two_fastq_files(after, before, outdir):
    comparisons_counted_snps = []
    sequence_records_with_different_sequences_after = []
    sequence_records_with_different_sequences_before = []
    for (
        title_after,
        sequence_after,
        quality_after,
        title_before,
        sequence_before,
        quality_before,
    ) in zip(FastqGeneralIterator(after), FastqGeneralIterator(before)):
        # Verify these two records have the same sequence header
        if str(title_after) != str(title_before):
            sys.stderr.write(
                'INFO: skipping {} in {} and {} in {} due to different'
                ' sequence headers.'.format(
                    title_after,
                    os.path.basename(after),
                    title_before,
                    os.path.basename(before),
                )
            )
            continue

        # Verify these two records have the same sequence length
        if len(sequence_after) != len(sequence_before):
            sys.stderr.write(
                'INFO: skipping {} in {} and {} in {} due to different'
                ' sequence lengths.'.format(
                    title_after,
                    os.path.basename(after),
                    title_before,
                    os.path.basename(before),
                )
            )
            continue

        # Compare the sequence composition between the two for SNPs
        cnt_SNPs = 0
        for nucleotide_after, nucleotide_before in zip(
            sequence_after, sequence_before
        ):
            if nucleotide_after != nucleotide_before:
                cnt_SNPs += 1
        comparisons_counted_snps.append(cnt_SNPs)

        if outdir is not None:
            # Construct new FastQ sequence record objects
            fastq_record_string = '@{}\n{}\n+\n{}\n'.format(
                title_after, sequence_after, quality_after
            )
            new_rec = SeqIO.read(StringIO(fastq_record_string), 'fastq')
            sequence_records_with_different_sequences_after.append(new_rec)
            fastq_record_string = '@{}\n{}\n+\n{}\n'.format(
                title_before, sequence_before, quality_before
            )
            new_rec = SeqIO.read(StringIO(fastq_record_string), 'fastq')
            sequence_records_with_different_sequences_before.append(new_rec)

    # Optionally, write only mutated output
    if outdir is not None:
        out_after = os.path.join(outdir, 'different_sequences.after.fastq')
        out_before = os.path.join(outdir, 'different_sequences.before.fastq')
        SeqIO.write(
            sequence_records_with_different_sequences_after,
            out_after,
            'fastq',
        )
        SeqIO.write(
            sequence_records_with_different_sequences_before,
            out_before,
            'fastq',
        )
        print(
            'INFO: saved all {} sequence records with different'
            ' sequences'.format(
                len(sequence_records_with_different_sequences_after)
            )
        )

    return comparisons_counted_snps


def calculate_stats(list_of_integers):
    stats = {}
    stats['q1'] = int(
        np.percentile(list_of_integers, 25, interpolation='midpoint')
    )
    stats['median'] = int(np.median(list_of_integers))
    stats['q3'] = int(
        np.percentile(list_of_integers, 75, interpolation='midpoint')
    )
    stats['mean'] = int(np.mean(list_of_integers))
    stats['stdev'] = int(np.std(list_of_integers))
    stats['total_sequence_records'] = len(list_of_integers)
    stats['total_SNPs'] = int(np.sum(list_of_integers))
    stats['most_SNPs_in_a_record'] = int(np.amax(list_of_integers))
    stats['fewest_SNPs_in_a_record'] = int(np.amin(list_of_integers))
    return stats


def main():
    opt = parse_args()

    # I/O handling
    after = os.path.realpath(os.path.expanduser(opt.after))
    before = os.path.realpath(os.path.expanduser(opt.before))
    if opt.outdir is None:
        outdir = None
    else:
        outdir = os.path.realpath(os.path.expanduser(opt.outdir))

    # Auto-handle gunzip compressed input
    tmp = mkdtemp()
    if after.endswith('.gz'):
        after = decompress_file(after, tmp)
    if before.endswith('.gz'):
        before = decompress_file(before, tmp)

    # Compare the two sequence files
    comparisons_counted_snps = compare_two_fastq_files(after, before, outdir)

    # Cleanup
    shutil.rmtree(tmp)

    # Report statistics on SNPs observed between the two sequence sets
    stats = calculate_stats(comparisons_counted_snps)
    for k, v in stats.items():
        print('{} = {}'.format(k.replace('_', ' '), v))


if __name__ == '__main__':
    main()
