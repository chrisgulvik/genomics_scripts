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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def parse_args():
    parser = ArgumentParser(
        description='Introduce SNPs into sequence records',
        epilog='NOTE: More mutations than the length of a sequence is'
        ' supported. Each mutation position is chosen at random, and it is'
        ' possible 2 mutations would occur at the same site.',
        add_help=False,
    )
    req = parser.add_argument_group('Required')
    req.add_argument(
        '-i',
        '--infile',
        required=True,
        metavar='FILE',
        help='input sequence file '
        '(FastA or FastQ, optionally gunzip compressed)',
    )
    opt = parser.add_argument_group('Input Options')
    opt.add_argument(
        '-k',
        '--keep-lowercase',
        action='store_true',
        default=False,
        help='keep mutated nucleotides lowercase [off]',
    )
    opt.add_argument(
        '-m',
        '--mutations-per-record',
        metavar='INT',
        type=require_int_nonnegative,
        default=1,
        help='number of mutation events to introduce for each sequence record'
        ' [%(default)s]',
    )
    opt.add_argument(
        '-o',
        '--outfile',
        metavar='FILE',
        type=str,
        default=None,
        help='output file [stdout]',
    )
    opt.add_argument(
        '-h',
        '--help',
        action='help',
        help='show this help message and exit'
    )
    return parser.parse_args()


def require_int_nonnegative(x):
    try:
        if int(x) < 0 or '.' in str(x):
            sys.stderr.write(
                'ERROR: {} must be a non-negative integer\n'.format(x)
            )
            sys.exit(1)
    except ValueError:
        sys.stderr.write('ERROR: {} must be an integer\n'.format(x))
        sys.exit(1)
    return int(x)


def decompress_file(infile, outdir):
    # Create a new output file name
    uncompressed_file = os.path.basename(infile).rstrip('.gz')
    outfile = os.path.join(outdir, uncompressed_file)

    # Decompress the gunzipped file
    with gzip.open(infile, 'rb') as ifh, open(outfile, 'wb') as ofh:
        shutil.copyfileobj(ifh, ofh)

    # Send the full path of the decompressed file back
    return outfile


def mutate_nucleotide_site(input_nucleotide, bool_keep_lowercase):
    nucleotides = ['a', 't', 'c', 'g']

    # Remove the input nucleotide from the list to guarantee we create a SNP
    # NOTE: even an 'N' or other non-ATCG nucleotide will be mutated
    if input_nucleotide.lower() in nucleotides:
        nucleotides.remove(input_nucleotide.lower())

    # Randomly choose a new (different) nucleotide
    random_value = np.random.randint(low=0, high=len(nucleotides))
    new_nucleotide = nucleotides[random_value]

    # Choose upper or lowercase for the mutation character
    if not bool_keep_lowercase:
        new_nucleotide = new_nucleotide.upper()

    return new_nucleotide


def mutate_nucleotide_string(input_string, position, bool_keep_lowercase):
    # Extract the individual nucleotide to be mutated
    old_nucleotide = input_string[position]

    # Randomly choose a different nucleotide
    new_nucleotide = mutate_nucleotide_site(
        old_nucleotide, bool_keep_lowercase
    )

    # Form the new full length nucleotide string
    mutated_nucleotide_string = (
        input_string[:position] + new_nucleotide + input_string[position+1:]
    )

    return mutated_nucleotide_string


def mutate_fastq(infile, mutations_per_record, bool_keep_lowercase, outfile):
    new_sequence_records = []
    for title, sequence, quality in FastqGeneralIterator(infile):
        seq = str(sequence)

        # Determine random positions within the nucleotide string to mutate
        random_positions = list(
            np.random.randint(
                low=0, high=len(sequence), size=mutations_per_record
            )
        )

        # Mutate the sequence; handles 1 or more mutations per record
        for pos in random_positions:
            seq = mutate_nucleotide_string(seq, pos, bool_keep_lowercase)

        # Verify we truly have mutated the sequence
        if seq == sequence:
            sys.stderr.write('ERROR: {} sequence not mutated\n'.format(seq))
            sys.exit(1)

        # Construct a new FastQ sequence record object
        fastq_record_string = '@{}\n{}\n+\n{}\n'.format(title, seq, quality)
        new_rec = SeqIO.read(StringIO(fastq_record_string), 'fastq')
        new_sequence_records.append(new_rec)

    # Write the mutated output and report some general info
    SeqIO.write(new_sequence_records, outfile, 'fastq')
    print(
        'INFO: mutated all {} sequence records'.format(
            len(new_sequence_records)
        )
    )


def mutate_fasta(infile, mutations_per_record, bool_keep_lowercase, outfile):
    new_sequence_records = []
    for record in SeqIO.parse(infile, 'fasta'):
        seq = str(record.seq)

        # Determine random positions within the nucleotide string to mutate
        random_positions = list(
            np.random.randint(low=0, high=len(seq), size=mutations_per_record)
        )

        # Mutate the sequence; handles 1 or more mutations per record
        for pos in random_positions:
            seq = mutate_nucleotide_string(seq, pos, bool_keep_lowercase)

        # Verify we truly have mutated the sequence
        if seq == record.seq:
            sys.stderr.write('ERROR: {} sequence not mutated\n'.format(seq))
            sys.exit(1)

        # Construct a new FastA sequence record object
        new_rec = SeqRecord(
            id=record.description, seq=Seq(seq), description='', name=''
        )
        new_sequence_records.append(new_rec)

    # Write the mutated output and report some general info
    SeqIO.write(new_sequence_records, outfile, 'fasta')
    print(
        'INFO: mutated all {} sequence records'.format(
            len(new_sequence_records)
        )
    )


def main():
    opt = parse_args()

    # I/O handling
    mutations_per_record = opt.mutations_per_record
    infile = os.path.realpath(os.path.expanduser(opt.infile))
    if opt.outfile is None:
        outfile = sys.stdout
    else:
        outfile = os.path.realpath(os.path.expanduser(opt.outfile))
    bool_keep_lowercase = opt.keep_lowercase

    # Auto-handle gunzip compressed input
    tmp = mkdtemp()
    if infile.endswith('.gz'):
        infile = decompress_file(infile, tmp)

    # Iteratate across each sequence record and mutate
    if infile.lower().endswith(('.fastq', '.fq', '.fsq')):
        mutate_fastq(
            infile, mutations_per_record, bool_keep_lowercase, outfile
        )
    elif infile.lower().endswith(('.fasta', '.fa', '.fsa', '.fas')):
        mutate_fasta(
            infile, mutations_per_record, bool_keep_lowercase, outfile
        )
    else:
        sys.stderr.write(
            'ERROR: unable to detect file extension as FastA or FastQ'
            ' format\n'
        )
        sys.exit(1)

    # Cleanup
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
