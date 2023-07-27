#!/usr/bin/env python3


import os
import sys
from argparse import ArgumentParser

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parseArgs():
    parser = ArgumentParser(add_help=False,
        description='splits multi-FastA file into individual FastA files of'
        ' non-ambiguous segments of specified minimum length',
        epilog='useful to get individual FastA files for each segment aligned'
        ' to a reference and have a `samtools consensus` file to split up')
    req = parser.add_argument_group('Required')
    req.add_argument('-i', '--infile', required=True, metavar='FILE',
        help='input multi-FastA file')
    opt = parser.add_argument_group('Optional')
    opt.add_argument('-e', '--ext', metavar='STR', type=str, default='.fasta',
        help='extension name for each output FastA file [\'.fasta\']')
    opt.add_argument('-h', '--help', action='help',
        help='show this help message and exit')
    opt.add_argument('-m', '--min-length', metavar='INT', type=int, default=31,
        help='minimum length of non-ambiguous sequences to output [31]')
    opt.add_argument('-o', '--outdir', metavar='PATH',
        help='output directory path [dir of outfiles]')
    return parser.parse_args()

def save_segment_fasta(non_ambiguous_segment, min_length, rec_id, coordinate_start, coordinate_stop, ext, infile, outdir):
    basename_input_filename = os.path.basename(infile)
    if len(non_ambiguous_segment) >= min_length:
        out_rec = SeqRecord(
            Seq(non_ambiguous_segment), 
            id=rec_id + ':' + str(coordinate_start) + '-' + str(coordinate_stop),
            description=basename_input_filename
            )
        SeqIO.write(
            out_rec, os.path.join(
                outdir,
                rec_id  + ':' + str(coordinate_start) + '-' + str(coordinate_stop) + ext
                ),
           'fasta'
           )

def main():
    opt = parseArgs()

    # I/O handling
    infile = os.path.abspath(os.path.expanduser(opt.infile))
    if opt.outdir is None:
        outdir = os.path.dirname(infile)
    else:
        outdir = os.path.abspath(os.path.expanduser(opt.outdir))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Define non-ambiguous characters we want
    non_ambiguous_nucleotides = ['A', 'T', 'C', 'G']

    # Iterate over each individual sequence record
    for rec in SeqIO.parse(infile, 'fasta'):

        # Individual sequence handling
        coordinate_start, coordinate_stop = 1, 1
        non_ambiguous_segment = ''
        record_sequence  = str(rec.seq).upper()

        # Iterate over each site or position within a record or contig
        for site in record_sequence:
            len_record_sequence = len(record_sequence)
            print(site, 'rec len:', str(len_record_sequence), str(coordinate_start), str(coordinate_stop))

            if site in non_ambiguous_nucleotides:
                print('site passed', str(coordinate_start), str(coordinate_stop))

                # When the last site is non-ambiguous and exceeds length to store, save it
                if coordinate_stop == len_record_sequence:
                    print('final segment saved', str(coordinate_start), str(coordinate_stop))
                    save_segment_fasta(
                        non_ambiguous_segment,
                        opt.min_length,
                        rec.id,
                        coordinate_start,
                        coordinate_stop,
                        opt.ext,
                        opt.infile,
                        outdir
                        )
                non_ambiguous_segment += site
                coordinate_stop += 1

            # Handles the next site when non-ambiguous
            else:

                # Only when the previous several non-ambiguous sites meet or exceed length, save it
                print('segment saved', str(coordinate_start), str(coordinate_stop))
                save_segment_fasta(
                    non_ambiguous_segment,
                    opt.min_length,
                    rec.id,
                    coordinate_start,
                    coordinate_stop,
                    opt.ext,
                    opt.infile,
                    outdir
                    )

                # Regardless of storing the segment or not, move onto the next site
                non_ambiguous_segment = ''
                coordinate_stop += 1
                coordinate_start = coordinate_stop

if __name__ == '__main__':
    main()
