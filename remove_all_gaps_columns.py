#!/usr/bin/env python

from __future__ import print_function

import argparse
import copy
import os.path

import Bio.SeqIO

def read_arguments () :
    """
    Read command line's arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Remove those columns with '
                                                  'all gaps in them')
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file (FASTA format) which '
                              'stores the aligned set of sequences'))

    return parser.parse_args()


def remove_gaps ( input_filename ) :
    """
    Remove alignment columns with a gap in every sequence and store on disk this
    modified alignment.

    Arguments :
        input_filename ( string )
            Aligned set of sequences (FASTA format).
    """
    # We are trying to keep mem consumption as low as possible thus we avoid
    # using Bio.AlignIO.read() method. Instead, we are reading one sequence
    # at a time which will increase execution time.
    ranges = []
    base_record = next(Bio.SeqIO.parse(input_filename, 'fasta'))
    len_alignment = len(base_record)
    start = -1
    # Check if every sequence has gaps in those columns where base_record has
    # them. We'll remove only those columns
    for col, bp in enumerate(base_record.seq) :
        if bp == '-' :
            all_gaps = True
            for record in Bio.SeqIO.parse(input_filename, 'fasta') :
                if record.seq[col] != '-' :
                    all_gaps = False

                    break

            if all_gaps :
                end = col
                ranges.append((start, end))

                start = -1
            else :
                start = col
        else :
            if start == -1 :
                start = col
 
    for record in Bio.SeqIO.parse(input_filename, 'fasta') :
        edited_record = copy.copy(record)
        edited_record.seq = ""
        for start, end in ranges :
            edited_record.seq += record.seq[start:end]

        yield edited_record

    
def main () :
    args = read_arguments()

    # Remove only those alignment columns which every sequence has gaps in it
    # and store the modified alignment on disk
    filename, extension = os.path.splitext(output_filename)
    output_filename = filename + '_without_gaps' + extension
    Bio.SeqIO.write(remove_gaps(input_filename), output_filename, 'fasta')


if __name__ == "__main__" :
    main()

