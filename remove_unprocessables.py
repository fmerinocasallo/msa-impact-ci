#!/usr/bin/env python

from __future__ import print_function

import argparse
import os.path

import Bio.Align
import Bio.SeqIO

def read_arguments () :
    """
    Read command line arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Remove unprocessable '
                                                  'sequences from the '
                                                  'given set of sequences')
                                    )
    parser.add_argument('-ss', '--set_sequences', required=True,
                        help=('filename of the file storing a set of sequences '
                              '(FASTA format) from which we want to remove '
                              'the unprocessable ones'))
    parser.add_argument('-up', '--unprocessables', required=True,
                        help=('filename of the file storing the unprocessable '
                              'sequences (FASTA format) to be remove'))

    return parser.parse_args()


def del_seqs_unprocessables ( filename_sequences, format_sequences,
                              filename_unprocessables, format_unprocessables ) :
    """
    Remove from the given set of sequences all those considered unprocessables.

    Arguments :
        - filename_sequences ( string )
            File storing a set of sequences
        - format_sequences ( string )
            File format of the set of sequences 
        - filename_unprocessables ( string )
            File storing unprocessable sequences
        - format_unprocessables ( string )
            File format of the unprocessable sequences
    """
    filename, extension = os.path.splitext(filename_sequences)
    filename_sequences_processables = filename + '_processables' + '.fasta'
    filename_sequences_unprocessables = filename + '_unprocessables' + '.fasta'

    with open(filename_sequences_processables, 'w') as handle_processables, \
         open(filename_sequences_unprocessables, 'w') as handle_unprocessables :

        num_processables = 0
        num_unprocessables = 0

        sequences = Bio.SeqIO.parse(filename_sequences, format_sequences)
        records_unprocessables = Bio.SeqIO.parse(filename_unprocessables,
                                                 format_unprocessables)

        for record_unprocessable in records_unprocessables :
            found = False
            for record in sequences :
                if record.id == record_unprocessable.id :
                    found = True
                
                    num_unprocessables += Bio.SeqIO.write(record,
                                                          handle_unprocessables,
                                                          'fasta')
                    break
                else : # record.id != record_unprocessable.id 
                    num_processables += Bio.SeqIO.write(record,
                                                        handle_processables,
                                                        'fasta')


        print('# Processables: {:d}'.format(num_processables))
        print('# Unprocessables: {:d}'.format(num_unprocessables))


def main () :
    args = read_arguments()

    del_seqs_unprocessables(args.set_sequences, 'fasta',
                            args.unprocessables, 'fasta')


if __name__ == "__main__" :
    main()

