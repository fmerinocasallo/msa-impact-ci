#!/usr/bin/env python
"""Convert GB file to FASTA file."""

from __future__ import print_function

import argparse
from os.path import exists, splitext
import sys

import Bio.SeqIO

__author__ = 'Francisco Merino'
__copyright__ = 'Copyright 2013, Final Project'
__credits__ = ['Francisco Merino', 'Jorge Alvarez', 'Elvira Mayordomo']
__license__ = 'GPL'
__version__ = '3.0'
__maintainer__ = 'Francisco Merino'
__email__ = 'fmerino@unizar.es'
__status__ = 'Development'

def read_arguments():
    """
    Read command line's arguments
    """
    parser = argparse.ArgumentParser(description=('Convert GB file to FASTA '
                                                  'file')
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file storing the set of '
                              'sequences (GB format)'))
    parser.add_argument('-of', '--output_filename',
                        help=('filename of the file storing the set of '
                              'sequences (FASTA format)'))
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='give detailed information about the process')
    return parser.parse_args()

def main():
    """
    Convert a given set of mtDNA sequences stored in a GB file to a FASTA file.
    """

    # Deal with command line's args
    args = read_arguments()

    # Set the new FASTA file's filename which is going to store the set of
    # sequences
    if args.output_filename:
        output_filename = args.output_filename
    else:
        input_filename, input_extension = splitext(args.input_filename)
        output_filename = input_filename + '.fasta'

    if args.verbose:
        print('Starting convertion to FASTA file...', end='')
        sys.stdout.flush()

    Bio.SeqIO.convert(args.input_filename, 'gb', output_filename, 'fasta')

    if args.verbose:
        print('done')
        print('\nThe new FASTA file is: {:s}'.format(output_filename))

if __name__ == '__main__':
    main()

