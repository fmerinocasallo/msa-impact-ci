#!/usr/bin/env python

from __future__ import print_function

import argparse
import collections
import os.path
import re

import Bio.SeqIO

def read_arguments () :
    """
    Read command line's arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Analyze gene fragments length '
                                                  'by haplogroups')
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file (GenBank format) which '
                              'stores the set of sequences'))

    return parser.parse_args()


def group_records_by_haplogroups ( input_filename, haplogroups ) :
    """
    Group records by haplogroups.

    Arguments :
        input_filename ( string )
            File storing a set of sequences ( GenBank format).
        haplogroups ( dictionary )
            Dictionary with haplogroups as keys and SeqRecords as values.
    """
    regex_main_haplogroups = '({:s})'.format('|'.join(haplogroups))

    for record in Bio.SeqIO.parse(input_filename, 'gb') :
        # get record haplogroup
        if record.features :
            record_haplogroups = \
                record.features[0].qualifiers.get('haplogroup', None)
            if record_haplogroups :
                record_haplogroup = re.split(regex_main_haplogroups,
                                             record_haplogroups[0])[1]

                haplogroups[record_haplogroup].append(record)
            else :
                print('WARNING!: {:s},'.format(record.id))
        else :
            print('WARNING!: {:s},'.format(record.id))


def store_records_by_haplogroups ( input_filename, haplogroups ) :
    """
    Store gene fragments in different files regarding its haplogroup.

    Arguments :
        input_filename ( string )
            File storing a set of sequences ( GenBank format).
        haplogroups ( dictionary )
            Dictionary with haplogroups as keys and SeqRecords as values.
    """
    for haplogroup, gene_fragments in haplogroups.items() :
        if ( haplogroup == 'mt-MRCA' ) or ( not gene_fragments ) :
            break
        else :
            filename, extension = os.path.splitext(input_filename)
            output_filename = filename + '_' + haplogroup + extension

            Bio.SeqIO.write(gene_fragments, output_filename, 'gb')


def
    """
    """
    
    filenames = "hmtDNA_2015_11_09_COX1_with_haplogroups_aligned_with_*_X*"
    command = ['diff', '-q',
               '--from-file', filenames,
               ';', 'echo', '$?']

    

def main () :
    args = read_arguments()

    tools = ['clustalo', 'mafft_auto', 'mafft_parttree']
    haplogroups = ['L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6',  'M',  'C',  'E',
                    'G',  'Q',  'Z',  'D',  'N',  'A',  'I',  'O',  'S',  'W',
                    'X',  'Y', 'R0', 'B6',  'F',  'J',  'P',  'K',  'B', 'HV',
                    'T',  'R',  'H',  'V',  'U',
                    '&',  'mt-MRCA']

    group_records_by_haplogroups(args.input_filename, haplogroups)

    store_records_by_haplogroups(args.input_filename, haplogroups)


if __name__ == "__main__" :
    main()

data/reports/hmtDNA_2015_11_09_COX1_with_haplogroups_aligned_with_clustalo_without_gaps_A_unweighted_entropy_less_1.0_detailed.txt
