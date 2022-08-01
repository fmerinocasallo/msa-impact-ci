#!/usr/bin/env python

from __future__ import print_function

import argparse
import collections
import fnmatch
import os
import os.path
import re
import subprocess

import Bio.SeqIO

def read_arguments () :
    """
    Read command line's arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Generate conservation '
                                                  'reports from given sets '
                                                  'of sequences')
                                    )
    # parser.add_argument('-if', '--input_filename', required=True,
    #                     help=('filename of the file (GenBank format) which '
    #                           'stores the set of sequences'))

    return parser.parse_args()


def main () :
    genes = { 'tRNA-Phe': (  577,   647),  'rRNA-12S': (  648,  1601),
              'tRNA-Val': ( 1602,  1670),  'rRNA-16S': ( 1671,  3229),
             'tRNA-Leu1': ( 3230,  3304),       'ND1': ( 3307,  4262),
              'tRNA-Ile': ( 4263,  4331),  'tRNA-Gln': ( 4329,  4400),
              'tRNA-Met': ( 4402,  4469),       'ND2': ( 4470,  5511),
              'tRNA-Trp': ( 5512,  5579),  'tRNA-Ala': ( 5587,  5655),
              'tRNA-Asn': ( 5657,  5729),  'tRNA-Cys': ( 5761,  5826),
              'tRNA-Tyr': ( 5826,  5891),      'COX1': ( 5904,  7445),
             'tRNA-Ser1': ( 7446,  7514),  'tRNA-Asp': ( 7518,  7585),
                  'COX2': ( 7586,  8269),  'tRNA-Lys': ( 8295,  8364),
                  'ATP8': ( 8366,  8572),      'ATP6': ( 8527,  9207),
                  'COX3': ( 9207,  9990),  'tRNA-Gly': ( 9991, 10058),
                   'ND3': (10059, 10404),  'tRNA-Arg': (10405, 10469),
                  'ND4L': (10470, 10766),       'ND4': (10760, 12137),
              'tRNA-His': (12138, 12206), 'tRNA-Ser2': (12207, 12265),
             'tRNA-Leu2': (12266, 12336),       'ND5': (12337, 14148),
                   'ND6': (14149, 14673),  'tRNA-Glu': (14674, 14742),
                  'CYTB': (14747, 15887),  'tRNA-Thr': (15888, 15953),
              'tRNA-Pro': (15956, 16023),    'D-loop': (16024,   576)}

    fms = ['unweighted', 'weighted']
    cms = ['entropy', 'variance']

    path = 'data/dbs/global/'
    args = read_arguments()

    for gene_id in genes :
        regex = 'hmtDNA_2014_12_04_aligned_mafft_*_' + gene_id + '.fasta'
        start_column = genes[gene_id][0]
        for filename in os.listdir(path) :
            if ( fnmatch.fnmatch(filename, regex) ) : 
                for fm in fms :
                    for cm in cms :
                        command = ['python3.4', 'generate_report.py',
                                   '-if', path + filename,
                                   '-st', 'dna',
                                   '-rt', 'detailed',
                                   '-fm', fm,
                                   '-cm', cm,
                                   '-co', 'less',
                                   '-th', '1.00',
                                   '-sc', str(start_column)]

                        if subprocess.call(command) != 0 :
                            print(('Unable to analyze: {:s}'
                                   ''.format(path + filename)))


if __name__ == "__main__" :
    main()

