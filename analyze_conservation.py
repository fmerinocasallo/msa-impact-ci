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
    fms = ['unweighted', 'weighted']
    cms = ['entropy', 'variance']
    cos = ['greater', 'less']
    # ths = ['0.5', '0.75', '0.90', '0.95', '0.99', '1.00']
    ths = ['1.00']

    # path = 'data/dbs/'
    path = '2017-09-14/'
    args = read_arguments()

    # regex = 'hmtDNA_2015_11_09_*aligned_with_*.fasta'
    regex = 'mmtDNA_2017-09-14_*aligned_with_*.fasta'
    for filename in os.listdir(path) :
        if ( fnmatch.fnmatch(filename, regex) ) : 
            for fm in fms :
                for cm in cms :
                    for co in cos :
                        for th in ths :
                            command = ['python3.4',
                                       'generate_report.py',
                                       '-if', path + filename,
                                       '-st', 'dna',
                                       '-rt', 'detailed',
                                       # '-rt', 'basic',
                                       '-fm', fm,
                                       '-cm', cm,
                                       '-co', co,
                                       '-th', th]

                            print(('Analyzing {:s} using {:s} and {:s}... '
                                   ''.format(filename, fm, cm)), end="")
                            if subprocess.call(command) != 0 :
                                print(('Unable to analyze: {:s}'
                                       ''.format(path + filename)))
                            else:
                                print('done')


if __name__ == "__main__" :
    main()

