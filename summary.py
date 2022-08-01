#!/usr/bin/env python

from __future__ import print_function

import argparse

import Bio.SeqIO
import numpy


def read_arguments () :
    """
    Read command line's arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Retrieve a statistical '
                                                  "summary about sequences' "
                                                  'length')
                                    )
    parser.add_argument('-fi', '--filename', required=True,
                        help=('filename of the file which stores a set of '
                              'sequences'))

    parser.add_argument('-fo', '--format', default='gb',
                        help='sequences file format')

    return parser.parse_args()


def main () :
    args = read_arguments()

    lengths = [len(rec) for rec in Bio.SeqIO.parse(args.filename, args.format)]

    print('Summary of {:d} seqs'.format(len(lengths)))
    print('-----------------------')
    print('min    : {:f}'.format(numpy.amin(lengths)))
    print('1st Qu : {:f}'.format(numpy.percentile(lengths, 25)))
    print('median : {:f}'.format(numpy.median(lengths)))
    print('mean   : {:f}'.format(numpy.average(lengths)))
    print('stdev  : {:f}'.format(numpy.std(lengths)))
    print('3rd Qu : {:f}'.format(numpy.percentile(lengths, 75)))
    print('max    : {:f}'.format(numpy.amax(lengths)))

if __name__ == "__main__" :
    main()

