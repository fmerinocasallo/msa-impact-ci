#!/usr/bin/env python

import argparse
import re


nucleotides = ['-', 'a', 'c', 'g', 't']

def read_arguments () :
    """
    Read command line's arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Show columns where the '
                                                  'consensus residue differs '
                                                  'from each report file')
                                    )

    parser.add_argument('-f1', '--filename1', required=True,
                        help=('filename of the first report file'))

    parser.add_argument('-f2', '--filename2', required=True,
                        help=('filename of the second report file'))

    parser.add_argument('-df', '--min_diff', default=1.0, type=float,
                        help=('min difference between the consensus residues of '
                              'a fixed position at both report files to consider '
                              'it significant (PERCENTAGE)'))

    return parser.parse_args()


def get_max (line) :
    """
    Returns a tuple which contains the residue which appears the most in the
    given position as well as its appearance frequency

    Arguments :
        line ( string )
            Line of text from a detailed conservation report

    Returns :
        Tuple
            Tuple containing both the most frequent residue of the given line
            and its frequency. It looks like this: ( max_residue, max_freq)
    """
    regex_residue = '[' + "".join(nucleotides) + ']'
    regex_float = "\d+\.\d{4}"
    regex_stats = "\'" + regex_residue + "\': +" + regex_float

    stats = re.findall(regex_stats, line)
    # Look for the residue w/ the highest frequency
    max_residue = ''
    max_freq = 0.0
    for stat in stats :
        residue_freq = float(re.findall(regex_float, stat)[0])
        if residue_freq > max_freq :
            max_residue = re.findall(regex_residue, stat)[0]
            max_freq = residue_freq

    return max_residue, max_freq
    
def main () : 
    args = read_arguments()
    with open(args.filename1, 'rU') as f1, open(args.filename2, 'rU') as f2:
        # We are only interested in the last section of the report file
        while 'There' not in next(f1):
            pass

        next(f1)

        # We are only interested in the last section of the report file
        while 'There' not in next(f2):
            pass

        next(f2)

        print('line,type of diff,residue X,residue Y,freq X,freq Y')
        for l1, l2 in zip(f1, f2):
            max_r1, max_f1 = get_max(l1)
            max_r2, max_f2 = get_max(l2)

            num_line = re.findall('\d+', l1.split(':')[0])[0]
            if ((max_r1 == '-') and (max_r2 != '-') or
                (max_r1 != '-') and (max_r2 == '-')):
                print(('{:s},gap,{:s},{:s},{:.4f},{:.4f}'
                       ''.format(num_line, max_r1, max_r2, max_f1, max_f2)))
            elif max_r1 == max_r2 :
                if abs(max_f1 - max_f2) >= args.min_diff :
                    print(('{:s},freq,{:s},{:s},{:.4f},{:.4f}'
                           ''.format(num_line, max_r1, max_r2, max_f1, max_f2)))
            else :
                print(('{:s},res,{:s},{:s},{:.4f},{:.4f}'
                       ''.format(num_line, max_r1, max_r2, max_f1, max_f2)))


if __name__ == "__main__" :
    main()

