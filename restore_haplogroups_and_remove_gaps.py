#!/usr/bin/env python

from __future__ import print_function
from six.moves import zip as _zip

import argparse
import copy
import os.path

import Bio.AlignIO
import Bio.Alphabet
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord

def read_arguments () :
    """
    Read command line's arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Restore haplogroup info and '
                                                  'remove columns with gaps in '
                                                  'the rCRS sequence')
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file (FASTA format) which '
                              'stores the aligned set of sequences'))
    parser.add_argument('-gi', '--gene_id', required=True,
                        help=('gene identifier'))

    return parser.parse_args()


def restore_haplogroups ( alignment_filename, haplogroups_filename ) :
    """
    [ Generator function ]
    Restore haplogroups information from a given set of sequences not aligned
    ( some of the alignment tools does not support GenBank format thus this
    information is lost when aligning sequences ).

    Arguments :
        alignment_filename ( string )
            File storing an aligned set of sequences (FASTA format).
        haplogroups_filename ( string )
            File storing same set of sequences but not aligned and with
            haplogroup info. 
    """
    alphabet = Bio.Alphabet.DNAAlphabet()
    for rec, rec_h in _zip(Bio.AlignIO.read(alignment_filename, 'fasta'),
                           Bio.SeqIO.parse(haplogroups_filename, 'gb')) :

        if ( rec.id == rec_h.id ) :
            record = copy.copy(rec_h)
            record.seq = rec.seq 
            record.seq.alphabet = alphabet

            yield record 
        else :
            print('WARNING!: {:s},{:s}'.format(rec.id, rec_h.id))

            return


def remove_gaps ( input_filename ) :
    """
    Remove alignment columns with a gap in every sequence and store on disk this
    new alignment.

    Arguments :
        input_filename ( string )
            Aligned set of sequences (GenBank format).
    """
    # Find rCRS sequence in the alignment
    found = False
    rCRS_id = 'NC_012920.1'
    rCRS = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''))
    for record in Bio.SeqIO.parse(input_filename, 'gb') :
        if record.id == rCRS_id :
            rCRS = record
            found = True

            break

    if not found :
        print('WARNING!: Unable to find rCRS in the given alignment')

        exit(1)

    # Traverse alignment removing columns where the rCRS sequence has gaps
    i = 0
    ranges = []
    len_rCRS = len(rCRS.seq) 
    while ( i < len_rCRS ) :
        start = i
        while ( i < len_rCRS ) and ( rCRS.seq[i] != '-' ) :
            i += 1
        end = i
        ranges.append((start, end))
        i += 1

    for record in Bio.SeqIO.parse(input_filename, 'gb') :
        edited_record = copy.copy(record)
        edited_record.seq = ""
        for start, end in ranges :
            edited_record.seq += record.seq[start:end]

        yield edited_record

    
def main () :
    args = read_arguments()

    # Restore haplogroup information
    split_filename = args.input_filename.split(args.gene_id)
    haplogroups_filename = (split_filename[0] + args.gene_id +
                            '_with_haplogroups.gb')
    filename, extension = os.path.splitext(haplogroups_filename)
    output_filename = filename + split_filename[1]
    filename, extension = os.path.splitext(output_filename)
    output_filename = filename + '.gb'
    Bio.SeqIO.write(restore_haplogroups(args.input_filename,
                                        haplogroups_filename),
                    output_filename, 'gb')

    # Remove alignment columns with gaps in rCRS sequence and store the updated
    # alignment on disk
    input_filename = output_filename
    filename, extension = os.path.splitext(output_filename)
    output_filename = filename + '_without_gaps' + extension
    Bio.SeqIO.write(remove_gaps(input_filename), output_filename, 'gb')


if __name__ == "__main__" :
    main()

