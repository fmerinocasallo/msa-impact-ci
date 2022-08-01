#!/usr/bin/env python

import argparse
import collections
import copy
import os.path

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
    parser = argparse.ArgumentParser(description=('Extract genes from the given '
                                                  'alignment')
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file (FASTA format) which '
                              'stores the aligned set of sequences'))

    return parser.parse_args()



def extract_gene ( alignment_filename, gene_section ) :
    """
    [ Generator function ]
    Extract given gene section from each of the sequences included on the
    given alignment.

    Arguments :
        alignment_filename ( string )
            File storing an aligned set of sequences (FASTA format).
        gene_section ( tuple )
            Tuple containing gene section's boundaries (lower and upper).
    """
    for record in Bio.SeqIO.parse(alignment_filename, 'fasta') :
        edited_record = copy.copy(record)
        edited_record.seq = record.seq[gene_section[0]:gene_section[1]]

        yield edited_record


def retrieve_new_genes_sections ( alignment_filename ) :
    """
    Retrieve new gene's sections based on rCRS's and the structure of the given
    alignment.

    Arguments :
        alignment_filename ( string )
            File storing an aligned set of sequences (FASTA format).

    Returns :
        dictionary
            Dictionary with gene identifiers as keys and its lower and upper
            boundaries as values (as tuples: (lower, upper))
    """
    genes = collections.OrderedDict([( 'tRNA-Phe', (  576,   646)),
                                     ( 'rRNA-12S', (  647,  1600)),
                                     ( 'tRNA-Val', ( 1601,  1669)),
                                     ( 'rRNA-16S', ( 1670,  3228)),
                                     ('tRNA-Leu1', ( 3229,  3303)),
                                     (      'ND1', ( 3306,  4261)),
                                     ( 'tRNA-Ile', ( 4262,  4330)),
                                     ( 'tRNA-Gln', ( 4328,  4399)),
                                     ( 'tRNA-Met', ( 4401,  4468)),
                                     (      'ND2', ( 4469,  5510)),
                                     ( 'tRNA-Trp', ( 5511,  5578)),
                                     ( 'tRNA-Ala', ( 5586,  5654)),
                                     ( 'tRNA-Asn', ( 5656,  5728)),
                                     ( 'tRNA-Cys', ( 5760,  5825)),
                                     ( 'tRNA-Tyr', ( 5825,  5890)),
                                     (     'COX1', ( 5903,  7444)),
                                     ('tRNA-Ser1', ( 7445,  7513)),
                                     ( 'tRNA-Asp', ( 7517,  7584)),
                                     (     'COX2', ( 7585,  8268)),
                                     ( 'tRNA-Lys', ( 8294,  8363)),
                                     (     'ATP8', ( 8365,  8571)),
                                     (     'ATP6', ( 8526,  9206)),
                                     (     'COX3', ( 9206,  9989)),
                                     ( 'tRNA-Gly', ( 9990, 10057)),
                                     (      'ND3', (10058, 10403)),
                                     ( 'tRNA-Arg', (10404, 10468)),
                                     (     'ND4L', (10469, 10765)),
                                     (      'ND4', (10759, 12137)),
                                     ( 'tRNA-His', (12137, 12205)),
                                     ('tRNA-Ser2', (12206, 12264)),
                                     ('tRNA-Leu2', (12265, 12335)),
                                     (      'ND5', (12336, 14147)),
                                     (      'ND6', (14148, 14672)),
                                     ( 'tRNA-Glu', (14673, 14741)),
                                     (     'CYTB', (14746, 15886)),
                                     ( 'tRNA-Thr', (15887, 15952)),
                                     ( 'tRNA-Pro', (15955, 16022))])
                                     # (   'D-loop': (16023,   575))])

    # Find rCRS sequence in the alignment
    found = False
    rCRS_id = 'NC_012920.1'
    rCRS = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''))
    for record in Bio.SeqIO.parse(alignment_filename, 'fasta') :
        if record.id == rCRS_id :
            rCRS = record
            found = True

            break

    if found :
        last_nc = -1 
        old_index = -1
        next_gene = True
        new_genes = collections.OrderedDict()
        genes_ids = (gene_id for gene_id in genes.keys())

        for new_index, bp in enumerate(rCRS) :
            if next_gene :
                try :
                    gene_id = next(genes_ids)
                except StopIteration :
                    break

                old_start, old_end = genes[gene_id]
                next_gene = False

            if bp != '-' :
                old_index += 1
                if old_index < old_start :
                    pass
                    # This could be the lower boundary of this gene
                    # last_nc = new_index
                elif old_index == old_start :
                    # This was the previous lower boundary of this gene
                    # new_start = last_nc + 1
                    new_start = new_index
                elif old_start < old_index <= old_end :
                    # We are traversing the previous gene's section
                    pass
                else : # old_end < old_index
                    # We've passed the upper boundary of this gene
                    new_end = new_index
                    new_genes[gene_id] = (new_start, new_end)
                    next_gene = True

            else : # bp == '-'
                pass

        return new_genes

    else : # not found
        print('WARNING!: Unable to find rCRS in the given alignment')

        exit(1)

def main () :
    args = read_arguments()

    genes_sections = retrieve_new_genes_sections(args.input_filename)

    for gene_id, gene_section in genes_sections.items() :
        filename, extension = os.path.splitext(args.input_filename)
        output_filename = filename + '_' + gene_id + extension
        Bio.SeqIO.write(extract_gene(args.input_filename, gene_section),
                        output_filename, 'fasta')

if __name__ == "__main__" :
    main()

