#!/usr/bin/env python

from __future__ import print_function

import argparse
import copy
import os.path

import Bio.SeqIO

# maximum number of deletions considered valid in any given fragment
MAX_NEG_DIFF = 3
# maximum number of insertions considered valid in any given fragment
MAX_POS_DIFF = 2
def read_arguments () :
    """
    Read command line's arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Filter for each gene '
                                                  'fragments only those '
                                                  'with similar length to '
                                                  "rCRS's gene fragment")
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file (GenBank format) which '
                              'stores the set of sequences'))

    return parser.parse_args()


def get_gene_id ( feature, record_id ) :
    """
    Get gene identifier

    Arguments :
        feature ( SeqFeature )
            a SeqFeature from a SeqRecord object

    Returns :
        string
            a standard gene id
    """
    # rRNA feature's qualifiers are not consistent in GenBank
    rRNA_dict = {'s-rRNA': 'rRNA-12S',
                 'small subunit ribosomal RNA': 'rRNA-12S',
                 'l-rRNA': 'rRNA-16S',
                 'large subunit ribosomal RNA': 'rRNA-16S'}
    # CDS feature's qualifiers are not consistent in GenBank
    CDS_dict = {'ATPASE6': 'ATP6',
                'ATPASE 6': 'ATP6',
                'ATPase 6': 'ATP6',
                'ATPase subunit 6': 'ATP6',
                'ATP synthase 6': 'ATP6',
                'ATP synthase F0 subunit 6': 'ATP6',
                'ATPASE8': 'ATP6',
                'ATPASE 8': 'ATP6',
                'ATPase 8': 'ATP8',
                'ATPase subunit 8': 'ATP8',
                'ATP synthase 8': 'ATP8',
                'ATP synthase F0 subunit 8': 'ATP8',
                'CO1': 'COX1',
                'COI': 'COX1',
                'cytochrome c oxidase subunit 1': 'COX1',
                'cytochrome c oxidase subunit I': 'COX1',
                'cytochrome oxidase subunit I': 'COX1',
                'CO2': 'COX2',
                'COII': 'COX2',
                'cytochrome c oxidase subunit 2': 'COX2',
                'cytochrome c oxidase subunit II': 'COX2',
                'cytochrome oxidase subunit II': 'COX2',
                'CO3': 'COX3',
                'COIII': 'COX3',
                'cytochrome c oxidase subunit 3': 'COX3',
                'cytochrome c oxidase subunit III': 'COX3',
                'cytochrome oxidase subunit III': 'COX3',
                'NADH1': 'ND1',
                'NADH dehydrogenase subunit 1': 'ND1',
                'NADH2': 'ND2',
                'NADH dehydrogenase subunit 2': 'ND2',
                'NADH3': 'ND3',
                'NADH dehydrogenase subunit 3': 'ND3',
                'NADH4': 'ND4',
                'NADH dehydrogenase subunit 4': 'ND4',
                'NADH4L': 'ND4L',
                'NADH dehydrogenase subunit 4L': 'ND4L',
                'NADH5': 'ND5',
                'NADH dehydrogenase subunit 5': 'ND5',
                'NADH6': 'ND6',
                'NADH dehydrogenase subunit 6': 'ND6',
                'CYB' : 'CYTB',
                'cytb' : 'CYTB',
                'cytochrome b' : 'CYTB'}

    if ( feature.type in ['D-loop'] ) :
        gene_id = feature.type
    elif ( feature.type in ['tRNA'] ) :
        gene_id = feature.qualifiers['product'][0]
    elif ( feature.type in ['rRNA'] ) :
        feature_id = feature.qualifiers['product'][0]
        if ( feature_id in rRNA_dict ) :
            gene_id = rRNA_dict[feature_id]
        else : # feature_id not in rRNA_dict
            gene_id = 'rRNA-{}'.format(feature_id.split(' ')[0])
    elif ( feature.type in ['CDS'] ) :
        if ( 'gene' in feature.qualifiers ) :
            feature_id = feature.qualifiers['gene'][0].upper()
        else :
            feature_id = feature.qualifiers['product'][0]

        if ( feature_id in CDS_dict ) :
            gene_id = CDS_dict[feature_id]
        elif ( feature_id in CDS_dict.values() ) :
            gene_id = feature_id
        else:
            print(record_id)
            print(feature_id)
            gene_id = 'unkn'
    else : # feature.type in ['misc_feature', 'gene', 'STS']
        gene_id = 'unkn'

    return gene_id


def extract_gene_section ( input_filename, gene_id ) :
    """
    [ Generator function ]

    Extract 'gene_id' section for each sequence stored in 'input_filename'.
    We are going to store on disk only those gene fragments with a length
    similar to the rCRS's gene fragment. They will be later aligned.

    Arguments :
        input_filename ( string )
            File storing a set of sequences ( GenBank format).
        gene_id ( string )
            Gene section to be extracted  from the sequences.
    """
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

    gene_length = genes[gene_id][1] - genes[gene_id][0] + 1
    for record in Bio.SeqIO.parse(input_filename, 'gb') :
        if record.features :
            record_haplogroups = \
                record.features[0].qualifiers.get('haplogroup', None)
            if record_haplogroups :
                # find gene_id section inside record's features
                record_genes = []
                for feature in record.features[1:] :
                    gene = get_gene_id(feature, record.id)
                    if gene in ['tRNA-Ser', 'tRNA-Leu'] :
                        if (gene + '1') in record_genes :
                            gene +=  '2'
                        else :
                            gene += '1'

                        record_genes.append(gene)

                    if gene == gene_id :
                        gene_record = copy.copy(record)
                        gene_record.seq = feature.extract(record).seq
        
                        # Check if the number of indel mutations is valid
                        num_indels = len(gene_record.seq) - gene_length
                        if ( -MAX_NEG_DIFF <= num_indels <= MAX_POS_DIFF ) :
                            yield gene_record

                        break
            else :
                print('WARNING!: {:s},'.format(record.id))
        else :
            print('WARNING!: {:s},'.format(record.id))


def main () :
    genes_id = ['COX1', 'COX2', 'COX3']

    args = read_arguments()

    for gene_id in genes_id :
        filename, extension = os.path.splitext(args.input_filename)
        records_to_align_filename = filename + '_' + gene_id + extension

        Bio.SeqIO.write(extract_gene_section(args.input_filename, gene_id),
                        records_to_align_filename, 'gb')


if __name__ == "__main__" :
    main()

