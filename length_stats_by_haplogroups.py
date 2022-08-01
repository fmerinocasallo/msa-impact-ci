#!/usr/bin/env python

from __future__ import print_function

import argparse
import collections
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


def get_gene_info_by_haplogroups ( input_filename, gene_id, haplogroups ) :
    """
    Get 'gene_id' fragments info for each sequence stored in 'input_filename'.
    We are going to generate some stats about it by haplogroups.

    Arguments :
        input_filename ( string )
            File storing a set of sequences ( GenBank format).
        gene_id ( string )
            Gene section to be extracted  from the sequences.
        haplogroups ( dictionary )
            Dictionary with haplogroups as keys and tuples as values. Each tuple
            as ( record_id, gene_fragment_length ).
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
                        haplogroups[record_haplogroup].\
                                setdefault(len(feature), []).append(record.id)

                        break
            else :
                print('WARNING!: {:s},'.format(record.id))
        else :
            print('WARNING!: {:s},'.format(record.id))


def print_stats ( gene_id, haplogroups ) :
    """
    Print previously generated stats about gene fragments lengths by
    haplogroups.

    Arguments :
        gene_id ( string )
            Gene section to be extracted  from the sequences.
        haplogroups ( dictionary )
            Dictionary with haplogroups as keys and a dictionary containing
            lengths as keys and a list of record ids as values.
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

    print('gene,haplogroup,length,diff,#records,%,record_ids,[...]')

    gene_length = genes[gene_id][1] - genes[gene_id][0] + 1
    for haplogroup, lengths in haplogroups.items() :
        # We're not interested in mt-MRCA haplogroup
        if ( haplogroup == 'mt-MRCA' ) or ( not lengths ) :
            break

        sum_records = sum(len(records_ids) for records_ids in lengths.values())
        for length, record_ids in sorted(lengths.items()) :
            num_records = len(record_ids)
            percentage = num_records / float(sum_records) * 100
            print(('{:s},{:s},{:d},{:+d},{:d},{:5.2f},{:s}'
                   ''.format(gene_id, haplogroup, length, length - gene_length,
                             num_records, percentage, ','.join(record_ids))))


def main () :
    genes_id = ['COX1', 'COX2', 'COX3']

    args = read_arguments()

    for gene_id in genes_id :
        # BEWARE! We are using this OrderedDict for regex matching and the order
        # of its items  matters. You should always include first those
        # haplogroups with longer name. That way, it will match them correctly.
        # For example, if we have both 'R' and 'R0' haplogroups, 'R0' must be
        # added first to the OrderedDict.
        haplogroups = collections.OrderedDict([('L0', {}), ('L1', {}),
                                               ('L2', {}), ('L3', {}),
                                               ('L4', {}), ('L5', {}),
                                               ('L6', {}),  ('M', {}),
                                                ('C', {}),  ('E', {}),
                                                ('G', {}),  ('Q', {}),
                                                ('Z', {}),  ('D', {}),
                                                ('N', {}),  ('A', {}),
                                                ('I', {}),  ('O', {}),
                                                ('S', {}),  ('W', {}),
                                                ('X', {}),  ('Y', {}),
                                               ('R0', {}), ('B6', {}),
                                                ('F', {}),  ('J', {}),
                                                ('P', {}),  ('K', {}),
                                                ('B', {}), ('HV', {}),
                                                ('T', {}),  ('R', {}),
                                                ('H', {}),  ('V', {}),
                                                ('U', {}),
                                                ('&', {}), ('mt-MRCA', {})])

        get_gene_info_by_haplogroups(args.input_filename, gene_id, haplogroups)

        print_stats(gene_id, haplogroups)


if __name__ == "__main__" :
    main()

