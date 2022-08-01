#i!/usr/bin/env python

from __future__ import print_function

import argparse
import copy
import os.path
import tempfile

import Bio.SeqIO
import Bio.Seq

import MEvoLib.Align

def read_arguments () :
    """
    Read command line's arguments.

    Returns :
        ArgumentParser
            Arguments passed by the user through the shell.
    """
    parser = argparse.ArgumentParser(description=('Extract gene sections for '
                                                  'each given sequence')
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file (GenBank format) which '
                              'stores the set of sequences'))

    parser.add_argument('-gi', '--gene_id', required=True,
                        help='gene identifier to be remove from sequences')

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
                 '12S rRNA': 'rRNA-12S',
                 'l-rRNA': 'rRNA-16S',
                 'large subunit ribosomal RNA': 'rRNA-16S',
                 '16S rRNA': 'rRNA-16S'}
    # CDS feature's qualifiers are not consistent in GenBank
    CDS_dict = {'MT-ATP6': 'ATP6',
                'ATPASE6': 'ATP6',
                'ATPase6': 'ATP6',
                'ATPASE 6': 'ATP6',
                'ATPase 6': 'ATP6',
                'ATPASE 6 GENE': 'ATP6',
                'ATPase subunit 6': 'ATP6',
                'ATP synthase 6': 'ATP6',
                'ATP synthase subunit 6': 'ATP6',
                'ATP synthase F0 subunit 6': 'ATP6',
                'F0-ATP synthase subunit 6': 'ATP6',
                'MT-ATP8': 'ATP8',
                'ATPASE8': 'ATP8',
                'ATPase8': 'ATP8',
                'ATPASE 8': 'ATP8',
                'ATPase 8': 'ATP8',
                'ATPASE 8 GENE': 'ATP8',
                'ATPase subunit 8': 'ATP8',
                'ATP synthase 8': 'ATP8',
                'ATP synthase subunit 8': 'ATP8',
                'ATP synthase F0 subunit 8': 'ATP8',
                'F0-ATP synthase subunit 8': 'ATP8',
                'MT-CO1': 'COX1',
                'CO1': 'COX1',
                'COI': 'COX1',
                'CO I': 'COX1',
                'COXI': 'COX1',
                'COX 1': 'COX1',
                'COX I': 'COX1',
                'cytochrome oxidase I': 'COX1',
                'cytochrome oxidase subunit 1': 'COX1',
                'cytochrome-c oxidase I': 'COX1',
                'cytochrome c oxidase subunit 1': 'COX1',
                'cytochrome c oxidase subunit I': 'COX1',
                'cytochrome-c oxidase subunit I': 'COX1',
                'cytochrome oxidase subunit I': 'COX1',
                'CYTOCHROME OXIDASE SUBUNIT I': 'COX1',
                'MT-CO2': 'COX2',
                'CO2': 'COX2',
                'COII': 'COX2',
                'CO II': 'COX2',
                'cCOII': 'COX2',
                'CCOII': 'COX2',
                'COXII': 'COX2',
                'COX 2': 'COX2',
                'COX II': 'COX2',
                'cytochrome oxidase II': 'COX2',
                'cytochrome oxidase subunit 2': 'COX2',
                'cytochrome-c oxidase II': 'COX2',
                'cytochrome c oxidase subunit 2': 'COX2',
                'cytochrome c oxidase subunit II': 'COX2',
                'cytochrome-c oxidase subunit II': 'COX2',
                'cytochrome oxidase subunit II': 'COX2',
                'CYTOCHROME OXIDASE SUBUNIT II': 'COX2',
                'MT-CO3': 'COX3',
                'CO3': 'COX3',
                'COIII': 'COX3',
                'CO III': 'COX3',
                'COXIII': 'COX3',
                'COX III': 'COX3',
                'COX 3': 'COX3',
                'cytochrome oxidase III': 'COX3',
                'cytochrome oxidase subunit 3': 'COX3',
                'cytochrome-c oxidase III': 'COX3',
                'cytochrome c oxidase subunit 3': 'COX3',
                'cytochrome c oxidase subunit III': 'COX3',
                'cytochrome-c oxidase subunit III': 'COX3',
                'cytochrome oxidase subunit III': 'COX3',
                'CYTOCHROME OXIDASE SUBUNIT III': 'COX3',
                'MT-ND1': 'ND1',
                'NAD1': 'ND1',
                'NADH1': 'ND1',
                'NADH subunit 1': 'ND1',
                'NADH dehydrogenase subunit 1': 'ND1',
                'MT-ND2': 'ND2',
                'NAD2': 'ND2',
                'NADH2': 'ND2',
                'NADH subunit 2': 'ND2',
                'NADH dehydrogenase subunit 2': 'ND2',
                'MT-ND3': 'ND3',
                'NAD3': 'ND3',
                'NADH3': 'ND3',
                'NADH subunit 3': 'ND3',
                'NADH dehydrogenase subunit 3': 'ND3',
                'MT-ND4': 'ND4',
                'NAD4': 'ND4',
                'NADH4': 'ND4',
                'NADH subunit 4': 'ND4',
                'NADH dehydrogenase subunit 4': 'ND4',
                'MT-ND4L': 'ND4L',
                'NAD4L': 'ND4L',
                'NADH4L': 'ND4L',
                'NADH subunit 4L': 'ND4L',
                'NADH dehydrogenase subunit 4L': 'ND4L',
                'MT-ND5': 'ND5',
                'NAD5': 'ND5',
                'NADH5': 'ND5',
                'NADH subunit 5': 'ND5',
                'NADH dehydrogenase subunit 5': 'ND5',
                'MT-ND6': 'ND6',
                'NAD6': 'ND6',
                'NADH6': 'ND6',
                'NADH subunit 6': 'ND6',
                'NADH dehydrogenase subunit 6': 'ND6',
                'MT-CYTB': 'CYTB',
                'CYB' : 'CYTB',
                'CTYB' : 'CYTB',
                'cob' : 'CYTB',
                'COB' : 'CYTB',
                'cytb' : 'CYTB',
                'cytochrome b' : 'CYTB'}

    if ( feature.type in ['D-loop'] ) :
        gene_id = feature.type
    elif ( feature.type in ['tRNA'] ) :
        if ( 'product' in feature.qualifiers ) :
            gene_id = feature.qualifiers['product'][0]
        elif ( 'gene' in feature.qualifiers ) :
            gene_id = feature.qualifiers['gene'][0]
        else :
            print(record_id)
            gene_id = 'unkn'
    elif ( feature.type in ['rRNA'] ) :
        if ( 'product' in feature.qualifiers ) :
            feature_id = feature.qualifiers['product'][0]
        else:
            feature_id = feature.qualifiers['gene'][0]

        if ( feature_id in rRNA_dict ) :
            gene_id = rRNA_dict[feature_id]
        else : # feature_id not in rRNA_dict
            gene_id = 'rRNA-{}'.format(feature_id.split(' ')[0])
    elif ( feature.type in ['CDS'] ) :
        if ( ( 'gene' not in feature.qualifiers ) and
             ( 'product' not in feature.qualifiers ) ) :

            print(record_id)
            gene_id = 'unkn'
        else :
            if ( 'gene' in feature.qualifiers ) :
                feature_id = feature.qualifiers['gene'][0].upper()
            else : # ( 'product' in feature.qualifiers )
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

def extract_gene ( input_filename, gene_id ) :
    """
    [ Generator function ]

    Try to extract 'gene_id' sections for each sequence stored in 'input_filename'.
    We return a generator of tuples, each tuple would looks like this:
    - If the sequence match the expected structure: (SeqRecord, None)
    - If the sequence doesn't match the expected structure: (None, SeqRecord)

    Arguments :
        input_filename ( string )
            File storing a set of sequences ( GenBank format).
        gene_id ( string )
            Gene section to be removed from the sequences.

    Returns :
        generator
            SeqRecord tuples, following the scheme mention above.
    """
    def gene_search_and_extraction ( record, gene_id ) :
        """
        Look for 'gene_id' sections in 'record'

        Arguments :
            record ( SeqRecord )
                Record where we are going to look for 'gene_id' sections.
            gene_id ( string )
                Gene section to be extracted from the sequences.

        Returns :
            SeqRecord
                The very same 'record' with only 'gene_id' sections.
                If 'gene_id' is not present in 'record', the returned value
                will be None.
        """
        invalid_info = False

        if record.features :
            gene_found = False

            # find gene_id section inside record's features
            record_genes = []
            for feature in record.features[1:] :
                gene = get_gene_id(feature, record.id)
                # Both 'tRNA-Ser' and 'tRNA-Leu' are split into two sections.
                if gene in ['tRNA-Ser', 'tRNA-Leu'] :
                    if (gene + '1') in record_genes :
                        gene +=  '2'
                    else :
                        gene += '1'

                    record_genes.append(gene)

                if gene == gene_id :
                    gene_found = True
                    new_record = copy.copy(record)
                    first_section = True
                    # Extract 'gene_id' sections from record
                    if feature.location :
                        # For cases with multiple locations such as D-loop gene,
                        # locations may look like this: (16027, 16570) and
                        # (0, 577). We need locations to be sorted in order to
                        # extract those sections.
                        sorted_locations = \
                                    sorted(feature.location.parts,
                                           key=lambda location: location.start)

                        # We expect every sequence with D-loop to have two
                        # sections, one at the beginning of the sequence and
                        # another one at the end. If that's not the case, we
                        # will have to extract them using rCRS scheme.
                        if ( (gene_id == 'D-loop') and
                             (len(feature.location.parts) != 2) ) :
                            print('{:s},{:s},,,,,'.format(gene_id, record.id))
                            invalid_info = True
                        else : # ( (gene_id != 'D-loop') or
                               #   (len(feature.location.parts) == 2) )
                            for location in sorted_locations :
                                if first_section :
                                    if gene_id == 'D-loop' :
                                        # The first D-loop section should always
                                        # end around position 576.
                                        if location.end < 500 :
                                            print(('{:s},,{:s},,,,'
                                                   ''.format(gene_id, record.id)))
                                            invalid_info = True
                                        else : # location.end > 500
                                            # The first D-loop section should
                                            # always start at the beginning of
                                            # the sequence.
                                            if location.start > 0 :
                                                print(('{:s},,{:s},,,,'
                                                       ''.format(gene_id,
                                                                 record.id)))
                                                invalid_info = True
                                                # new_record.seq = \
                                                #    record.seq[0:location.end]
                                            else : # location.start = 0
                                                new_record.seq = \
                                                    record.seq[0:location.end]
                                    else : # gene_id != 'D-loop'
                                        new_record.seq = \
                                            record.seq[location.start:location.end]

                                    first_section = False
                                else : # first_section == False
                                    if gene_id == 'D-loop' :
                                        # The second D-loop section should always
                                        # end at the end of the sequence.
                                        if location.end < len(record.seq) - 1:
                                            print(('{:s},,{:s},,,,'
                                                   ''.format(gene_id, record.id)))
                                            invalid_info = True

                                        else : # location.end = len(record.seq) - 1
                                            pass

                                        new_record.seq += \
                                            record.seq[location.start:]
                                    else : # gene_id != 'D-loop'
                                        new_record.seq += \
                                            record.seq[location.start:location.end]

                                # If we have found any invalid information
                                # we have to stop the process and mark the
                                # sequence as unproccesable
                                if invalid_info :
                                    break
                                else : # invalid_info == False
                                    pass

                        # If we have found any invalid information
                        # we have to stop the process and mark the
                        # sequence as unproccesable
                        if invalid_info :
                            break
                        else : # invalid_info == False
                            return new_record

                    else : # feature.location == None
                        # There are no location info for this feature
                        print('{:s},,,{:s},,,'.format(gene_id, record.id))
                        invalid_info = True

                else : # gene != gene_id
                    pass

            if not gene_found :
                # There are no 'gene_id' sections in this record
                print('!{:s},,,,{:s},,'.format(gene_id, record.id))
                invalid_info = True

        else : # record.features == None
            # There are no features associated with this record
            print('{:s},,,,,{:s},'.format(gene_id, record.id))
            invalid_info = True

        if invalid_info :
            return None
        else :
            print('UNKNOWN: {:s}'.format(record.id))

    print(('gene id, D-loop w/ 2 sections,D-loop wrong boundaries,'
           'no locations info, no gene id features, no features, unable to '
           'extract'))

    for record in Bio.SeqIO.parse(input_filename, 'gb') :
        new_record = gene_search_and_extraction(record, gene_id)
        if new_record :
            # 'gene_id' sections have been successfully extracted
            yield new_record, None
        else :
            # Unable to extract 'gene_id' sections because 'record' doesn't
            # match the expected scheme
            yield None, record


def main () :
    args = read_arguments()

    filename, extension = os.path.splitext(args.input_filename)

    # Look for processable and unprocessable sequences and store them
    filename_processables = (filename + '_' + args.gene_id + '_' +
                             'processables.fasta')
    filename_unprocessables = (filename + '_' + args.gene_id + '_' +
                       'unprocessables.fasta')

    with open(filename_processables, 'w') as handle_processables, \
         open(filename_unprocessables, 'w') as handle_unprocessables :

        num_processables = 0
        num_unprocessables = 0

        sequences = extract_gene(args.input_filename, args.gene_id)
        for sequence in sequences :
            if sequence[0] != None :
                num_processables += Bio.SeqIO.write(sequence[0],
                                                    handle_processables,
                                                    'fasta')
            else :
                num_unprocessables += Bio.SeqIO.write(sequence[1],
                                                      handle_unprocessables,
                                                      'fasta')


        print(('We\'ve extracted {:s} gene fragments from {:d} sequences '
               'and store them in {:s}'
               ''.format(args.gene_id, num_processables,
                         filename_processables)))

        print(('We\'ve found {:d} sequences where we\'ve been unable to extract '
               '{:s} gene fragments and store them in {:s}'
               ''.format(num_unprocessables, args.gene_id,
                         filename_unprocessables)))


if __name__ == "__main__" :
    main()
