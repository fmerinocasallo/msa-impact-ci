#! /usr/bin/python3.3
#-------------------------------------------------------------------------------
# File :  length_analysis.py
# Description :  Script to analyze the length of a set of sequences
#
# Author :  F. Merino-Casallo  ( fmerino@unizar.es )
# Last version :  v0.1 ( 22/Aug/2014 )
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  02/Sep/2014
#   VERSION :  v0.1
#   AUTHOR(s) :  F. Merino-Casallo
#
#-------------------------------------------------------------------------------
import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser(description=('Length analysis'))
parser.add_argument('-sf', '--seqs_filename', required=True,
                    help='filename of the file storing the set of sequences')
parser.add_argument('-sfo', '--seqs_format', required=True,
                    help='format of the file storing the set of sequences')
args = parser.parse_args()

seqs = SeqIO.parse(args.seqs_filename, args.seqs_format)
lengths = {}
for seq in seqs:
	length = len(seq)
	if not length in lengths:
		lengths[length] = 1
	else:
		lengths[length] += 1

for length in sorted(lengths.keys()):
	print('Length: {:d}'.format(length))
	print('Reps: {:d}'.format(lengths[length]))

