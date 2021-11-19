#!/usr/bin/python3
"""
Does insilico conversion: converts all C-->T on R1
takes reverse complement of R2 and converts C-->T
"""

import sys
import os
import time
import itertools
import dnaio
import argparse
from Bio.Seq import Seq

def insilico_conversion(sequence, read):
	if read == 'R2': # if R2 read, converts to reverse complement and replaces C-->T
		return( str(Seq(sequence).reverse_complement()).replace('C','T') )
	else: # if R1, replaces Cs with Ts
		return(sequence.replace('C', 'T'))


def convert_file(infile, outfile, read):
	num=0
	start = time.time()
	with dnaio.open(infile) as f, dnaio.open(outfile, mode='w') as w:
		for entry in f:
			entry.sequence = insilico_conversion(entry.sequence, read)
			w.write(entry)
			num+=1
	end = time.time()
	print(f'Converted {num} reads in {end-start:.3f} seconds')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i',"--filename_in", type=str, help="R1 or R2 fasta/fastq input file", nargs=1)
	parser.add_argument('-o',"--filename_out", type=str,help="Name of fasta/fastq output file", nargs = 1)
	parser.add_argument('-r',"--read_type", type=str,help="Read type, R1 or R2", choices = ['R1','R2'], nargs = 1)
	args = parser.parse_args()
	convert_file(args.filename_in[0], args.filename_out[0], args.read_type[0])
	print(f'Completed insilico conversion on {args.filename_in[0]}, saved to {args.filename_out[0]}')
