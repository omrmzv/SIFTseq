# Does insilico conversion ...
import sys
import os
import itertools
import dnaio
from Bio.Seq import Seq

def read_conversion(sequence, read):
	sequence = sequence.upper()
	if read == 'R2':
		endl = ''
		if (sequence[-1]=='\n'):
			endl = '\n'
			sequence = sequence.strip('\n')
		sequence = str(Seq(sequence).reverse_complement()) + endl
	sequence = sequence.replace('C', 'T')
	return(sequence)

def write_files(infile, outfile, read):
	with dnaio.open(infile) as f, dnaio.open(outfile, mode = 'w') as w:
		for record in f:
			record.sequence = read_conversion(record.sequence, read)
			w.write(record)


def parse_args():
	infile=sys.argv[1]
	outfile=sys.argv[2]
	read = sys.argv[3]
	return(infile, outfile, read)

if __name__ == '__main__':
	infile, outfile, read = parse_args()
	write_files(infile, outfile, read)
