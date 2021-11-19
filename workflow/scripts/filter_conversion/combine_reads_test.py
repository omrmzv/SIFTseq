
import os
import dnaio
import ctypes
from ctypes import *
from Bio.Seq import Seq
import time


#recompile the C code
os.system('gcc -std=c99 -I/usr/local/include/python2.7 -lpython2.7 -shared -Wl,-soname,combine -o combine.so -fPIC combine_reads_edited.c')

#load the shared object file
combine = CDLL('./combine.so')


# Initialize some variables
max_overlap = 75
min_overlap = 10
max_mismatch_density = 0.25
allow_outies = 1


combine_reads = combine.combine_reads
combine_reads.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_int, c_int, c_float, c_bool]
combine_reads.restype = c_char_p #make sure the output is a char


R1_file = 'R1short.fq.gz'
R2_file = 'R2short.fq.gz'

merged = open('merged_reads.txt', 'w')

start = time.time()
i = 0
comb=0
with dnaio.open(R1_file) as r1, dnaio.open(R2_file) as r2:
	for record1, record2 in zip(r1, r2):
		if i < 10:
			#Strings going into the C function must be bytes. Read 2 and qual 2 strings also must be reversed
			seq1 = record1.sequence.encode()
			qual1 = record1.qualities.encode()
			seq2 = str(Seq(record2.sequence).reverse_complement()).encode()
			qual2 = record2.qualities[::-1].encode()
			
			#allocate memory on the Python side for the returned string to go into
			combined_seq = ctypes.create_string_buffer(len(record2.sequence)+len(record1.sequence))
			
			combined_seq = combine_reads(seq1, seq2, qual1, qual2, combined_seq, min_overlap, max_overlap, max_mismatch_density, allow_outies).decode()
			if combined_seq != 'Not Combined':
				comb+=1
				print(record1.sequence)
				print(record2.sequence[::-1])
				print(combined_seq+'\n')
				merged.write(record1.name)
				merged.write('\n')
				merged.write(combined_seq)
				merged.write('\n')
				
				i+=1
		else:
			break
			
merged.close()
end = time.time()
print(f'{comb} out of {i} sequences combined')
print(f'Took {end-start} seconds to run')


