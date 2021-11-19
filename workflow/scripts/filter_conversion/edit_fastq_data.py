

import pandas as pd
import numpy as np
import time
import dnaio
import ctypes
from ctypes import *
from Bio.Seq import Seq


### IMPORT DATA FROM FASTQ AND PUT INTO DATAFRAME 
start = time.time()

# Parameters for FLASH
max_overlap = 75 
min_overlap = 10
max_mismatch_density = 0.25
allow_outies = 1

# Import FLASH code:
combine = CDLL('./combine.so')
combine_reads = combine.combine_reads
combine_reads.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_int, c_int, c_float, c_bool]
combine_reads.restype = c_char_p  #make sure the output is a char

df = pd.DataFrame(columns = ['name','C_count','T_count','length','avg_prob_error'])

with dnaio.open("BK_MIX_R1.fastq.gz") as r1, dnaio.open("BK_MIX_R2.fastq.gz") as r2:
	i = 1
	for record1, record2 in zip(r1,r2):
		if i < 15000: # include only 30,000 records for dataset
			# Only take records which match the BK patients
			if record1.name.startswith('NB500947:293:HGNVLBGX2:1:11101') | record1.name.startswith('NS500503:259:HF52FBGXY:1:11101'):
				i+=1
				# add converted data to dataframe
				# convert sequences to correct format
				seq2 = str(Seq(record2.sequence).reverse_complement()).encode()
				qual2 = record2.qualities[::-1].encode()
				# Implement overlapping reads merger to get merged sequence:
				combined_seq = ctypes.create_string_buffer(len(record2.sequence)+len(record1.sequence))
				combined_seq = combine_reads(record1.sequence.encode(), seq2, record1.qualities.encode(), qual2, combined_seq, min_overlap, max_overlap, max_mismatch_density, allow_outies).decode()
				
				# Get count of Cs and Ts in the overlapped sequence
				c = combined_seq.count("C")
				t = combined_seq.count("T")
				length = len(combined_seq)
			
				# Create dictionary and append to dataframe
				d = {'name':record1.name, 'C_count':c, 'T_count':t,'length': length, 'avg_prob_error':np.mean([10**(-1*(ord(q)-33)/10) for q in record1.qualities])}
				df = df.append(d, ignore_index=True)
		else:
			break
		

# Add label for bisulfite vs standard sequenced reads
mask_bs = df['name'].str.startswith('NB500947:293:HGNVLBGX2:1:11101')
mask_std = df['name'].str.startswith('NS500503:259:HF52FBGXY:1:11101')
df.loc[mask_bs,'seqtype'] = 0 #the rows of reads that were bisulfite treated
df.loc[mask_std,'seqtype'] = 1 #rows of reads that were standard sequenced

print(df.head())

end = time.time()
print("Size of dataset: ",len(df))

df.to_csv("FASTQ_dataset.csv", index=False)
print("Time to create dataset: ",end-start)

