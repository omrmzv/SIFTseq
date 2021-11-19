
import dnaio
import time
import gzip
import argparse
import os
import pysam
import subprocess
from Bio.Seq import Seq
import ctypes
from ctypes import *

# Script to edit FASTQ, FASTA, and BAM files for the number of Cs in each sequence by first merging the overlapping reads, saving the edited file in a separate .gz compressed file called 'edited.original_filename.gz'
# Utilizes repos found here along with biopython: https://github.com/marcelm/dnaio
#							 					  https://github.com/betteridiot/bamnostic

# Inputs: 
# -R1,R2 file name pairs
# -filter: either 'n' for number of Cs or 'p' for percentage of Cs
# -cutoff: for n, the number of Cs and for p, a number between 0 and 1 of the percentage for cutoff
# -output_dir: directory to save edited files to, optional
# -output_tag: string added to beginning of edited filename, optional (default='edited')
# -ignore_motif: motif containing Cs to ignore the count of (e.g. CG, CHH, CHG), user can define their own motif where the counted C is in capitals 
# 		+CCGG = ignore both cytosines in this motif
# 		+cCGG = only ignore the second cytosine
# 		+CcGG = only ignore the first cytosine
# -remove_or_tag: either 'r' or 't', indicating whether to remove reads that do not meet the threshold or to simply tag them with 'highc'


def editfastqfile(r1_in, r2_in, filter_type, cutoff, remove_or_tag, output_tag, output_dir, ignore_motif):
	r1_out = output_dir+'/'+output_tag+r1_in
	r2_out = output_dir+'/'+output_tag+r2_in
	r1output = dnaio.open(r1_out, mode='w')
	r2output = dnaio.open(r2_out, mode='w')

	max_overlap = 75 #parameters for overlapping read merge function
	min_overlap = 10
	max_mismatch_density = 0.25
	allow_outies = 1
	
	with dnaio.open(r1_in) as r1, dnaio.open(r2_in) as r2:
		ent = 0
		cut = 0
		for record1, record2 in zip(r1,r2):
			# R2 is the reverse complement- take reverse to edit for Cs
			seq2 = str(Seq(record2.sequence).reverse_complement()).encode()
			qual2 = record2.qualities[::-1].encode()
			# Implement overlapping reads merger to get merged sequence:
			combined_seq = ctypes.create_string_buffer(len(record2.sequence)+len(record1.sequence))
			combined_seq = combine_reads(record1.sequence.encode(), seq2, record1.qualities.encode(), qual2, combined_seq, min_overlap, max_overlap, max_mismatch_density, allow_outies).decode()
			
			count = count_cs_fastaq(combined_seq, ignore_motif)
			if within_cutoff(count, len(combined_seq), filter_type, cutoff):
				r1output.write(record1)
				r2output.write(record2)
			else: # Entry has a high number of Cs, either remove from output or add highc tag
				if remove_or_tag == 't':
					record1.name = 'highc'+record1.name[1:]
					record2.name = 'highc'+record2.name[1:]
					r1output.write(record1)
					r2output.write(record2)
				cut+=1
			ent += 1
	print(f'Processed {ent} paired-end reads')
	if remove_or_tag =='t':
		print(f'Tagged {cut}({(cut/ent)*100:.2f}%) of reads')
	else:
		print(f'Removed {cut}({(cut/ent)*100:.2f}%) of reads')
	r1output.close()
	r2output.close()
				 
				 
def edit_bam(bam_file, filter_type, cutoff, remove_or_tag, output_tag, output_dir, ignore_motif):
	num = 0 #total number of paired-end reads
	cut = 0 #number of paired-end reads cut from file
	
	max_overlap = 75 #parameters for overlapping read merge function
	min_overlap = 10
	max_mismatch_density = 0.25
	allow_outies = 1
	
	#create output file
	bam = pysam.AlignmentFile(bam_file, 'rb') #read in bam file
	bam_out = pysam.AlignmentFile(output_dir+'/'+output_tag+bam_file, 'wb', template = bam) 
	if bam.header['HD']['SO'] == 'unsorted':
		# if the file is unsorted, sort by name and save to temporary file
		pysam.sort('-o', 'temp.bam', '-n', bam_file)
		bam = pysam.AlignmentFile('temp.bam', 'rb')
		
	# Check to see if paired-end or single sequencing:
	# if first two reads are paired, assume the whole file is
	if next(bam).is_paired and next(bam).is_paired:
		#treatment for paired reads 
		for i,read in enumerate(bam):
			read1 = read
			read2 = next(bam)
			if i % 2 == 0: #reading pairs of reads together
				#check for overlap based on the alignment positions  
				xm1 = read1.get_tags()[2][1]
				xm2 = read2.get_tags()[2][1]
				seq1 = read1.query_sequence
				seq2 = str(Seq(read2.query_sequence).reverse_complement())
				if read1.reference_end > read2.reference_start: #reads are overlapping
					olap = read1.reference_end - read2.reference_start
					#remove the end of sequence and xm string for read 1 to prevent double-counting
					seq1 = seq1[0:len(seq1)-olap-1]
					xm1 = xm1[0:len(xm1)-olap-1]
					count = count_cs_bam(xm1+xm2, seq1+seq2, ignore_motif)
				else: #not overlapping, treat the two strings differently
					count = count_cs_bam(xm1+xm2, seq1+seq2, ignore_motif)
					
				#check if it is within the cutoff:
				if within_cutoff(count, len(seq1+seq2), filter_type, cutoff):
					bam_out.write(read1)
					bam_out.write(read2)
				else: # If not within C cutoff, remove from output or tag with 'highc'
					if remove_or_tag == 't': 
						read1.query_name = 'highc'+read1.query_name[1:]
						read2.query_name = 'highc'+read2.query_name[1:]
						bam_out.write(read1)
						bam_out.write(read2)
					cut+=1		
			num+=1
		
	else: #treatment for single reads 
		for read in bam:
			xm = read.get_tags()[2][1]
			count = count_cs_bam(xm, read.query_sequence, ignore_motif)
					
			#check if it is within the cutoff:
			if within_cutoff(count, len(xm), filter_type, cutoff):
				bam_out.write(read)
			else: # If not within C cutoff, remove from output or tag with 'highc'
				if remove_or_tag == 't': 
					read.query_name = 'highc'+read.query_name[1:]
					bam_out.write(read)
				cut+=1		
			num+=1
			
	print(f'Processed {num} paired-end reads')
	if remove_or_tag =='t':
		print(f'Tagged {cut}({(cut/num)*100:.2f}%) of reads')
	else:
		print(f'Removed {cut}({(cut/num)*100:.2f}%) of reads')
	bam_out.close()
		
	if os.path.exists('./temp.bam'):
		os.remove('temp.bam') # delete temporary sorted BAM file 
			
				 
def within_cutoff(count, sequence_len, filter_type, cutoff):
	if filter_type == 'n': # if editing the file by total number of Cs:
		if count < cutoff:
			return True
		else:
			return False
	else: # if editing the file by percentage of Cs:
		if ( count / sequence_len ) < cutoff:
			return True
		else:
			return False

def count_cs_fastaq(sequence, ignore_motif):
	if ignore_motif is None:
		return sequence.count('C')
	else:
		#number of upper case Cs = the number to count from the motif
		num_to_count = ignore_motif.count('C') 
		ignore_motif = ignore_motif.upper() # change motif to upper case for counting
		count = sequence.count('C')
		motif_count = num_to_count * ( sequence.count(ignore_motif) ) # multiply the motif count by the number of Cs to count
		return(count-motif_count) # subtract out the motifs from the count
	 
	 
def count_cs_bam(xm_str, sequence, ignore_motif):
	'''
	xm_str: the full XM string for both sequences put together
	sequence: the full sequence for both reads put together
	ignore_motif: None of the string specifying which motifs to ignore qhen counting Cs 
	Counts cs using the full (not overlapped) sequence and XM strings based on the ignore_motif'''
	if ignore_motif is None: #count Cs based on XM flag
		return xm_str.count('c')+xm_str.count('C')+xm_str.count('x')+xm_str.count('X')+xm_str.count('h')+xm_str.count('H')
	else:					#for x/X tag		
		if ignore_motif == 'CG': #for c/C tag
			ignore_cs = xm_str.count('c')+xm_str.count('C')
			
		else: #otherwise, use the sequence to count cs and motifs 
			#number of upper case Cs = the number to count from the motif
			num_to_count = ignore_motif.count('C') 
			ignore_motif = ignore_motif.upper() # change motif to upper case for counting
			# multiply the motif count by the number of Cs to count
			ignore_cs = num_to_count * ( sequence.count(ignore_motif) ) 
		
		return(sequence.count('C')-ignore_cs) # subtract out the motifs from the count
	 
			


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-fn',"--filenames", help="r1 r2 paired or BAM input filenames", nargs='+')
	parser.add_argument('-ft',"--filter", type=str,help="type of filter used to detect contamination, either n (number) or p (percent)", choices = ['n','p'], nargs = 1)
	parser.add_argument('-c',"--cutoff", type=int,help="percentage of Cs or total # of Cs used for cutoff", nargs = 1)
	parser.add_argument('-rt','--remove_or_tag',type=str, help = 'Remove entry or tag with highc', nargs = 1, choices = ['r','t'])
	parser.add_argument('-d',"--output_dir", type=str, help = 'Path to output, optional', nargs = '?', default = os.getcwd())
	parser.add_argument('-t',"--output_tag", type=str, help = 'String added to beginning of edited filenames', nargs = '?', default = 'edited.')
	parser.add_argument('-m', "--ignore_motif", type=str, help = 'Motif including C of pattern to ignore within C count for filtering', nargs = '?', default = None)

	args = parser.parse_args()
	rs = args.filenames
	
	# Import edited FLASH read overlapping code:
	combine = CDLL('./combine.so')
	combine_reads = combine.combine_reads
	combine_reads.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_int, c_int, c_float, c_bool]
	combine_reads.restype = c_char_p #make sure the output is a char
	
	if'.bam' not in rs[0]: 
		if len(rs) % 2 != 0:
			raise ValueError('Number of R1 and R2 inputs do not match.')
		r1in = rs[0:len(rs)-1:2]
		r2in = rs[1:len(rs):2]
		
		for r1i, r2i in zip(r1in, r2in):
			start = time.time()
			editfastqfile(r1i, r2i, args.filter[0], args.cutoff[0], args.remove_or_tag[0], args.output_tag, args.output_dir, args.ignore_motif)
			end = time.time()
			print(f'Files {r1i} and {r2i} took {end-start:.3f} seconds to run')
		
	else:
		for bam_file in rs: 
			start = time.time()
			edit_bam(bam_file, args.filter[0], args.cutoff[0], args.remove_or_tag[0], args.output_tag, args.output_dir, args.ignore_motif)
			end = time.time()
			print(f'File {bam_file} took {end-start:.3f} seconds to run')


