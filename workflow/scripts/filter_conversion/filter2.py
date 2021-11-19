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

from custom_tblat_filt import rev_complement, get_cytosines, get_reference_seq, mismatch

from pyfaidx import Fasta

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

			count = count_cs(combined_seq, ignore_motif)
			if within_cutoff(count, len(combined_seq), filter_type, cutoff):
				r1output.write(record1)
				r2output.write(record2)
			else: # Entry has a high number of Cs, either remove from output or add highc tag
				if remove_or_tag == 't':
					record1.name = 'highc_'+record1.name[0:]
					record2.name = 'highc_'+record2.name[0:]
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

	for i, read in enumerate(bam):
		read1 = read
		read2 = next(bam)
		remainder = 0
		if i % 2 == remainder: #read lines in pairs
			if read1.query_name == read2.query_name: # check if reads are paired
				# Change second read to reverse compliment
				seq1 = str(read1.query_sequence).encode()
				seq2 = str(Seq(read2.query_sequence).reverse_complement()).encode()
				qual1 = pysam.qualities_to_qualitystring(read1.query_qualities).encode()
				qual2 = pysam.qualities_to_qualitystring(read2.query_qualities)[::-1].encode()


				# Implement overlapping reads merger to get merged sequence:
				combined_seq = ctypes.create_string_buffer(len(read2.query_sequence)+len(read1.query_sequence))
				combined_seq = combine_reads(seq1, seq2, qual1, qual2, combined_seq, min_overlap, max_overlap, max_mismatch_density, allow_outies).decode()

				count = count_cs(combined_seq, ignore_motif)
				if within_cutoff(count, len(combined_seq), filter_type, cutoff):
					bam_out.write(read1)
					bam_out.write(read2)
				else: # If not within C cutoff, remove from output or tag with 'highc'
					if remove_or_tag == 't':
						read1.query_name = 'highc'+read1.query_name[1:]
						read2.query_name = 'highc'+read2.query_name[1:]
						bam_out.write(read1)
						bam_out.write(read2)
					cut+=1
			else: # If only one read is available, do the same but for 1/2 the cutoff
				if read1.is_read2: #change to reverse compliment
					seq = str(Seq(read1.query_sequence).reverse_complement())
				else:
					seq = read1.query_sequence
				count = count_cs(seq, ignore_motif)
				if within_cutoff(count, len(seq), filter_type, cutoff/2):
					bam_out.write(read1)
				else:
					if remove_or_tag == 't': #either remove or tag with 'highc'
						read1.query_name = '@highc'+read1.query_name[1:]
						bam_out.write(read1)
					cut+=1
				#Change remainder that governs the loop to take pairs
				if remainder == 0:
					remainder = 1
				else:
					remainder = 0
			num+=1
	print(f'Processed {num} paired-end reads')
	if remove_or_tag =='t':
		print(f'Tagged {cut}({(cut/num)*100:.2f}%) of reads')
	else:
		print(f'Removed {cut}({(cut/num)*100:.2f}%) of reads')
	bam_out.close()

	if os.path.exists('./temp.bam'):
		os.remove('temp.bam') # delete temporary sorted BAM file


def edit_custom_tblat(r1, r2, tblatpe, output, ref):
	"""
	tblat files typically contain metagenomic data where we can use
	species level inferences to remove low-cytosine reads. This is handled later
	as part of deeper analysis. This module here only creates a mismatch String
	.....T....C....T...... where:
	. = non C nucleotide
	T = "unmethylated" (ie converted) cytosine
	C = "methylated" (ie unconverted) cytosine
	"""
	fasta = Fasta(ref)
	with open(tblatpe) as f, open(output, 'w') as w, dnaio.open(r1) as fqr1, dnaio.open(r2) as fqr2:
		w.write('\t'.join(['taxid', 'qseqid', 'strand', 'sseqid', \
							'pident_R1', 'length_R1', 'mismatch_R1', 'gapopen_R1', 'qstart_R1', \
							'qend_R1', 'sstart_R1', 'send_R1', 'evalue_R1', 'bitscore_R1', 'qlen_R1', \
							'pident_R2', 'length_R2', 'mismatch_R2', 'gapopen_R2', 'qstart_R2', \
							'qend_R2', 'sstart_R2', 'send_R2', 'evalue_R2', 'bitscore_R2', 'qlen_R2', \
							'genome_len', 'effective_length', 'mismatch_string', 'molecule', 'reference_sequence'])+ '\n')
		for entry, record1, record2 in zip(f, fqr1, fqr2):
			w.write(mismatch(entry, record1, record2, fasta))

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


def count_cs(sequence, ignore_motif):
	if ignore_motif is None:
		return sequence.count('C')
	else:
		#number of upper case Cs = the number to count from the motif
		num_to_count = ignore_motif.count('C')
		ignore_motif = ignore_motif.upper() # change motif to upper case for counting
		count = sequence.count('C')
		motif_count = num_to_count * ( sequence.count(ignore_motif) ) # multiply the motif count by the number of Cs to count
		return(count-motif_count) # subtract out the motifs from the count





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-fn',"--filenames", help="r1 r2 paired or BAM or CUSTOM_TBLAT (r1, r2 & tblat file) input filenames", nargs='+')
	parser.add_argument('-ft',"--filter", type=str,help="type of filter used to detect contamination, either n (number) or p (percent)", choices = ['n','p'], nargs = 1)
	parser.add_argument('-c',"--cutoff", type=int,help="percentage of Cs or total # of Cs used for cutoff", nargs = 1)
	parser.add_argument('-rt','--remove_or_tag',type=str, help = 'Remove entry or tag with highc', nargs = 1, choices = ['r','t'])
	parser.add_argument('-d',"--output_dir", type=str, help = 'Path to output, optional', nargs = '?', default = os.getcwd())
	parser.add_argument('-t',"--output_tag", type=str, help = 'String added to beginning of edited filenames', nargs = '?', default = 'edited.')
	parser.add_argument('-m', "--ignore_motif", type=str, help = 'Motif including C of pattern to ignore within C count for filtering', nargs = '?', default = None)
	parser.add_argument("--output_tblat", type=str, help = 'output filename for filtered tblat', nargs = '?', default = None)
	parser.add_argument("--metagenomic_reference", type=str, help = 'metagenomic reference used (for tblat)', nargs = '?', default = None)

	args = parser.parse_args()
	rs = args.filenames

	# Import edited FLASH read overlapping code:
	combine = CDLL('/workdir/apc88/WGBS_pipeline/scripts/filter_conversion/combine.so')
	combine_reads = combine.combine_reads
	combine_reads.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_int, c_int, c_float, c_bool]
	combine_reads.restype = c_char_p #make sure the output is a char

	if'.bam' not in rs[0]:
		if (len(rs) == 2): #implies paired-end FQ
			r1i = rs[0]
			r2i = rs[1]
			start = time.time()
			editfastqfile(r1i, r2i, args.filter[0], args.cutoff[0], args.remove_or_tag[0], args.output_tag, args.output_dir, args.ignore_motif)
			end = time.time()
			print(f'Files {r1i} and {r2i} took {end-start:.3f} seconds to run')
		if (len(rs) == 3): #implies paired-end FQ with tblat file
			r1i = rs[0]
			r2i = rs[1]
			tblat = rs[2]
			start = time.time()
			edit_custom_tblat(r1i, r2i, tblat, args.output_tblat, args.metagenomic_reference)
			end = time.time()

			print(f'Files {r1i} and {r2i} took {end-start:.3f} seconds to run')

	else:
		bam_file = rs
		start = time.time()
		edit_bam(bam_file, args.filter[0], args.cutoff[0], args.remove_or_tag[0], args.output_tag, args.output_dir, args.ignore_motif)
		end = time.time()
		print(f'File {bam_file} took {end-start:.3f} seconds to run')
