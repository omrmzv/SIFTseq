#!/usr/bin/env python3

import dnaio
import time
import gzip
import argparse
import os
import pysam
import subprocess
import joblib
import random
import textwrap
import ctypes
import pandas as pd
import numpy as np
from pyfaidx import Fasta
from Bio.Seq import Seq
from ctypes import *
import time

from joblib import Parallel, delayed
import multiprocessing
from functools import partial
from Bio import pairwise2

# Script to edit and/or provide C and T counts for FASTQ, FASTA, BAM, and TBLAT files.
# Utilizes repos found here along with biopython: https://github.com/marcelm/dnaio
#												  https://github.com/betteridiot/bamnostic

# Inputs:
# -file name (one name for BAM and TBLAT, two for FASTA or FASTQ)

# -file_type: either bam, tblat, or fastq

# -filter_or_count: either 'f' or 'c', determines whether the file is either filtered for contaminants using a pre-trained SVM model or simply counted for Cs and Ts

# -remove_or_tag: either 'r' or 't', indicating whether to remove reads that do not meet the threshold or to simply tag them.
#  => FASTA/Q files are tagged in their name with 'contam' if the read is contaminated
#  => BAM files have an added tag 'CM' with either 0 or 1
#  => TBLAT files have an added column labeled 'pass' with either 0 or 1

# -model: name of a joblib object containing a pre-trained model used for classification. Required if filter_or_count == 'f'.

# -output_dir: directory to save edited files to, optional
# -output_name: string used to name the edited filename(s)
# -ignore_motif: motif containing Cs to ignore the count of (e.g. CG, CHH, CHG), user can define their own motif where the counted C is in capitals
# 		+CCGG = ignore both cytosines in this motif
# 		+cCGG = only ignore the second cytosine
# 		+CcGG = only ignore the first cytosine
# -report_filename: string used to name the CSV file of a dataframe which contains the read name, C count, and T count for each read from the original file.
#  => if filter_or_count = 'f', this will also include another column, 'classification' which shows if the read was labeled as a contaminant or clean
#  => if it is a TBLAT file, this will also include the column 'species_contam_proba' which is the average score from 'classification' grouped by species.


def editfastq(r1_in, r2_in, model, remove_or_tag, count_or_filter, ignore_motif, output_name, output_dir, report_filename):

	list_of_rows = []
	row = {'read_id': '','C_count':0, 'T_count':0, 'classification':10000}

	if count_or_filter == 'f':
		svc = joblib.load(model)
		r1_out = output_dir+'/R1.'+output_name
		r2_out = output_dir+'/R2.'+output_name
		r1output = dnaio.open(r1_out, mode='w')
		r2output = dnaio.open(r2_out, mode='w')
		if remove_or_tag == 'r':
			r1_out_contam = output_dir+'/removed.R1.'+output_name
			r2_out_contam = output_dir+'/removed.R2.'+output_name
			r1output_contam = dnaio.open(r1_out_contam, mode='w')
			r2output_contam = dnaio.open(r2_out_contam, mode='w')

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

			# Count number of Ts using normal count()
			count_T = combined_seq.count("T")

			# Count number of Cs using the proper ignore_motif
			count_C = count_cs_fastaq(combined_seq, ignore_motif)

			# Add data to list for dataframe:
			row['C_count'] = count_C
			row['T_count'] = count_T
			row['read_id'] = record1.name

			if count_or_filter == 'f':
				if is_contam(svc, count_C, count_T):
					row['classification'] = 1 #classified as contaminant
					if remove_or_tag == 't': #if tagging, add 'contam' to name
						record1.name = 'contam'+record1.name[1:]
						record2.name = 'contam'+record2.name[1:]
						r1output.write(record1)
						r2output.write(record2)
					else: #if removing, write the removed reads into contam files
						r1output_contam.write(record1)
						r2output_contam.write(record2)
					cut+=1
				else:
					row['classification'] = 0 #classified as clean
					r1output.write(record1)
					r2output.write(record2)
			list_of_rows.append(row.copy())
			ent += 1

	if count_or_filter == 'f':
		print(f'Processed {ent} paired-end reads')
		if remove_or_tag =='t':
			print(f'Tagged {cut}({(cut/ent)*100:.2f}%) of reads')
		else:
			print(f'Removed {cut}({(cut/ent)*100:.2f}%) of reads')
		r1output.close()
		r2output.close()
		if remove_or_tag == 'r': #close files with removed reads
			r1output_contam.close()
			r2output_contam.close()
		#save reporting file with predictions
		pd.DataFrame(list_of_rows).to_csv(report_filename, sep = '\t', index=False)

	else:
		# If just counting, remove the predictions column from the dataframe
		df = pd.DataFrame(list_of_rows)
		df[['read_id','C_count','T_count']].to_csv(report_filename, sep = '\t', index=False)



def edit_bam(bam_file, model, remove_or_tag, count_or_filter, output_name, output_dir, ignore_motif, report_filename):
	num = 0 #total number of paired-end reads
	cut = 0 #number of paired-end reads cut from file

	# Initialize list to create dataframe
	list_of_rows = []
	row = {'read_id': '','C_count':0, 'T_count':0, 'classification':10000}

	bam = pysam.AlignmentFile(bam_file, 'rb') #read in bam file

	if count_or_filter =='f':
		# load model
		svc = joblib.load(model)
		#create output file
		bam_out = pysam.AlignmentFile(output_dir+'/'+output_name, 'wb', template = bam)
		if remove_or_tag == 'r': #create second output file for removed reads
			bam_out_contam = pysam.AlignmentFile(output_dir+'/removed.'+output_name, 'wb', template = bam)
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
			#if read1.query_name == read2.query_name: #reading pairs of reads together
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
				count_C = count_cs_bam(xm1+xm2, seq1+seq2, ignore_motif)

			else: #not overlapping, treat the two strings differently
				count_C = count_cs_bam(xm1+xm2, seq1+seq2, ignore_motif)

			# Count Ts using the full sequence:
			count_T = (seq1+seq2).count("T")

			# Add counts to dataframe:
			row['C_count'] = count_C
			row['T_counts'] = count_T
			row['read_id'] = read1.query_name

			if count_or_filter == 'f':
				#check if it is marked as a contaminant:
				if is_contam(svc, count_C, count_T):
					row['classification'] = 1 # read marked as contaminant
					cut+=1
					if remove_or_tag == 't':
						read1.tags += [("CM", 1)]
						read2.tags += [("CM", 1)]
						bam_out.write(read1)
						bam_out.write(read2)
					else: #if removing from the original file, write into the removed file
						bam_out_contam.write(read1)
						bam_out_contam.write(read2)
				else:
					read1.tags += [("CM", 0)]
					read2.tags += [("CM", 0)]
					bam_out.write(read1)
					bam_out.write(read2)
					row['classification'] = 0 # read marked as clean
			# Write into dataframe
			list_of_rows.append(row.copy())
			num+=1

	else: #treatment for single reads
		for read in bam:
			xm = read.get_tags()[2][1]
			count_C = count_cs_bam(xm, read.query_sequence, ignore_motif)
			count_T = read.query_sequence.count("T")

			row['C_count'] = count_C
			row['T_counts'] = count_T
			row['read_id'] = read.query_name

			if count_or_filter == 'f':
				#check if it is within the cutoff:
				if is_contam(svc, count_C, count_T):
					row['classification'] = 1 # read marked as contaminant
					if remove_or_tag == 't':
						read.tags += [("CM", 1)]
						bam_out.write(read)
					else:
						bam_out_contam.write(read)
					cut+=1
				else:
					row['classification'] = 0 # read marked as clean
					bam_out.write(read)

			# Write into dataframe
			list_of_rows.append(row.copy())
			num+=1
	print(f'Processed {num} paired-end reads')

	if count_or_filter =='f':
		if remove_or_tag =='t':
			print(f'Tagged {cut} ({(cut/num)*100:.2f}%) of reads')
		else:
			print(f'Removed {cut} ({(cut/num)*100:.2f}%) of reads')
		bam_out.close()
		pd.DataFrame(list_of_rows).to_csv(report_filename, sep = '\t', index=False)
		if remove_or_tag == 'r':
			bam_out_contam.close()

	else: # If counting, remove the predictions column from the dataframe
		df = pd.DataFrame(list_of_rows)
		df[['read_id','C_count','T_count']].to_csv(report_filename, sep = '\t', index=False)

	if os.path.exists('./temp.bam'):
		os.remove('temp.bam') # delete temporary sorted BAM file

####### Functions for reading in the tblat + R1 and R2 files

def rev_complement(string):
	d = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	reverse = string[::-1]
	rev_comp = ''
	for i in reverse:
		rev_comp += d.get(i, 'N')
	return(rev_comp)

def get_reference_seq(sseqid, ref_start, ref_end, strand, fasta):
	if strand == 'Plus/Plus':
		reference_sequence = fasta[sseqid][(ref_start-1):(ref_end)].seq
	if strand == 'Plus/Minus':
		reference_sequence = fasta[sseqid][(ref_end-1):(ref_start)].seq
		reference_sequence = rev_complement(reference_sequence)

	return(reference_sequence)

def gap_getter(read_seq, ref_seq):
	read_seqCT = read_seq.replace('C', 'T')
	ref_seqCT = ref_seq.replace('C', 'T')

	alignment = pairwise2.align.localxs(read_seqCT, ref_seqCT, -0.5, -0.1, gap_char="X", one_alignment_only=True)

	read_seqCT = alignment[0][0]
	ref_seqCT = alignment[0][1]
	read_seq_real=''

	read_split = read_seqCT.split('X')
	read_seq_good = ''
	i=0
	for sp in read_split:
		read_seq_good += read_seq[i:i+len(sp)]+'X'
		i=len(sp)
	read_seq_good = read_seq_good.rstrip('X')

	ref_split = ref_seqCT.split('X')
	ref_seq_good = ''
	i=0
	for sp in ref_split:
		ref_seq_good += ref_seq[i:i+len(sp)]+'X'
		i=len(sp)
	ref_seq_good = ref_seq_good.rstrip('X')

	return(read_seq_good, ref_seq_good)

def mismatch(entry, record1, record2, fasta):
	entry_split = entry.split('\t')
	read_id = entry_split[1]
	strand = entry_split[2]
	sseqid = entry_split[3]

	r1_start = int(entry_split[8])
	r1_end = int(entry_split[9])
	ref_r1_start = int(entry_split[10])
	ref_r1_end = int(entry_split[11])

	r2_start = int(entry_split[19])
	r2_end = int(entry_split[20])
	ref_r2_start = int(entry_split[21])
	ref_r2_end = int(entry_split[22])

	#frag_length = int(float(entry.strip().split('\t')[27]))
	#genome_length = int(float(entry.split('\t')[26]))

	R1_seq = record1.sequence[r1_start-1:r1_end]
	R2_seq = rev_complement(record2.sequence)[r2_start-1:r2_end]

	R1_ref_seq = get_reference_seq(sseqid, ref_r1_start, ref_r1_end, strand, fasta)
	R2_ref_seq = get_reference_seq(sseqid, ref_r2_start, ref_r2_end, strand, fasta)

	if (len(R1_seq) != len(R1_ref_seq)):
		R1_seq, R1_ref_seq = gap_getter(R1_seq, R1_ref_seq)

	if (len(R2_seq) != len(R2_ref_seq)):
		R2_seq, R2_ref_seq = gap_getter(R2_seq, R2_ref_seq)

	R1_mismatch = ''
	for bp, ref in zip(R1_seq, R1_ref_seq):
		if bp == 'C' and ref == 'C':
			R1_mismatch+='C'
		elif bp == 'T' and ref == 'C':
			R1_mismatch += 'T'
		else:
			R1_mismatch += '.'

	R2_mismatch = ''
	for bp, ref in zip(R2_seq, R2_ref_seq):
		if bp == 'C' and ref == 'C':
			R2_mismatch+='C'
		elif bp == 'T' and ref == 'C':
			R2_mismatch += 'T'
		else:
			R2_mismatch += '.'

	new_entry = entry.strip('\n').split('\t')+[R1_mismatch, R2_mismatch]
	return(new_entry)

def mismatch2(entry, record, fasta, read):
	entry_split = entry.split('\t')
	read_id = entry_split[1]
	strand = entry_split[2]
	sseqid = entry_split[3]
	if read == 'R1':
		r_start = int(entry_split[8])
		r_end = int(entry_split[9])
		ref_r_start = int(entry_split[10])
		ref_r_end = int(entry_split[11])
		R_seq = record.sequence[r_start-1:r_end]

	if read == 'R2':
		r_start = int(entry_split[19])
		r_end = int(entry_split[20])
		ref_r_start = int(entry_split[21])
		ref_r_end = int(entry_split[22])
		R_seq = rev_complement(record.sequence)[r_start-1:r_end]
	#frag_length = int(float(entry.strip().split('\t')[27]))
	#genome_length = int(float(entry.split('\t')[26]))

	R_ref_seq = get_reference_seq(sseqid, ref_r_start, ref_r_end, strand, fasta)

	if (len(R_seq) != len(R_ref_seq)):
		R_seq, R_ref_seq = gap_getter(R_seq, R_ref_seq)

	R_mismatch = ''
	for bp, ref in zip(R_seq, R_ref_seq):
		if bp == 'C' and ref == 'C':
			R_mismatch+='C'
		elif bp == 'T' and ref == 'C':
			R_mismatch += 'T'
		else:
			R_mismatch += '.'

	return(R_mismatch)

def process_read(l):
	tblat_file, read_file, metagenomic_reference, read, TMPNAME = l
	fasta = Fasta(metagenomic_reference)
	TMPNAME+=read
	with open(tblat_file) as f, dnaio.open(read_file) as fq, open(TMPNAME, 'w') as w:
		header = f.readline()
		i=0
		for entry, record in zip(f,fq):
			#rows.append(mismatch(entry, record1, record2, fasta))
			w.write(mismatch2(entry, record, fasta, read)+'\n')
			if i % 100000 ==0:
				print('processed another 100k in ' + str(time.time()-s))
				s = time.time()
				#print('processed '+ str(i) + ' in ' + str(time.time()-start))
			i+=1

def edit_tblat(tblat_file, r1, r2, R1_model, R2_model, count_or_filter, remove_or_tag, ignore_motif, output_name, output_dir, report_filename, metagenomic_reference, LUT,n_jobs):
	"""
	tblat files typically contain metagenomic data where we can use
	species level inferences to remove low-cytosine reads. It contains a mismatch string
	.....T....C....T...... where:
	. = non C nucleotide
	T = "unmethylated" (ie converted) cytosine
	C = "methylated" (ie unconverted) cytosine
	"""

	num = 0 #total number of paired-end reads
	cut = 0 #number of paired-end reads cut/tagged from file
	fasta = Fasta(metagenomic_reference)

	LUT = pd.read_csv(LUT, sep='\t').fillna(value=-1)
	LUT = LUT.rename(columns = {'Taxid':'taxid'})
	LUT = LUT[["taxid", "species"]]
	LUT['species'] = LUT['species'].astype(int)
	#LUT['taxid'] = LUT['taxid'].astype(str)
	# Check if matching R files are in the directory
	if not os.path.exists(r1):
		raise ValueError("Fie {} not in path".format(r1))
	elif not os.path.exists(r2):
		raise ValueError("Fie {} not in path".format(r2))

	print(f"Reading file {tblat_file}")

	start = time.time()
	s = start
	TMPNAME = 'TEMP_FILE' + ''.join(random.choice('0123456789ABCDEF') for i in range(16))

	#Parallel(n_jobs=n_jobs)(delayed(process_files)(l) for l in [[tblat_file, r1, metagenomic_reference, 'R1', TMPNAME], [tblat_file, r2, metagenomic_reference, 'R2', TMPNAME]])
	with open(tblat_file) as f, dnaio.open(r1) as fqr1, dnaio.open(r2) as fqr2, open(TMPNAME, 'w') as w:
		header = f.readline()
		i=0
		for entry, record1, record2 in zip(f,fqr1, fqr2):
			#rows.append(mismatch(entry, record1, record2, fasta))
			w.write('\t'.join(mismatch(entry, record1, record2, fasta))+'\n')
			if i % 100000 ==0:
				print('processed another 100k in ' + str(time.time()-s))
				s = time.time()
				#print('processed '+ str(i) + ' in ' + str(time.time()-start))
			i+=1
		print(i)

	df = pd.read_csv(TMPNAME, sep = '\t', header=None, names = ['taxid', 'qseqid', 'strand', 'sseqid', \
					'pident_R1', 'length_R1', 'mismatch_R1', 'gapopen_R1', 'qstart_R1', \
					'qend_R1', 'sstart_R1', 'send_R1', 'evalue_R1', 'bitscore_R1', 'qlen_R1', \
					'pident_R2', 'length_R2', 'mismatch_R2', 'gapopen_R2', 'qstart_R2', \
					'qend_R2', 'sstart_R2', 'send_R2', 'evalue_R2', 'bitscore_R2', 'qlen_R2', \
					'genome_len', 'effective_length', 'R1_mismatch', 'R2_mismatch'])
	df['qseqid'] = df['qseqid'].astype(str) + '-1'
	print(len(df.index))
	print(df.head())
	#rows = None # very lazy way to clear some memory
	df['R1_C'] = df['R1_mismatch'].str.count('C')/df['R1_mismatch'].str.len()
	df['R1_T'] = df['R1_mismatch'].str.count('T')/df['R1_mismatch'].str.len()
	df['R2_C'] = df['R2_mismatch'].str.count('C')/df['R2_mismatch'].str.len()
	df['R2_T'] = df['R2_mismatch'].str.count('T')/df['R2_mismatch'].str.len()

	print('C and T count done')

	df = pd.merge(df, LUT, on= 'taxid')
	print(len(df.index))
	if count_or_filter == 'f':
		# Load model

		R1_svc = joblib.load(R1_model)

		R2_svc = joblib.load(R2_model)
		# Make predictions, #0 = clean/true positive, 1 = contaminant/false positive:
		df['R1_classification'] = R1_svc.predict(df[['R1_C', 'R1_T']])
		df['R2_classification'] = R2_svc.predict(df[['R2_C', 'R2_T']])
		df['classification'] = (df['R1_classification'] + df['R2_classification'] >=1).astype(int)
		print('classified')

		species = df.drop_duplicates(subset = ['qseqid', 'species'])[["species", "classification"]] #random enough
		species = species.groupby(['species'])['classification'].mean().reset_index()
		species = species.rename(columns = {'classification' : 'species_contam_proba'})
		print(len(df.index))
		df = pd.merge(df, species, on = 'species')
		print(len(df.index))
		# Create new column showing whether or not the read is contaminated
		# 1= contaminated, average species score >0.5
		# 0= clean, average species score <= 0.5
		df['pass'] = ((df['species_contam_proba']>0.2) | (df['classification'] == 1)).astype(int)
		#df['pass'] = (df['species_contam_proba']>0.2).astype(int)

	# Print amounts of contaminants
		fp = len(df[df['pass']==1])
		perc = fp/len(df.index)*100

		print(f"Processed {len(df.index)} reads")
		if remove_or_tag == 't':
			print(f"Tagged {fp} ({perc}%) reads as contaminants")
			df.to_csv(output_dir+'/tagged.'+output_name, index=False, sep = '\t')
		else: # if remove, remove all entries where the species are contaminated
			print(f"Removed {fp} ({perc}%) contaminated reads")
			fail = df[df['pass']==1]
			pas = df[df['pass']==0]
			cols_in_order = ['qseqid', 'sseqid', \
				'pident_R1', 'length_R1', 'mismatch_R1', 'gapopen_R1', 'qstart_R1', \
				'qend_R1', 'sstart_R1', 'send_R1', 'evalue_R1', 'bitscore_R1', 'qlen_R1', \
				'pident_R2', 'length_R2', 'mismatch_R2', 'gapopen_R2', 'qstart_R2', \
				'qend_R2', 'sstart_R2', 'send_R2', 'evalue_R2', 'bitscore_R2', 'qlen_R2', \
				'genome_len', 'effective_length', 'R1_mismatch', 'R2_mismatch', \
				'R1_C', 'R1_T', 'R2_C', 'R2_T', 'R1_classification', 'R2_classification', \
				'classification', 'species_contam_proba', 'strand', 'taxid', 'species']
			fail[cols_in_order].to_csv(output_dir+'/removed.'+output_name, index=False, sep = '\t')
			pas[cols_in_order].to_csv(output_dir+'/'+output_name, index=False, sep = '\t')

			# Save counts to the report dataframe

			cols = ['qseqid','taxid','R1_C', 'R1_T', 'R2_C', 'R2_T', 'R1_classification', 'R2_classification', 'classification', 'species_contam_proba','pass']
			df[cols].to_csv(report_filename, sep = '\t', index=False)

	else:
		df[['qseqid','R1_C','R1_T', 'R2_C', 'R2_T', 'species']].to_csv(report_filename, sep = '\t', index=False)

def is_contam(model_obj, C_count, T_count):
	row = [[C_count, T_count]]
	return( int(model_obj.predict(row)[0]) )

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
	else:
		if ignore_motif == 'CG':
			ignore_cs = xm_str.count('c')+xm_str.count('C')

		else: #otherwise, use the sequence to count cs and motifs
			#number of upper case Cs = the number to count from the motif
			num_to_count = ignore_motif.count('C')
			ignore_motif = ignore_motif.upper() # change motif to upper case for counting
			# multiply the motif count by the number of Cs to count
			ignore_cs = num_to_count * ( sequence.count(ignore_motif) )

		return(sequence.count('C')-ignore_cs) # subtract out the motifs from the count



if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=textwrap.dedent("""\
	This script takes input reads from either a FASTQ, BAM, or TBLAT + FASTQ file and will run throughh one of two analyses:
	1) 'count': Moves through the file and counts Cs and Ts after merging paired-end reads, outputting a reporting CSV file containing a dataframe of counts associated with each read ID. The counting algorithm depends on the filetype:
		a) FASTQ: Cs are counted as normal but ignoring Cs as specified in the ignore_motif. Ts are counted as normal.
		b) BAM: Cs are counted using the XM string if there is no ignore_motif using the symbols cCxXhH. Otherwise, Cs are counted using the motif and the sequence itself. Ts are counted as normal.
		c) TBLAT: Cs and Ts are counted using the mismatch string created within the tblat code.
	2) 'filter': Moves through the files and counts Cs and Ts as described above, then uses this data in conjunction with a pre-trained SVM model to predict whether each read is clean or a contaminant. This information is shown in one of two ways:
		a) 'remove': Writes reads which were classified as contaminants to one file with the prefix 'removed.' and writes clean reads to the original output file.
		b) 'tag': Tags each read based on the filetype:
			i) FASTQ: adds 'contam' to the beginning of the read ID
			ii) BAM: adds a new flag 'CM' which has a value of either 0: clean or 1: contaminated
			iii) TBLAT: adds a column to the new tblat file called 'pass' which has a value of 0: clean or 1: contaminated.
		Note: both the FASTQ and BAM reads are filtered by reads individually based on the model's prediction. TBLAT reads are classified using the model, then grouped and averaged by species. If the average species score is > 0.5, all species reads are classified as contaminants. Otherwise, they are classified as clean.

	Example commands to run the program:
	"python3 filter3.py -fn BK_short_R1.fastq.gz BK_short_R2.fastq.gz BK_short.tblat -ft tblat -mo model1.joblib -cf f -rt t -o test1.tblat -r report_test1.csv -mr NCBIGenomes06.fna"
	=> filters the above filenames for contaminants by tagging using model1.joblib

	"python3 filter3.py -fn BK_short_R1.fastq.gz BK_short_R2.fastq.gz -ft fastq -cf c -r report_test1.csv"
	=> counts Cs and Ts from the above filenames and returns them to report_test1.csv
	"""))
	parser.add_argument('-fn',"--filenames", help="path to r1 r2 paired FASTQ, BAM, or tblat input filenames", nargs='+')
	parser.add_argument('-ft',"--file_type", type=str,help="type of input file", choices = ['bam','fastq','tblat'], nargs = 1)
	parser.add_argument('-mo',"--model", type=str,help="path to joblib file containing trained model for C count and T count contamination prediction, not required if counting only is selected", nargs = '+', default = None)
	parser.add_argument('-cf','--count_or_filter', type=str, help='Filter data in input or just return counts',nargs=1, choices=['f','c'])
	parser.add_argument('-rt','--remove_or_tag',type=str, help = 'Remove entry or tag, optional if choosing to count', nargs = '?', choices = ['r','t'], default = 'r')
	parser.add_argument('-m', "--ignore_motif", type=str, help = 'Motif including C of pattern to ignore within C count for filtering', nargs = '?', default = None )
	parser.add_argument('-d',"--output_dir", type=str, help = 'Path to output, optional', nargs = '?', default = os.getcwd())
	parser.add_argument('-o',"--output_name", type=str, help = 'Path for edited files, required for filtering',nargs='?', default=None)
	parser.add_argument('-r',"--report_filename", type=str, help = 'Name of file containing counts and contamination classification (if filter selected)', nargs = '?')
	parser.add_argument('-mr', "--metagenomic_reference", type=str, help = 'metagenomic reference used (for tblat)', nargs = '?', default = None)
	parser.add_argument('-LUT', "--lookup_table", type=str, help = 'metagenomic reference lookup table (for tblat)', nargs = '?', default = None)
	parser.add_argument('-j', "--n_jobs", type=int, default=1)

	args = parser.parse_args()
	rs = args.filenames
	file_type = args.file_type[0]
	R1_model = args.model
	remove_or_tag = args.remove_or_tag
	count_or_filter = args.count_or_filter[0]
	ignore_motif = args.ignore_motif
	output_dir = args.output_dir
	output_name = args.output_name
	report_filename = args.report_filename
	metagenomic_reference = args.metagenomic_reference
	LUT = args.lookup_table
	n_jobs = args.n_jobs
	#fasta=Fasta(metagenomic_reference)
	#pandarallel.initialize(nb_workers=30, progress_bar=)

	if len(R1_model) == 2:
		R2_model = R1_model[1]
		R1_model = R1_model[0]


	if file_type == 'fastq':
		# Import edited FLASH read overlapping code:
		combine = CDLL('./combine.so')
		combine_reads = combine.combine_reads
		combine_reads.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_int, c_int, c_float, c_bool]
		combine_reads.restype = c_char_p #make sure the output is a char

		if len(rs) != 2:
			raise ValueError('Incorrect number of FASTQ input files (require 2)')
		r1in = rs[0]
		r2in = rs[1]

		start = time.time()
		editfastq(r1in, r2in, model, remove_or_tag, count_or_filter, ignore_motif, output_name, output_dir, report_filename)
		end = time.time()
		print(f'Files {r1i} and {r2i} took {end-start:.3f} seconds to run')

	elif file_type == 'bam':
		if len(rs) != 1:
			raise ValueError('Incorrect number of BAM input files (require 1)')

		start = time.time()
		edit_bam(rs[0], model, remove_or_tag, count_or_filter, output_name, output_dir, ignore_motif, report_filename)
		end = time.time()
		print(f'File {bam_file} took {end-start:.3f} seconds to run')

	else: #will be tblat file
		if len(rs) != 3:
			raise ValueError('Incorrect number of TBLAT input files (require 3)')

		r1 = rs[0]
		r2 = rs[1]
		tblat = rs[2]
		start = time.time()
		edit_tblat(tblat, r1, r2, R1_model, R2_model, count_or_filter, remove_or_tag, ignore_motif, output_name, output_dir, report_filename, metagenomic_reference, LUT, n_jobs)
		end = time.time()
		print(f'File {tblat} took {end-start:.3f} seconds to run')
