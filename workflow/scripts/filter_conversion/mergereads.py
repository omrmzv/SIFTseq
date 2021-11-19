

import subprocess
import os
import dnaio
import string
import time
from Bio.Seq import Seq
			
def is_match(bp1, bp2):
	matches = {'A':['A','N'], 'T':['T','N'], 'G':['G','N'], 'C':['C','N'], 'N':['A','T','G','C','N']}
	return( bp2 in matches[bp1] )
			
def find_best_overlap(f_entry, r_entry, file_type, min_overlap, qual):
	''' 
	Takes a forward and reverse entry from dnaio-read FASTA or FASTQ file and returns dictionary with the 
score of the overlap, the locations of the start of the forward and end of the reverse sequences, the quality score of the overlap, and both sequence segments of overlap. 

	A score of 0 means that there are no mismatches between the two segments
	'''

	# Save reads and quality scores
	r1f = f_entry.sequence
	r1r = str(Seq(r_entry.sequence).reverse_complement()) # Take reverse complement of sequence
	q1f = f_entry.qualities
	q1r = r_entry.qualities[::-1] #reverse the order of the quality scores


	if len(r1f)>len(r1r): # If the forward section is larger, start from the middle of the forward segment
						  # and end at the length of the forward segment
		start_f = len(r1f)-len(r1r)-1
		end_r = len(r1f)-1
	else: # if the reverse is longer, start at 0 and end at the length of the reverse segment
		start_f = 0
		end_r = len(r1r)-1
		

	# Initialize dictionary of the best score, location of reverse strand, quality score, and segments matching
	best = {
	'score':1000000000000, 
	'olap':0, # locations of the start of the forward strand and end of reverse strand that lines the two up
	'q':0 , # quality score of best overlap for comparison. 
	'segments':('','') } # tuple of forward, reverse sequences

	#Calculate average length not including Ns
	good_overlap = (len(r1f)+len(r1r) - r1f.count('N') - r1r.count('N')) / 2
	
	while good_overlap >= min_overlap: 
		# Get overlapping sections: end of forward is always the length, start of reverse is always 0
		f = r1f[start_f:len(r1f)-1]
		r = r1r[0:end_r] 
		lengthf = len(f)-f.count('N')
		lengthr = len(r)-r.count('N') # Get length not counting Ns
		
		if lengthf<10 or lengthr<10: #If there are less than 10 non-N bps, skip the rest of the code
			continue 
		
		# Calculate the score = number of mismatches / total length
		mismatches =  [i for i in range(len(f)) if not is_match(f[i],r[i])]
		score = len(mismatches) / (lengthf + lengthr)
		
		# Keep the longest overlap for which there are no mismatches and there are more than 10 non-Ns 
		if score == 0 and lengthf > 10: 
			best['score'] = score
			best['olap'] = (start_f, end_r)
			#best['q'] = 'NA' # because there are no mismatches there is no quality score to be saved
			best['segments'] = (f,r)
			break
		elif score < best['score'] and lengthf > 10: # Keep this overlap as the running best
			best['score'] = score
			best['olap'] = (start_f, end_r)
			best['segments'] = (f,r)
		
		elif score == best['score']: # Compare q scores to determine best
			# Calculate the average quality score of all mismatches but only for FASTQ file types: 
			# For FASTQ files with quality scores:
			if file_type == 'fq':
				q_score = (sum([qual[q1f[i]] for i in mismatches]) + sum([qual[q1r[i]] for i in mismatches])) / (4*len(mismatches))
				if  q_score > best['q']:
					best['score'] = score
					best['olap'] = (start_f, end_r)
					best['q'] = q_score
					best['segments'] = (f,r)
		
		# Move one bp over:
		start_f+=1
		end_r = end_r-1
		
		# Save number of non-N bps:
		good_overlap = (len(f)+len(r) - f.count('N') - r.count('N')) / 2
	return(best)


R1_file = 'R1short.fq.gz'
R2_file = 'R2short.fq.gz'

# Create dictionary of ASCII characters to quality scores:
qual = {}
for s in string.printable:
	qual[s] = ord(s)-33
	

i = 0
perfect_match = 0
start = time.time()
with dnaio.open(R1_file) as r1, dnaio.open(R2_file) as r2:
	for entry1, entry2 in zip(r1, r2):
		best = find_best_overlap(entry1, entry2, 'fq', 10, qual)
		if i < 50:
			if best['score'] == 0:
				perfect_match+=1
		else:
			break
		i+=1
end = time.time()
print(f'Ran in {end-start:.2f} seconds')
print(f'There were {perfect_match} matches out of {i}')
			
		
		
		
		




