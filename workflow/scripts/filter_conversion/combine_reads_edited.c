/*
 * combine_reads.c:  This file contains the code implementing the core algorithm
 * to combine reads in FLASH.
 */

/*
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013, 2014 Eric Biggers
 *
 * Edited for implementation in Python by Sami Smalling.
 */

 #include "combine_reads_edited.h"

 #include <assert.h>
 #include <stdlib.h>
 #include <string.h>
 #include <limits.h>
 #include <Python.h>





void compute_mismatch_stats(const char * restrict seq_1,
	 const char * restrict seq_2,
 	 const char * restrict qual_1,
 	 const char * restrict qual_2,
 	 int * restrict len_p,
 	 unsigned * restrict num_mismatches_ret,
 	 unsigned * restrict mismatch_qual_total_ret,
 	 bool is_fasta)
  {
 	int num_uncalled = 0;
 	unsigned num_mismatches = 0;
 	unsigned mismatch_qual_total = 0;
 	int len = *len_p;

 	for (int i = 0; i < len; i++) {
 		if (seq_1[i] == 'N' || seq_2[i] == 'N') {
 			num_uncalled++;
 		} else {
 			if (seq_1[i] != seq_2[i])  {
 				num_mismatches++;
				if (!is_fasta) { //if fastq file, use quality string
					mismatch_qual_total += min(qual_1[i], qual_2[i]);
				}
 			}
 		}
 	}


 	/*Return results in pointer arguments  */
	/*If the length of the quality strings is 0, then return null to
	mismatch_qual_total_ret to indicate that the file used is a FASTA */
	*mismatch_qual_total_ret = mismatch_qual_total;
	*num_mismatches_ret = num_mismatches;
 	*len_p -= num_uncalled;
 }


 #define NO_ALIGNMENT INT_MIN

 int * pair_align(const char * restrict read_1, const char * restrict read_2,
 	const char * restrict qual_1, const char * restrict qual_2,
 	   int min_overlap, int max_overlap, float max_mismatch_density,
 	   bool allow_outies, bool * was_outie)
 {

 	/* Best (smallest) mismatch density that has been found so far in an
 	 * overlap. */
 	float best_mismatch_density = max_mismatch_density + 1.0f;
 	float best_qual_score = 0.0f;
 	int best_position = NO_ALIGNMENT;
 	bool best_was_outie = false;
 	bool doing_outie = false;
 	int start;
 	int end;
 	int read_offset = 0;
 	int best_offset = 0;
 	int seq_len1 = strlen(read_1);
 	int seq_len2 = strlen(read_2);
  	bool is_fasta;

  if (strlen(qual_1) == 0 & strlen(qual_2)==0) {
		is_fasta = true;
	} else {
		is_fasta = false;
	}

 	/* Require at least min_overlap bases overlap, and require that the
 	 * second read is not overlapped such that it is completely contained in
 	 * the first read.  */
 	start = max(0, seq_len1 - seq_len2);
 	end = seq_len1 - min_overlap + 1;


 	if (start == 0) {
 		read_offset = (seq_len2 - seq_len1);

 		if (read_offset < 0) {
 			read_offset *= -1;
 		}
 	}

	for (int i = start; i < end; i++) {
		unsigned num_mismatches;
		unsigned mismatch_qual_total;
		int overlap_len = seq_len1 - i;

		/*This modification allows for engulf cases of the read*/
		compute_mismatch_stats(read_1 + i,
			read_2 + read_offset,
			qual_1 + i,
			qual_2 + read_offset,
			&overlap_len,
			&num_mismatches,
			&mismatch_qual_total,
			is_fasta);

  again:
		/*Logic to make sure read makes minimum requirements*/
		if (((!doing_outie && overlap_len >= min_overlap) || (doing_outie && overlap_len >= min_overlap))) {
			// if reads come from a FASTQ file, use qual strings
			if (!is_fasta) {
				float score_len = (float)min(overlap_len, max_overlap);
				float qual_score = mismatch_qual_total / score_len;
				float mismatch_density = num_mismatches / score_len;

				if (mismatch_density <= best_mismatch_density &&
				    (mismatch_density < best_mismatch_density ||
				     qual_score < best_qual_score))
				{
					best_qual_score       = qual_score;
					best_mismatch_density = mismatch_density;
					best_position         = i;
					best_offset           = read_offset;
					best_was_outie        = doing_outie;
				}
			if (read_offset != 0) {
				read_offset -= (i + 1);
				i--;
				}
			} else { // if reads come from FASTA file, pick randomly between ties
				float score_len = (float)min(overlap_len, max_overlap);
				float mismatch_density = num_mismatches / score_len;
        float qual_score = 0;

				if (mismatch_density <= best_mismatch_density)
				{
					best_mismatch_density = mismatch_density;
					best_position         = i;
					best_offset           = read_offset;
					best_was_outie        = doing_outie;
          best_qual_score       = 0;
				}
			if (read_offset != 0) {
				read_offset -= (i + 1);
				i--;
				}
			}
		}
	}

 	if (allow_outies) {
 		const char *tmp;
 		read_1 = read_2;
 		read_2 = tmp;
 		allow_outies = false;
 		doing_outie = true;
 		goto again;
 	}

 	int *position_and_offset = (int*)malloc(2 * sizeof(int));

 	position_and_offset[0] = best_position;
 	position_and_offset[1] = best_offset;

 	if (best_mismatch_density > max_mismatch_density) {
 		position_and_offset[0] = NO_ALIGNMENT;
 		position_and_offset[1] = NO_ALIGNMENT;
 	}

 	*was_outie = best_was_outie;
 	return position_and_offset;
  free(position_and_offset);
 }


char * generate_combined_read(const char * restrict read_1,
 	const char * restrict read_2,
 	const char * restrict q_1,
 	const char * restrict q_2,
  char * combined_seq,
 	int overlap_begin,
 	int read_offset,
 	bool was_outie)
 {
 	int seq_len1 = strlen(read_1);
 	int seq_len2 = strlen(read_2);

	bool is_fasta = false;
	if (strlen(q_1)==0 & strlen(q_2)==0){
		is_fasta = true;
	}

 	/* Length of the overlapping part of two reads.  */

 	int overlap_len = seq_len1 - overlap_begin;
 	/* Length of the part of the second read not overlapped with the first
 	 * read.  */
 	int remaining_len = seq_len2 - overlap_len ;
 	int combined_seq_len = seq_len1 + remaining_len;


  const char * restrict seq_1 = read_1;
  const char * restrict seq_2 = read_2;
  const char * restrict qual_1 = q_1;
  const char * restrict qual_2 = q_2;


 	if (was_outie) {
    
 		//Case of outie, Engulf case
 		//So, read 2 is engulfed by read 1
 		if (read_offset != 0) {
			//switch the two reads
		    seq_1 = read_2;
		    seq_2 = read_1;
		    qual_1 = q_2;
		    qual_2 = q_1;
			
			overlap_begin = read_offset;
 			remaining_len = 0;
 			combined_seq_len = seq_len2+read_offset;
 			overlap_len = seq_len2;
 		} else {
 			combined_seq_len = seq_len1 - overlap_begin;
 			seq_1 += overlap_begin;
 			overlap_begin = 0;
 			remaining_len = 0;
 		}

     /*Typical Case*/
   } else {
 		//Innie, engulf case
      if (read_offset != 0) {
        //R1 is shorter than R2 (take first part of read R2)

        /*ignore first part of read 2*/
        seq_2 += read_offset;
        qual_2 += read_offset;

        /*don't take any part of the first read*/
        overlap_begin = 0;

        /*Seq_2 minus first part of read 1*/
        combined_seq_len = seq_len2 - read_offset;
 			  remaining_len = seq_len2 - seq_len1 - read_offset;

 		} 
 	}


    char * combined_seq_ptr = combined_seq;
    char * combined_qual_ptr = (char *)malloc((combined_seq_len+1)*sizeof(char));
    char * combined_qual = combined_qual_ptr;
    

 	/* Copy the beginning of read 1 (not in the overlapped region).  */

 	while (overlap_begin--) {
 		*combined_seq++ = *seq_1++;
 		*combined_qual++ = *qual_1++;
 	}


  /* Copy the overlapped region.  */
 	while (overlap_len-- > 0) {
 		if (*seq_1 == *seq_2) {
 			/* Same base in both reads. If FASTQ file take the higher quality read. */
 			*combined_seq = *seq_1;
			if (!is_fasta){
				*combined_qual = max(*qual_1, *qual_2);
			}
 		} else {
 			/* Different bases in the two reads. If its a FASTA file, pick randomly.
			If FASTQ file, pick the higher quality read.
 			 */
			 if (!is_fasta) {
				 *combined_qual = max(abs(*qual_1 - *qual_2), 2);

	  			if (*qual_1 > *qual_2) {
	  				*combined_seq = *seq_1;
	  			} else if (*qual_1 < *qual_2) {
	  				*combined_seq = *seq_2;
	  			} else {
	  				/* Same quality value; take the base from the
	  				 * first read if the base from the second read
	  				 * is an 'N'; otherwise take the base from the
	  				 * second read. */

	  				if (*seq_2 == 'N'){
              *combined_seq = *seq_1;
            } else {
              *combined_seq = *seq_2;
            }
	  			}
			 }

       else { //is a fasta file
         if (*seq_1 =='N' || *seq_2 == 'N') { //if one of them is an N take the other
           if (*seq_2 == 'N'){
             *combined_seq = *seq_1;
           } else {
             *combined_seq = *seq_2;
           }
         } else { // otherwise, pick randomly
           srand(time(0)); //set seed for random number
  				 double num = (double)rand() / (double)RAND_MAX; //pick random number between 0 and 1
  				 if (num <= 0.5) {
  					 *combined_seq = *seq_1;
  				 } else {
  					 *combined_seq = *seq_2;
  				 }
         }
			 }
 		}
 		combined_seq++;
 		combined_qual++;
 		seq_1++;
 		seq_2++;
 		qual_1++;
 		qual_2++;
 	}

 	/* Copy the end of read 2 (not in the overlapped region).  */
  while (remaining_len--) {
    if (read_offset != 0 && was_outie) {
      *combined_seq++ = *seq_2++;
      *combined_qual++ = *qual_2++;
    } else {
      *combined_seq++ = *seq_2++;
      *combined_qual++ = *qual_2++;
    }
  }

  while ( strlen(combined_seq_ptr) > combined_seq_len) {
    combined_seq[strlen(combined_seq)-1] = '\0';
      	}

  return(combined_seq_ptr);
 }


char * combine_reads(const char * restrict read_1, const char * restrict read_2,
  const char * restrict qual_1, const char * restrict qual_2,
  char * combined_seq,
  int min_overlap, int max_overlap, float max_mismatch_density,
  bool allow_outies)
 {
  int overlap_begin, read_offset, *overlap_and_offset;
  enum combine_status status;
  bool was_outie = false;

  /* Do the alignment.  */
  overlap_and_offset = pair_align(read_1, read_2, qual_1, qual_2,
 	 min_overlap, max_overlap, max_mismatch_density,
 	 allow_outies, &was_outie);

  overlap_begin = overlap_and_offset[0];
  read_offset = overlap_and_offset[1];

  /*
 	* If overlap_begin == NO_ALIGNMENT, then no sufficient overlap between
 	* the reads was found.
 	*
 	* If !@was_outie, then the pair forms an "innie" overlap, and
 	* overlap_begin is the 0-based position in read_1 at which read_2
 	* begins.  (Shown below with read 2 already reverse complemented!)
 	*
 	*		      0	        overlap_begin
 	*	        |         |
 	*	Read 1: ------------------>
 	*	Read 2:           ---------------------->
 	*
 	* If @was_outie, then the pair forms an "outie" overlap, and
 	* overlap_begin is the 0-based position in read_2 at which read_1
 	* begins. (Shown below with read 2 already reverse complemented!)
 	*
 	*	        0         overlap_begin
 	*	        |         |
 	*	Read 2: ------------------>
 	*	Read 1:           ---------------------->
 	*/

  if (overlap_begin == NO_ALIGNMENT){
	combined_seq = "Not Combined";
 	return combined_seq;
  }

  if (!was_outie) {
 	 status = COMBINED_AS_INNIE;
  } else {
 	 const char *tmp;

 	 /* Simplify generation of the combined read by turning the outie
 		* case into the innie case.  */

 	 tmp = read_1;
 	 read_1 = read_2;
 	 read_2 = tmp;

 	 status = COMBINED_AS_OUTIE;
 	 /*
 		* Now it's just:
 		*
 		*		0	  overlap_begin
 		*	        |         |
 		*	Read 1: ------------------>
 		*	Read 2:           ---------------------->
 		*
 		* The same as the "innie" case.
 		*/
  }

 		/* Fill in the combined read.  */
 	combined_seq = generate_combined_read(read_1, read_2, qual_1, qual_2,
      combined_seq, overlap_begin, read_offset, was_outie);
    return combined_seq;
 }



char * combine_reads(const char * restrict, const char * restrict,
  const char * restrict, const char * restrict, char *,
  int, int, float, bool);

