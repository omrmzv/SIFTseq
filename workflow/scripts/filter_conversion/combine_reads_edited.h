#ifndef _FLASH_COMBINE_READS_H_
#define _FLASH_COMBINE_READS_H_

#include <stdbool.h>
#include <stddef.h>
#include <pthread.h>
#include <stdio.h>
#include <Python.h>
#define ARRAY_LEN(A) (sizeof(A) / sizeof((A)[0]))

#ifdef __GNUC__
#	if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
#	endif
#	define __noreturn __attribute__((noreturn))
#	define __format(type, format_str, args_start) \
			__attribute__((format(type, format_str, args_start)))
#	define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#	define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })
#	define inline inline __attribute__((always_inline))
#else
#	define __noreturn
#	define __format(type, format_str, args_start)
#	define max(a,b) (((a) > (b)) ? (a) : (b))
#	define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


/* Parameters for the core algorithm of FLASH.  See the help output for more
 * information.  */
struct combine_params {
	/* --min-overlap  */
	int min_overlap;

	/* --max-overlap  */
	int max_overlap;

	/* --max-mismatch-density  */
	float max_mismatch_density;

	/* --cap-mismatch-quals  */
	bool cap_mismatch_quals;

	/* --allow-outies  */
	bool allow_outies;
};

/* Result of a call to combine_reads()  */
enum combine_status {
	/* The reads could not be combined.  */
	NOT_COMBINED = 0,

	/* The reads were combined in "innie" orientation, like the following:
	 *
	 * ---------->
	 *     <------------
	 *
	 * (Note: read_2 is reverse complemented before the call to
	 * combine_reads()).  */
	COMBINED_AS_INNIE,

	/* The reads were combined in "outie" orientation, like the following:
	 *
	 * <----------
	 *     ------------>
	 *
	 * (Note: read_2 is reverse complemented before the call to
	 * combine_reads()).  */
	COMBINED_AS_OUTIE,
};

//enum combine_seq;
char * combine_reads(const char * restrict read_1, const char * restrict read_2,
	const char * restrict qual_1, const char * restrict qual_2, char * combined_seq, 
	int min_overlap, int max_overlap, float max_mismatch_density,
	bool allow_outies);

//int free_pointer(char * pointer);


#endif /* _FLASH_COMBINE_READS_H_  */

