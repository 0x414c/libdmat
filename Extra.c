#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "Const.h"
#include "Extra.h"


/**
 \fn	int in (long double *c, size_t s, long double x) 
 \brief	Check if element is exists in array.			 
 \date	20-May-14										 
 \param [in]	c		If non-null, the long double * to process.
 \param	s			 	The array size to process.
 \param	x			 	The element value.				 
 \return	1 if success, else 0.
 */
bool exists_d (long double *c, size_t s, long double x) {
	for (size_t i = 0; i < s; i++) {
		if (equal_d(c[i], x)) {
			return 1; 
		}
	}

	return 0;
}

bool exists_u (size_t *c, size_t s, size_t e, size_t x) {
	for (size_t i = e; i < s; i++) {
		if (c[i] == x) {
			return 1;
		}
	}

	return 0;
}

/**
\fn	int nextCombination (int *index, int k, int n)	 
\brief	Makes next combination over the index array while it can be made.  	 
\date	15-May-14															 
\param [in,out]	index		If non-null, * to index array.
\param	k				 	The K.
\param	n				 	The N.											 
\return	true if next combination is generated, false if no possible combinations left.
*/
bool nextCombination (size_t *index, ptrdiff_t k, ptrdiff_t n) { //TODO: 
	for (ptrdiff_t i = k - 1; i >= 0; --i) {
		if (index[i] < n - k + i + 1) {	//warning C4018: '<' : signed/unsigned mismatch
			++index[i];
			for (ptrdiff_t j = i + 1; j < k; ++j) {
				index[j] = index[j - 1] + 1;
			}
			return true;
		}
	}

	return false;
}

void fillIndex(size_t *index, size_t Size) {
	for (size_t i = 0; i < Size; i++) {
		index[i] = i + 1;
	}

	return;
}

int64_t *iAllocVec (size_t Size) {
	int64_t *v;
	v = ((int64_t *) malloc(sizeof(int64_t) * (Size)));
	//v = ((int *) calloc(Size, sizeof(int)));
	Assert(v != NULL, "Cannot allocate.");
	
	return v;
}

size_t *uAllocVec (size_t Size) {
	size_t *v;
	v = ((size_t *) malloc(sizeof(size_t) * (Size)));
	//v = ((int *) calloc(Size, sizeof(int)));
	Assert(v != NULL, "Cannot allocate.");
	
	return v;
}

long double *ldAllocVec (size_t Size) {
	long double *v;
	v = ((long double *) malloc(sizeof(long double) * Size));
	//v = ((long double *) calloc(Size, sizeof(long double)));
	Assert(v != NULL, "Cannot allocate.");
	
	return v;
}

bool equal (double a, double b) {
	return (a==b);
}
