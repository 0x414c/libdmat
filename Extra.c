#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "Const.h"
#include "Extra.h"


/**
\fn	bool nextCombination (size_t *index, ptrdiff_t k, ptrdiff_t n)
\brief	Makes next combination over the index array while it can be made.
\date	15-May-14
\param [in,out]	index	If non-null, * to index array.
\param	k				The K.
\param	n				The N.
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


#pragma region "Search in array"

/**
 \fn	bool exists_d (long double *c, size_t s, long double x)

 \brief	Check$ if element x is exists in array c of size s.

 \date	20-May-14

 \param [in]	c	If non-null, the * to array process.
 \param	start		The array size to process.
 \param	value		The element value.

 \return	true if success, else false.
 */
bool exists_d (long double *c, size_t start, long double value) {
	for (size_t i = 0; i < start; i++) {
		if (equals_ld(c[i], value)) {
			return true;
		}
	}

	return false;
}

bool exists_u (size_t *c, size_t start, size_t end, size_t value) {
	for (size_t i = end; i < start; i++) {
		if (c[i] == value) {
			return true;
		}
	}

	return false;
}
#pragma endregion "Search in array"


#pragma region "Allocation routines"

/**
 \fn	int64_t *AllocVec_i (size_t Size)

 \brief	Allocate vector.

 \param	Size	The size.

 \return	null if it fails, else an int64_t*.
 */
int64_t *AllocVec_i(size_t Size) {
	int64_t *v;
	v = ((int64_t *) malloc(sizeof(int64_t) * (Size)));
	//v = ((int *) calloc(Size, sizeof(int)));
	Assert$(v != NULL, "Cannot allocate.");

	return v;
}

/**
 \fn	size_t *AllocVec_u (size_t Size)

 \brief	Allocate vector.

 \param	Size	The size.

 \return	null if it fails, else a size_t*.
 */
size_t *AllocVec_u(size_t Size) {
	size_t *v;
	v = ((size_t *) malloc(sizeof(size_t) * (Size)));
	//v = ((int *) calloc(Size, sizeof(int)));
	Assert$(v != NULL, "Cannot allocate.");

	return v;
}

/**
 \fn	long double *AllocVec_ld (size_t Size)

 \brief	Ld allocate vector.

 \param	Size	The size.

 \return	null if it fails, else a double*.
 */
long double *AllocVec_ld(size_t Size) {
	long double *v;
	v = ((long double *) malloc(sizeof(long double) * Size));
	//v = ((long double *) calloc(Size, sizeof(long double)));
	Assert$(v != NULL, "Cannot allocate.");

	return v;
}
#pragma endregion "Allocation routines"

/**
 \fn	void fillSequential_u(size_t *index, size_t size, size_t start)

 \brief	Fills index array.

 \param [in,out]	index	If non-null, * to index array.
 \param	size			 	The array size.
 \param	start			 	The starting value.
 */
void fillSequential_u(size_t *index, size_t size, size_t start) {
	for (size_t i = 0; i < size; i++) {
		index[i] = start++;
	}

	return;
}

int64_t GCD_euclid (int64_t a, int64_t b) {
	int64_t t = 0;
	while (b) {
		t = a;
		a = b;
		b = t % b;
	}

	return llabs(a);
}
