#include <stdint.h>
#include <stddef.h>

#include "Vector.h"


bool nextPermutation (size_t *index, ptrdiff_t k, ptrdiff_t n) { //TODO:
	Assert$(index != NULL, "index should not be NULL.");

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
bool exists_d (double *c, size_t start, double value) {
	Assert$(c != NULL, "c should not be NULL.");

	for (size_t i = 0; i < start; i++) {
		if (equals_d(c[i], value)) {
			return true;
		}
	}

	return false;
}

bool exists_u (size_t *c, size_t start, size_t end, size_t value) {
	Assert$(c != NULL, "c should not be NULL.");

	for (size_t i = end; i < start; i++) {
		if (c[i] == value) {
			return true;
		}
	}

	return false;
}

/**
 \fn	int64_t *AllocVec_i (size_t size)

 \brief	Allocate vector.

 \param	size	The size.

 \return		null if it fails, else an int64_t pointer.
 */
int64_t *AllocVec_i (size_t size) {
	Assert$(size > 0, "size should be greater than 0");

	int64_t *v = NULL;
	v = ((int64_t *) malloc(sizeof(int64_t) * (size)));
	//v = ((int64_t *) calloc(Size, sizeof(int64_t)));
	Assert$(v != NULL, "Cannot allocate memory.");

	return v;
}

/**
 \fn	size_t *AllocVec_u (size_t size)

 \brief	Allocate vector.

 \param	size	The size.

 \return	null if it fails, else a size_t pointer.
 */
size_t *AllocVec_u (size_t size) {
	Assert$(size > 0, "size should be greater than 0");

	size_t *v = NULL;
	v = ((size_t *) malloc(sizeof(size_t) * (size)));
	//v = ((size_t *) calloc(Size, sizeof(size_t)));
	Assert$(v != NULL, "Cannot allocate memory.");

	return v;
}

/**
 \fn	long double *AllocVec_d (size_t size)

 \brief	Ld allocate vector.

 \param	size	The size.

 \return	null if it fails, else a double pointer.
 */
double *AllocVec_d (size_t size) {
	Assert$(size > 0, "size should be gtreater than 0");

	double *v = NULL;
	v = ((double *) malloc(sizeof(double) * size));
	//v = ((double *) calloc(Size, sizeof(double)));
	Assert$(v != NULL, "Cannot allocate memory.");

	return v;
}

/**
 \fn	void fillSequential_u(size_t *index, size_t size, size_t start)

 \brief	Fills index array.

 \param [in,out]	index	If non-null, * to index array.
 \param	size			 	The array size.
 \param	start			 	The starting value.
 */
void fillSequential_u (size_t *index, size_t size, size_t start) {
	Assert$(index != NULL, "index should not be NULL.");

	for (size_t i = 0; i < size; i++) {
		index[i] = start++;
	}

	return;
}
