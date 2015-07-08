#include <stdint.h>
#include <stddef.h>

#include "Vector.h"


bool nextPermutation (size_t *index, ptrdiff_t k, ptrdiff_t n) { //TODO:
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

/**
 \fn	int64_t *AllocVec_i (size_t size)

 \brief	Allocate vector.

 \param	size	The size.

 \return		null if it fails, else an int64_t pointer.
 */
int64_t *AllocVec_i (size_t size) {
    int64_t *v = NULL;
    v = ((int64_t *) malloc(sizeof(int64_t) * (size)));
    //v = ((int64_t *) calloc(Size, sizeof(int64_t)));
    Assert$(v != NULL, "Cannot allocate.");

    return v;
}

/**
 \fn	size_t *AllocVec_u (size_t size)

 \brief	Allocate vector.

 \param	size	The size.

 \return	null if it fails, else a size_t pointer.
 */
size_t *AllocVec_u (size_t size) {
    size_t *v = NULL;
    v = ((size_t *) malloc(sizeof(size_t) * (size)));
    //v = ((size_t *) calloc(Size, sizeof(size_t)));
    Assert$(v != NULL, "Cannot allocate.");

    return v;
}

/**
 \fn	long double *AllocVec_ld (size_t size)

 \brief	Ld allocate vector.

 \param	size	The size.

 \return	null if it fails, else a double pointer.
 */
long double *AllocVec_ld (size_t size) {
    long double *v = NULL;
    v = ((long double *) malloc(sizeof(long double) * size));
    //v = ((long double *) calloc(Size, sizeof(long double)));
    Assert$(v != NULL, "Cannot allocate.");

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
    for (size_t i = 0; i < size; i++) {
        index[i] = start++;
    }

    return;
}
