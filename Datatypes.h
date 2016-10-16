#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "Config.h"


/**
 \typedef	double_t entry_type

 \brief		Defines an alias representing the Matrix entry type.
 */
#ifdef WITH_DOUBLE
typedef double_t entry_type;
#else //WITH_DOUBLE
typedef float_t entry_type;
#endif //WITH_DOUBLE

/**
 \struct	_mat_s

 \brief	A matrix structure.
 */
struct _mat_s {
	/**
	 \brief	The entry_type**.
	 */
	entry_type **data;

	/**
	 \brief	Number of rows.
	 */
	size_t rowsCount;

	/**
	 \brief	Number of columns.
	 */
	size_t colsCount;

	/**
	 \brief	The rank of Matrix.
	 */
	size_t rank;

	/**
	 \brief	The trace of Matrix (sum of all the entries on main diagonal).
	 */
	entry_type trace;

	/**
	 \brief	The determinant of a Matrix.
	 */
	entry_type det;

	/**
	 \brief	true if this Matrix is singular.
	 */
	//TODO: to bitfield
	bool isSingular;

	/**
	 \brief	true if this Matrix is identity matrix.
	*/
	//TODO: to bitfield
	bool isIdentity;

	/**
	 \brief	true if this Matrix is Symmetric positive-definite.
	 */
	bool isSPD;

	/**
	 \brief	true if this Matrix is rank deficient.
	 */
	bool isRankDeficient;

	/**
	 \brief	The permutation sign (switches from 1 to -1 every time when rows being swapped).
	 */
	int8_t permutationSign;
};

/**
 \typedef	struct _Mat_struct *Mat

 \brief	Defines an alias representing the matrix.
		The `Mat` 'itself' is a * to `_mat_s`.
 */
typedef struct _mat_s *Mat;
