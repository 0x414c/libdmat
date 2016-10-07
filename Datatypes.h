#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "Config.h"


/**
 \typedef	double_t entry_t

 \brief		Defines an alias representing the Matrix entry type.
 */
#ifdef DOUBLE_PRECISION
typedef double_t entry_t;
#else
typedef float_t entry_t;
#endif

/**
 \struct	_Mat_struct

 \brief	A matrix structure.
 */
struct _Mat_struct {
	/**
	 \brief	The entry_t**.
	 */
	entry_t **data;

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
	entry_t trace;

	/**
	 \brief	The determinant of a Matrix.
	 */
	entry_t det;

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
		The `Mat` 'itself' is a ptr to _Mat_struct.
 */
typedef struct _Mat_struct *Mat;
