#pragma once

#include <stddef.h>
#include <stdbool.h>


/**
 \typedef	double entry_t

 \brief	Defines an alias representing the Matrix entry type.
 */
typedef double entry_t; //TODO:

/**
 \struct	_Mat_struct

 \brief	A matrix structure.
 */
struct _Mat_struct {
	/**
	 \brief	The entry_t**.
	 */
	entry_t **a;	

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
	 \brief	true if this Matrix is singular.
	 */
	bool isSingular;

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
	int permutationSign;

	/**
	 \brief	The trace of Matrix (sum of all the entries on main diagonal).
	 */
	entry_t trace;

	/**
	 \brief	The determinant of a Matrix.
	 */
	entry_t det;
};

/**
 \typedef	struct _Mat_struct *Mat

 \brief	Defines an alias representing the matrix.
		The Matrix 'itself' is a ptr to _Mat_struct.
 */
typedef struct _Mat_struct *Mat;
