#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <limits.h>

#include "Gauss.h"
#include "Matrix.h"
#include "MatrixOps.h"
#include "Maths.h"
#include "Extras.h"
#include "Vector.h"


#pragma region "Determinant computation"

/**
 \fn		double Det_Gauss (Mat A)

 \brief		Calculates matrix determinant by Gauss' method for n&gt;3 and with rule of Sarrus for n&lt;=3.

 \date		13-May-14

 \param	A	The Mat to process.

 \return	Determinant value.
 */
entry_t Det_Gauss (Mat A) {
	Assert$(IsSquare$(A), "Matrix A should be square.");

	Mat T = NULL;
	entry_t **a = A->mat;
	entry_t det = 1.0;

	switch (A->rowsCount) {
		case 1:
			det = (a[0][0]);
			break;
		case 2:
			det = (((a[0][0] * a[1][1]) - (a[0][1] * a[1][0])));
			break;
		case 3:
			det = (
				a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
				a[0][1] * (a[2][2] * a[1][0] - a[1][2] * a[2][0]) +
				a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
			);
			break;
		case 4:
			det = (
				a[0][0] * (
					a[1][1] * (a[2][2] * a[3][3] - a[2][3] * a[3][2]) -
					a[2][1] * (a[1][2] * a[3][3] - a[1][3] * a[3][2]) +
					a[3][1] * (a[1][2] * a[2][3] - a[1][3] * a[2][2])
				) -
				a[1][0] * (
					a[0][1] * (a[2][2] * a[3][3] - a[2][3] * a[3][2]) -
					a[2][1] * (a[0][2] * a[3][3] - a[0][3] * a[3][2]) +
					a[3][1] * (a[0][2] * a[2][3] - a[0][3] * a[2][2])
				) +
				a[2][0] * (
					a[0][1] * (a[1][2] * a[3][3] - a[1][3] * a[3][2]) -
					a[1][1] * (a[0][2] * a[3][3] - a[0][3] * a[3][2]) +
					a[3][1] * (a[0][2] * a[1][3] - a[0][3] * a[1][2])
				) -
				a[3][0] * (
					a[0][1] * (a[1][2] * a[2][3] - a[1][3] * a[2][2]) -
					a[1][1] * (a[0][2] * a[2][3] - a[0][3] * a[2][2]) +
					a[2][1] * (a[0][2] * a[1][3] - a[0][3] * a[1][2])
				)
			);
			break;
		default:
			T = DeepCopy(A);

			toRowEchelonForm_reference(T);

			if (!(T->isSingular)) {
				for (size_t i = 0; i < A->rowsCount; i++) {
					det *= T->mat[i][i];
				}
				det *= T->permutationSign;
			} else {
				det = 0.0;
				A->isSingular = true;
			}
			freeMat$(T);
            A->det = det;

            return det;
	}

	if (iszero(det)) { A->isSingular = true; }
	A->det = det;

	return det;
}

/**
 \fn		double Det_Bareiss (Mat A)

 \brief		Computes Matrix determinant by Bareiss' algorithm.
			The Bareiss Algorithm is a fraction-free method for determinant computation.
			However, it can also be thought of as a sophisticated form of row reduction.
			Note that the divisions computed at any step are exact; thus Bareiss' Algorithm is
			indeed fraction-free. Entry a[n][n] is the determinant of A (after `pivoting` and `main` steps).

 \param	A	The Mat to process.

 \return	Determinant of Matrix A.
 */
//TODO: move to another file
entry_t Det_Bareiss (Mat A) {
  Assert$(IsSquare$(A), "Matrix A should be square.");

	Mat T = DeepCopy(A); //TODO: replace w/ Copy()
	entry_t **t = T->mat;

	// Pivotize
	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		for (size_t i = k; i < A->rowsCount; i++) {
			if (abs(t[i][k]) > abs(t[pivot][k])) {
				pivot = i;
			}
		}
		if (pivot != k) {
			entry_t *a_pivot = t[pivot];
			t[pivot] = t[k];
			t[k] = a_pivot;
			T->permutationSign *= -1; //-V127
		}
	}

	// Bareiss algorithm main step
	for (size_t i = 0; i < T->rowsCount - 1; i++) {
		Check$(isnotzero(t[i][i]), "Singular matrix.");
		for (size_t j = i + 1; j < T->rowsCount; j++) {
			for (size_t k = i + 1; k < T->rowsCount; k++) {
				t[j][k] = t[j][k] * t[i][i] - t[j][i] * t[i][k];
				if (i != 0) {
					t[j][k] /= t[i-1][i-1];
				}
			}
		}
	}

	entry_t det = t[T->rowsCount - 1][T->rowsCount - 1] * T->permutationSign;

	freeMat$(T);

	return det;
}
#pragma endregion "Determinant computation"


#pragma region "Transforming routines"

/**
 \fn		void toRowEchelonForm (Mat A)

 \brief		Transforms matrix A into a row echelon form. Optimized implementation.

 \param	A	The Matrix to process.
 */
void toRowEchelonForm (Mat A) {
	entry_t **a = A->mat;

	// W/ immediate rows swapping
	for (size_t k = 0; k < A->rowsCount; k++) {
		// Pivotize
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (abs(a[i][k]) > abs(a[k][k]))
			if (isnotzero(a[i][k])) {
				entry_t *a_i = a[i];
				a[i] = a[k];
				a[k] = a_i;
				A->permutationSign *= -1; //-V127
			}
		}
		// Eliminate
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (isnotzero(a[k][k])) {
				entry_t factor = a[i][k] / a[k][k];
				for (size_t j = k; j < A->colsCount; j++) {
					a[i][j] -= factor * a[k][j];
				}
			} else {
				A->isSingular = true;
				Check$(A->isSingular == false, "Singular matrix.");
				return;
			}
		}
	}

	return;
}

/**
 \fn		void toRowEchelonForm_reference (Mat A)

 \brief		Transforms Matrix A into a row echelon form.
 			Reference implementation.

 \param	A	The Matrix to process.
 */
void toRowEchelonForm_reference (Mat A) {
	entry_t **a = A->mat;

	// Reference implementation of pivoting algorithm
	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		// Find pivot
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (abs(a[i][k]) > abs(a[pivot][k])) {
				pivot = i;
			}
		}
		if (iszero(a[pivot][k])) {
			A->isSingular = true;
			Check$(A->isSingular == false, "Singular matrix.");
			return;
		}
		// Swap rows
		if (pivot != k) {
			entry_t *a_pivot = a[pivot];
			a[pivot] = a[k];
			a[k] = a_pivot;
			A->permutationSign *= -1; //-V127
		}
		// Calculate factor, eliminate k-th column
		if (isnotzero(a[k][k])) {
			for (size_t i = k + 1; i < A->rowsCount; i++) {
				entry_t f = a[i][k] / a[k][k];
//				a[i][k] /= a[k][k];
				for (size_t j = k; j < A->colsCount; j++) {
					a[i][j] -= f * a[k][j];
//					a[i][j] -= a[i][k] * a[k][j];
//					a[i][j] -= a[k][j] * (a[i][k] / a[k][k]);
				}
			}
		}
	}

	return;
}

/**
 \fn		void toReducedRowEchelonForm (Mat A)

 \brief		Transforms matrix into reduced row echelon form (aka row canonical form).
			The reduced row echelon form of A is unique, the pivot positions are
			uniquely determined and do not depend on whether or not row interchanges
			are performed in the reduction process.

 \date		13-May-14

 \param	A	The Matrix to process.
 */
void toReducedRowEchelonForm (Mat A) {
	entry_t **a = A->mat;
	size_t pivotCol = 0;
	size_t rowsCount = A->rowsCount, colsCount = A->colsCount;
	size_t i, j, k, r;

	for (r = 0; r < rowsCount; r++) {
		if (colsCount < pivotCol + 1) {
			break;
		}
		i = r;
		while (iszero(a[i][pivotCol])) {
			i++;
			if (i == rowsCount) {
				i = r;
				pivotCol++;
				if (colsCount == pivotCol) {
					pivotCol--;
					break;
				}
			}
		}
		for (j = 0; j < colsCount; j++) {
			swap(a[i][j], a[r][j]);
		}
		if (isnotzero(a[r][pivotCol])) {
			entry_t divisor = a[r][pivotCol];
			for (j = 0; j < colsCount; j++) {
				a[r][j] /= divisor;
			}
		}
		for (j = 0; j < rowsCount; j++) {
			if (j != r)	{
				entry_t sub = a[j][pivotCol];
				for (k = 0; k < colsCount; k++) {
					a[j][k] -= (sub * a[r][k]);
				}
			}
		}
		pivotCol++;
	}

	return;
}
#pragma endregion "Transforming routines"


#pragma region "Solving routines"

/**
 \fn		Mat Solve_GaussJordan (Mat A, Mat B)

 \brief		Solves system of linear equations using Gauss-Jordan method.

 \param	A	Coefficients matrix.
 \param	B	Right hand side as column-vector.

 \return	Solution as column-vector.
 */
Mat Solve_GaussJordan (Mat A, Mat B) {
	Assert$(B->colsCount == 1, "Matrix B must be column-vector.");
	Assert$(A->rowsCount == B->rowsCount, "Number of equations doesn't equal to number of unknown variables.");

	Mat X = AllocMat(A->rowsCount, 1);
	Mat AU = DeepCopy(A);

	concat(AU, B);
    entry_t **au = AU->mat;
    entry_t **x = X->mat;

	toReducedRowEchelonForm(AU);

	for (size_t i = 0; i < X->rowsCount; i++) {
		x[i][0] = au[i][AU->colsCount - 1];
	}

	freeMat$(AU);

	return X;
}

/**
 \fn		Mat Solve_Gauss (Mat A, Mat B)

 \brief		Solves system of linear equations using Gauss elimination.

 \param	A	Coeffs matrix.
 \param	B	Right hand side as column-vector.

 \return	Solution as column-vector.
 */
Mat Solve_Gauss (Mat A, Mat B) {
	Assert$(B->colsCount == 1, "Matrix B must be column-vector.");
	Assert$(A->rowsCount == B->rowsCount, "Number of equations must be equal to number of unknown variables.");

	Mat X = AllocMat(B->rowsCount, B->colsCount);
	Mat AU = DeepCopy(A);

	concat(AU, B);

	//Forward step (elimination with row pivoting)
	toRowEchelonForm_reference(AU);

	//Back-substitution
	for (ssize_t i = AU->rowsCount - 1; i >= 0; i--) {
		X->mat[i][0] = AU->mat[i][AU->colsCount-1];
		for (size_t j = i + 1; j < AU->rowsCount; j++) {
			X->mat[i][0] -= AU->mat[i][j] * X->mat[j][0];
		}
		X->mat[i][0] /= AU->mat[i][i];
	}

	freeMat$(AU);

	return X;
}
#pragma endregion "Solving routines"
