#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <limits.h>

#include "Gauss.h"
#include "Matrix.h"
#include "MatrixOps.h"
#include "Maths.h"
#include "Extras.h"


#pragma region "Determinant computation"

/**
 \fn		double Det_Gauss (Mat A)

 \brief		Calculates matrix determinant by Gauss' method for n&gt;3 and with rule of Sarrus for n&lt;=3.

 \date		13-May-14

 \param	A	The Mat to process.

 \return	Determinant value.
 */
entry_t Det_Gauss (Mat A) {
	Mat T = NULL;
	entry_t det = 1.0;
	entry_t **a = A->a;

	switch (A->rowsCount) {
		case 1:
			det = (a[0][0]);
			break;
		case 2:
			det = (((a[0][0] * a[1][1]) - (a[0][1] * a[1][0])));
			break;
		case 3:
			det = (
				(a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])) -
				(a[0][1] * (a[2][2] * a[1][0] - a[1][2] * a[2][0])) +
				(a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
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

			toRowEchelonForm(T);

			if (!(T->isSingular)) {
				for (size_t i = 0; i < A->rowsCount; i++) {
					det *= T->a[i][i];
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

	if (fabs(det) <= EPS) { A->isSingular = true; }
	A->det = det;

	return det;
}

/**
 \fn		double Det_Bareiss (Mat A)

 \brief		Computes Matrix determinant by Bareiss' algorithm.
			The Bareiss Algorithm is a fraction-free method for determinant computation.
			However, it can also be thought of as a sophisticated form of row reduction.
			Note that the divisions computed at any step are exact; thus Bareiss’ Algorithm is
			indeed fraction-free. Entry a[n][n] is the determinant of A (after `pivoting` and `main` steps).

 \param	A	The Mat to process.

 \return	Determinant of Matrix A.
 */
entry_t Det_Bareiss (Mat A) { //TODO: move to another file
	Mat T = DeepCopy(A);
	entry_t **a = T->a;

	// Pivotize
	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		for (size_t i = k; i < A->rowsCount; i++) {
			if (fabs(a[i][k]) > fabs(a[pivot][k])) {
				pivot = i;
			}
		}
		if (pivot != k) {
//			for (size_t j = 0; j < A->colsCount; j++) {
//				swap_d(T->a[k][j], T->a[pivot][j]);
//			}
			entry_t *a_pivot = a[pivot];
			a[pivot] = a[k];
			a[k] = a_pivot;
			T->permutationSign *= -1; //-V127
		}
	}

	// Bareiss algorithm main step
	for (size_t i = 0; i < T->rowsCount - 1; i++) {
		Check$(fabs(a[i][i]) > EPS, "Singularity...");
		for (size_t j = i + 1; j < T->rowsCount; j++) {
			for (size_t k = i + 1; k < T->rowsCount; k++) {
				a[j][k] = a[j][k] * a[i][i] - a[j][i] * a[i][k];
				if (i != 0) {
					a[j][k] /= a[i-1][i-1];
				}
			}
		}
	}

	entry_t det = a[T->rowsCount - 1][T->rowsCount - 1] * T->permutationSign;
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
	entry_t **a = A->a;
//	size_t _ic = 0;

	// W/ immediate rows swapping
	for (size_t k = 0; k < A->rowsCount; k++) {
		// Pivotize
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (fabs(a[i][k]) > fabs(a[k][k]))
			if (fabs(a[i][k]) > EPS) {
//				for (size_t j = 0; j < A->colsCount; j++) {
//					swap_d(a[k][j], a[i][j]);
//				}
//				_ic++;
				entry_t *a_i = a[i];
				a[i] = a[k];
				a[k] = a_i;
				A->permutationSign *= -1; //-V127
			}
		}
		// Eliminate
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (fabs(a[k][k]) > EPS) {
				entry_t factor = a[i][k] / a[k][k];
				for (size_t j = k; j < A->colsCount; j++) {
					a[i][j] -= factor * a[k][j];
				}
			} else {
				A->isSingular = true;
				Check$(false, "Singular matrix.");
				return;
			}
		}
	}

//	printf("O:%zu\n", _ic);

	return;
}

/**
 \fn		void toRowEchelonForm_reference (Mat A)

 \brief		Transforms Matrix A into a row echelon form.
 			Reference implementation.

 \param	A	The Matrix to process.
 */
void toRowEchelonForm_reference (Mat A) {
	entry_t **a = A->a;
//	size_t _ic = 0;

	// Reference implementation of pivoting algorithm
	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		// Find pivot
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (fabs(a[i][k]) > fabs(a[pivot][k])) {
				pivot = i;
			}
		}
		if (fabs(a[pivot][k]) <= EPS) {
			A->isSingular = true;
			Check$(false, "Singular matrix.");
			return;
		}
		// Swap rows
		if (pivot != k) {
//			for (size_t j = 0; j < A->colsCount; j++) {
//				swap_d(a[k][j], a[pivot][j]);
//			}
//			_ic++;
			entry_t *a_pivot = a[pivot];
			a[pivot] = a[k];
			a[k] = a_pivot;
			A->permutationSign *= -1; //-V127
		}
		// Calculate factor, eliminate k-th column
		if (fabs(a[k][k]) > EPS) {
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

//	printf("R:%zu\n", _ic);

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
	entry_t **a = A->a;
	size_t pivotCol = 0;
	size_t rowsCount = A->rowsCount, colsCount = A->colsCount;
	size_t i, j, k, r;

	for (r = 0; r < rowsCount; r++) {
		if (colsCount < pivotCol + 1) {
			break;
		}
		i = r;
		while (fabs(a[i][pivotCol]) <= EPS) {
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
			swap_d(a[i][j], a[r][j]);
		}
		if (fabs(a[r][pivotCol]) > EPS) {
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
 \param	B	Right hand side.

 \return	Solution as column-vector.
 */
Mat Solve_GaussJordan (Mat A, Mat B) {
	Assert$(B->colsCount == 1, "Matrix B must be column-vector.");
	Assert$(A->rowsCount == B->rowsCount, "Number of equations doesn't equal to number of unknown variables.");

	Mat X = AllocMat(A->rowsCount, 1);
	Mat AU = DeepCopy(A);

	concat(AU, B);
    entry_t **au = AU->a;
    entry_t **x = X->a;

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
 \param	B	Right hand side.

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
		X->a[i][0] = AU->a[i][AU->colsCount-1];
		for (size_t j = i + 1; j < AU->rowsCount; j++) {
			X->a[i][0] -= AU->a[i][j] * X->a[j][0];
		}
		X->a[i][0] /= AU->a[i][i];
	}

	freeMat$(AU);

	return X;
}
#pragma endregion "Solving routines"




//-----------------------------somewhat outdated, but working...----------------
#pragma region old
/**
 \fn	void simpleSolver (double **a, size_t Size, double *x)

 \brief	Simple back-substitution routine for undeterminedSolver.

 \date	22-May-14

 \param [in] 	a		If non-null, the double ** to matrix in rref form.
 \param			Size	The matrix size.
 \param [out]	x		* to array to write the solution to.
 */
void simpleSolver_h (entry_t **a, size_t size, entry_t *x) {
	for (ssize_t i = size - 1; i >= 0; i--) {
		for (size_t j = i + 1; j < size; j++) {
			x[i] -= a[i][j];
		}
		if ((fabs(x[i]) > EPS) && (fabs(a[i][i]) > EPS)) {
			x[i] /= a[i][i];
		}
	}

	return;
}

/**
 \fn	void undeterminedSolver (dMat RREF, dMat A, dMat R)

 \brief	Solves undetermined system of linear equations.

		Gives fundamental system of solution vectors as rows of matrix.
 \date	22-May-14

 \param	[in]	RREF	The matrix A in rref form.
 \param	[in]	A		The matrix A.
 \param	[out]	R   	The matrix with solution vectors.
 */
void undeterminedSolver_h (Mat RREF, Mat A, Mat R) {
    entry_t **rref = RREF->a;
	size_t *f = AllocVec_u(A->rowsCount);
	size_t i, j, c = 0, row = 0;

	for (i = 0; i < RREF->rowsCount; i++) {
		for (j = 0; j < RREF->colsCount; j++) {
			if (fabs(rref[i][j]) > EPS) {
				f[c++] = i;
				break;
			}
		}
	}
	c = 0;

	for (i = 0; i < RREF->rowsCount; i++) {
		if (!(exists_u(f, RREF->rowsCount, 0, i))) {
			f[c++] = i;
		}
	}
	j = RREF->rowsCount-1;

	do {
		Mat Copy = DeepCopy(A);

        entry_t **copy = Copy->a;
		for (i = 0; i < c-1; i++) {
			copy[(f[i])][(f[i])] += 1;
		}
		toReducedRowEchelonForm(Copy);
		simpleSolver_h(copy, Copy->rowsCount, (R->a[row]));
		for (i = RREF->rowsCount - c; i < RREF->rowsCount; i++) {
			*(R->a[row] + i) = 0.0;
		}
		*(R->a[row] + (j--)) = 1.0;
		freeMat$(Copy);
		row++;
	} while (nextPermutation(f, c - 1, RREF->rowsCount - 1));

	free(f); f = NULL;

	return;
}

/**
 \fn		double *GaussianSolve (dMat A)

 \brief	Solves homogeneous system of linear equations defined by NxN matrix.
		Theorem: Every homogeneous system has either exactly one solution or infinitely many solutions.
		If a homogeneous system has more unknowns than equations, then it has infinitely many solutions.

 \date	15-May-14

 \param	[in,out]	A  	The double-valued matrix to process.

 \return	Matrix with solution vectors written as rows.
*/
Mat GaussianSolve_h (Mat A)	{
    entry_t **a = A->a;
	Mat copyA = DeepCopy(A);

	toReducedRowEchelonForm(A);
	A->rank = Rank(A);

	if (A->rank < A->rowsCount) {
		Mat RES = AllocMat(A->rowsCount - A->rank, A->rowsCount);
		undeterminedSolver_h(A, copyA, RES);
		freeMat$(copyA);
		freeMat$(A);

		return RES;
	} else {
		Mat RES = AllocMat(1, A->rowsCount);
		simpleSolver_h(a, A->rowsCount, RES->a[0]);
		freeMat$(A);

		return RES;
	}
}
#pragma endregion old
