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
	entry_t **a = A->data;
	entry_t det;

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
				det = (entry_t) T->permutationSign;

				for (size_t i = 0; i < A->rowsCount; i++) {
					det *= T->data[i][i];
				}
			} else {
				det = 0.0;

				A->isSingular = true;
			}

			freeMat$(T);

            A->det = det;

            return det;
	}

	if (iszero(det)) {
		A->isSingular = true;
	}

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
	entry_t **t = T->data;

	// Pivotize
	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t piv = k;

		for (size_t i = k + 1; i < A->rowsCount; i++) {
			// Reference implementation of pivoting algorithm
			if (abs(t[i][k]) > abs(t[piv][k])) {
				piv = i;
			}
		}

		// Swap rows
		if (piv > k) {
			swapRows(T, piv, k);

			T->permutationSign *= -1; //-V127
		}
	}

	// Bareiss algorithm main step
	for (size_t i = 0; i < T->rowsCount - 1; i++) {
		if (isnotzero(t[i][i])) {
			for (size_t j = i + 1; j < T->rowsCount; j++) {
				for (size_t k = i + 1; k < T->rowsCount; k++) {
					t[j][k] = t[j][k] * t[i][i] - t[j][i] * t[i][k];

					if (i != 0) {
						t[j][k] /= t[i - 1][i - 1];
					}
				}
			}
		} else {
			A->isSingular = true;
			Check$(!(A->isSingular), "Singular matrix.");

			freeMat$(T);

			A->det = 0.0;

			return 0.0;
		}
	}

	entry_t det = t[T->rowsCount - 1][T->rowsCount - 1] * T->permutationSign;

	freeMat$(T);

	A->det = det;

	return det;
}
#pragma endregion "Determinant computation"


#pragma region "Transforming routines"
/**
 \fn		void toRowEchelonForm (Mat A)

 \brief		Transforms matrix A into a row echelon form.
            Optimized implementation.

 \param	A	The Matrix to process.
 */
void toRowEchelonForm (Mat A) {
	entry_t **a = A->data;

	for (size_t k = 0; k < A->rowsCount; k++) {
		// Find pivot.
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			// Immediate rows swapping
			if (abs(a[i][k]) > abs(a[k][k])) {
				swapRows(A, i, k);

				A->permutationSign *= -1; //-V127
			}
		}

		// Compute multipliers and eliminate k-th column.
		if (isnotzero(a[k][k])) {
			for (size_t i = k + 1; i < A->rowsCount; i++) {
				entry_t f = a[i][k] / a[k][k];

				for (size_t j = k; j < A->colsCount; j++) {
					a[i][j] -= f * a[k][j];
				}
			}
		} else {
			A->isSingular = true;
			Check$(!(A->isSingular), "Singular matrix.");
		}
	}

	return;
}

/**
 \fn		void toRowEchelonForm_reference (Mat A)

 \brief		Transforms Matrix A into a row echelon form.
 			Reference implementation w/ partial pivoting.

 \param	A	The Matrix to process.
 */
void toRowEchelonForm_reference (Mat A) {
	entry_t **a = A->data;

	for (size_t k = 0; k < A->rowsCount; k++) {
		// Find pivot
		size_t piv = k;

		for (size_t i = k + 1; i < A->rowsCount; i++) {
			// Reference implementation of pivoting algorithm
			if (abs(a[i][k]) > abs(a[piv][k])) {
				piv = i;
			}
		}

		if (isnotzero(a[piv][k])) {
			// Swap rows
			if (piv > k) {
				swapRows(A, piv, k);

				A->permutationSign *= -1; //-V127
			}

			// Compute multipliers and eliminate k-th column.
			for (size_t i = k + 1; i < A->rowsCount; i++) {
				entry_t f = a[i][k] / a[k][k];

				for (size_t j = k + 1; j < A->colsCount; j++) {
					a[i][j] -= f * a[k][j];
				}

				a[i][k] = 0.0;
			}
		} else {
			A->isSingular = true;
			Check$(!(A->isSingular), "Singular matrix.");
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
	entry_t **a = A->data;
	size_t piv_col = 0;

	for (size_t k = 0; k < A->rowsCount; k++) {
		if (piv_col == A->colsCount) {
			break;
		} else {
			size_t piv_row = k;

			while (true) {
				for (size_t i = k + 1; i < A->rowsCount; i++) {
					if (abs(a[i][piv_col]) > a[piv_row][piv_col]) {
						piv_row = i;
					}
				}

				if (iszero(a[piv_row][piv_col])) {
					A->isSingular = true;
					Check$(!(A->isSingular), "Singular matrix.");

					piv_col++;

					if (piv_col == A->colsCount) {
						return;
					}
				} else {
					break;
				}
			}

			if (piv_row > k) {
				swapRows(A, piv_row, k);

				A->permutationSign *= -1; //-V127
			}

			if (!equals(a[k][piv_col], (entry_t) 1.0)) {
				entry_t d = a[k][piv_col];

				for (size_t j = 0; j < A->colsCount; j++) {
					a[k][j] /= d;
				}
			}

			for (size_t i = 0; i < k; i++) {
				entry_t f = a[i][piv_col];

				for (size_t j = 0; j < A->colsCount; j++) {
					a[i][j] -= f * a[k][j];
				}
			}

			for (size_t i = k + 1; i < A->rowsCount; i++) {
				entry_t f = a[i][piv_col];

				for (size_t j = 0; j < A->colsCount; j++) {
					a[i][j] -= f * a[k][j];
				}
			}

			piv_col++;
		}
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
	Assert$(B->colsCount == 1, "Matrix B should be column-vector.");
	Assert$(A->rowsCount == B->rowsCount, "Number of equations must be equal to the number of unknown variables.");

	Mat X = AllocMat(A->rowsCount, 1);
	Mat AB = DeepCopy(A);

	join(AB, B);
	entry_t **au = AB->data;
	entry_t **x = X->data;

	toReducedRowEchelonForm(AB);

	Check$(!(AB->isSingular), "Cannot solve for singular matrix.");

	for (size_t i = 0; i < X->rowsCount; i++) {
		x[i][0] = au[i][AB->colsCount - 1];
	}

	freeMat$(AB);

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
	Assert$(B->colsCount == 1, "Matrix B should be column-vector.");
	Assert$(A->rowsCount == B->rowsCount, "Number of equations must be equal to number of unknown variables.");

	Mat X = AllocMat(B->rowsCount, B->colsCount);
	Mat AB = DeepCopy(A);

	join(AB, B);

	//Forward step (elimination with row pivoting)
	toRowEchelonForm_reference(AB);

	Check$(!(AB->isSingular), "Cannot solve for singular matrix.");

	//Back-substitution step
	for (ssize_t i = AB->rowsCount - 1; i >= 0; i--) {
		X->data[i][0] = AB->data[i][AB->colsCount-1];

		for (size_t j = (size_t) (i + 1); j < AB->rowsCount; j++) {
			X->data[i][0] -= AB->data[i][j] * X->data[j][0];
		}

		X->data[i][0] /= AB->data[i][i];
	}

	freeMat$(AB);

	return X;
}
#pragma endregion "Solving routines"
