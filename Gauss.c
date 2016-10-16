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

 \brief		Calculates matrix determinant by Gauss' method for n > 4 and with rule of Sarrus for n <= 4.

 \date		13-May-14

 \param	A	The Mat to process.

 \return	Determinant value.
 */
entry_type Det_Gauss (Mat A) {
	Assert$(IsSquare$(A), "Matrix A should be square.");

	Mat A1 = NULL;
	entry_type **a = A->data;
	entry_type det;

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
			A1 = DeepCopy(A);

			toUpperTriangularForm_ref(A1);

			if (!(A1->isSingular)) {
				det = (entry_type) A1->permutationSign;

				for (size_t i = 0; i < A->rowsCount; i++) {
					det *= A1->data[i][i];
				}

				A->isSingular = false;
			} else {
				det = 0.0;

				A->isSingular = true;
			}

			freeMat$(A1);

			A->det = det;

			return det;
	}

	if (iszero(det)) {
		A->isSingular = true;
	} else {
		A->isSingular = false;
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
entry_type Det_Bareiss (Mat A) {
	Assert$(IsSquare$(A), "Matrix A should be square.");

	Mat A1 = DeepCopy(A); //TODO: replace w/ Copy()
	entry_type **a1 = A1->data;

	//Pivotize
	for (size_t k = 0; k < A->rowsCount; k++) {
		//Find pivot
		size_t piv_row = k;

		for (size_t i = k + 1; i < A->rowsCount; i++) {
			//Reference implementation of pivoting algorithm
			if (abs(a1[i][k]) > abs(a1[piv_row][k])) {
				piv_row = i;
			}
		}

		//Swap rows
		if (piv_row > k) {
			swapRows(A1, piv_row, k);

			A1->permutationSign *= -1; //-V127
		}
	}

	//Bareiss algorithm main step
	for (size_t i = 0; i < A1->rowsCount - 1; i++) {
		if (isnotzero(a1[i][i])) {
			for (size_t k = i + 1; k < A1->rowsCount; k++) {
				for (size_t j = i + 1; j < A1->colsCount; j++) {
					a1[k][j] = a1[k][j] * a1[i][i] - a1[i][j] * a1[k][i];

					if (i > 0) {
						a1[k][j] /= a1[i - 1][i - 1];
					}
				}
			}
		} else {
			A1->isSingular = true;
			Check$(!(A1->isSingular), "Singular matrix.");
			A->isSingular = A1->isSingular;
			freeMat$(A1);
			A->det = 0.0;

			return 0.0;
		}
	}

	entry_type det = a1[A1->rowsCount - 1][A1->rowsCount - 1] * A1->permutationSign;

	freeMat$(A1);

	A->det = det;

	return det;
}
#pragma endregion "Determinant computation"


#pragma region "Transforming routines"
/**
 \fn		void toUpperTriangularForm (Mat A)

 \brief		Transforms matrix A into an upper triangular form.

 \param	A	The Matrix to process.
 */
void toUpperTriangularForm (Mat A) {
	Check$(IsSquare$(A), "Matrix A should be square.");

	entry_type **a = A->data;

	//Pivotize
	for (size_t k = 0; k < A->rowsCount; k++) {
		//Find pivot.
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			//Immediate rows swapping
			if (abs(a[i][k]) > abs(a[k][k])) {
				swapRows(A, i, k);

				A->permutationSign *= -1; //-V127
			}
		}

		//Compute multipliers and eliminate k-th column.
		if (isnotzero(a[k][k])) {
			for (size_t i = k + 1; i < A->rowsCount; i++) {
				entry_type f = a[i][k] / a[k][k];

				for (size_t j = k; j < A->colsCount; j++) {
					a[i][j] -= f * a[k][j];
				}
			}
		} else {
			A->isSingular = true;
		}
	}

	Check$(!(A->isSingular), "Singular matrix.");

	return;
}

/**
 \fn		void toUpperTriangularForm_ref (Mat A)

 \brief		Transforms Matrix A into an upper triangular form.
 			Reference implementation w/ partial pivoting.

 \param	A	The Matrix to process.
 */
void toUpperTriangularForm_ref (Mat A) {
	Check$(IsSquare$(A), "Matrix A should be square.");

	entry_type **a = A->data;

	//Pivotize and eliminate
	for (size_t k = 0; k < A->rowsCount; k++) {
		//Find pivot
		size_t piv_row = k;

		for (size_t i = k + 1; i < A->rowsCount; i++) {
			//Reference implementation of pivoting algorithm
			if (abs(a[i][k]) > abs(a[piv_row][k])) {
				piv_row = i;
			}
		}

		if (isnotzero(a[piv_row][k])) {
			//Swap rows
			if (piv_row > k) {
				swapRows(A, piv_row, k);

				A->permutationSign *= -1; //-V127
			}

			//Compute multipliers and eliminate k-th column
			for (size_t i = k + 1; i < A->rowsCount; i++) {
				entry_type f = a[i][k] / a[k][k];

				for (size_t j = k + 1; j < A->colsCount; j++) {
					a[i][j] -= f * a[k][j];
				}

				a[i][k] = 0.0;
			}
		} else {
			A->isSingular = true;
		}
	}

	Check$(!(A->isSingular), "Singular matrix.");

	return;
}

/**
 \fn		void toRowEchelonForm (Mat A)

 \brief		Transforms Matrix A into a row echelon form.
 			Reference implementation w/ partial pivoting.

 \param	A	The Matrix to process.
 */
void toRowEchelonForm (Mat A) {
	entry_type **a = A->data;
	size_t piv_col = 0;

	//Pivotize and eliminate
	for (size_t k = 0; k < A->rowsCount; k++) {
		if (piv_col == A->colsCount) {
			break;
		} else {
			//Find pivot
			size_t piv_row = k;

			while (true) {
				for (size_t i = k + 1; i < A->rowsCount; i++) {
					//Reference implementation of pivoting algorithm
					if (abs(a[i][piv_col]) > a[piv_row][piv_col]) {
						piv_row = i;
					}
				}

				if (iszero(a[piv_row][piv_col])) {
					A->isSingular = true;

					piv_col++;

					if (piv_col == A->colsCount) {
						goto ref_end;
					}
				} else {
					break;
				}
			}

			//Swap rows
			if (piv_row > k) {
				swapRows(A, piv_row, k);

				A->permutationSign *= -1; //-V127
			}

			//Compute multipliers and eliminate piv_col-th column
			for (size_t i = k + 1; i < A->rowsCount; i++) {
				entry_type f = a[i][piv_col] / a[k][piv_col];

				for (size_t j = piv_col + 1; j < A->colsCount; j++) {
					a[i][j] -= f * a[k][j];
				}

				a[i][piv_col] = 0.0;
			}

			piv_col++;
		}
	}

ref_end:
	Check$(!(A->isSingular), "Singular matrix.");

	return;
}

void toRowEchelonForm_Bareiss (Mat A) {
	entry_type **a = A->data;

	//Pivotize
	for (size_t k = 0; k < A->rowsCount; k++) {
		//Find pivot
		size_t piv_row = k;

		for (size_t i = k + 1; i < A->rowsCount; i++) {
			//Reference implementation of pivoting algorithm
			if (abs(a[i][k]) > abs(a[piv_row][k])) {
				piv_row = i;
			}
		}

		//Swap rows
		if (piv_row > k) {
			swapRows(A, piv_row, k);

			A->permutationSign *= -1; //-V127
		}
	}

	entry_type prev_piv = 1.0;

	for (size_t i = 0; i < A->rowsCount; i++) {
		if (isnotzero(a[i][i])) {
			entry_type piv = a[i][i];

			for (size_t k = i + 1; k < A->rowsCount; k++) {
				for (size_t j = i + 1; j < A->colsCount; j++) {
					a[k][j] = (a[k][j] * piv - a[i][j] * a[k][i]) / prev_piv;
				}

				a[k][i] = 0.0;
			}

			prev_piv = piv;
		} else {
			A->isSingular = true;
		}
	}

	Check$(!(A->isSingular), "Singular matrix.");
}

void toReducedRowEchelonForm_Bareiss (Mat A) {
	entry_type **a = A->data;

	//Pivotize
	for (size_t k = 0; k < A->rowsCount; k++) {
		//Find pivot
		size_t piv_row = k;

		for (size_t i = k + 1; i < A->rowsCount; i++) {
			//Reference implementation of pivoting algorithm
			if (abs(a[i][k]) > abs(a[piv_row][k])) {
				piv_row = i;
			}
		}

		//Swap rows
		if (piv_row > k) {
			swapRows(A, piv_row, k);

			A->permutationSign *= -1; //-V127
		}
	}

	entry_type prev_piv = 1.0;

	for (size_t i = 0; i < A->rowsCount; i++) {
		if (isnotzero(a[i][i])) {
			entry_type piv = a[i][i];

			for (size_t k = 0; k < A->rowsCount; k++) {
				if (k != i) {
					for (size_t j = 0; j < A->colsCount; j++) {
						if (j != i) {
							a[k][j] = (a[k][j] * piv - a[i][j] * a[k][i]) / prev_piv;
						}
					}

					a[k][i] = 0.0;
				}
			}

			prev_piv = piv;
		} else {
			A->isSingular = true;
		}
	}

	Check$(!(A->isSingular), "Singular matrix.");
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
	entry_type **a = A->data;
	size_t piv_col = 0;

	//Pivotize and eliminate
	for (size_t k = 0; k < A->rowsCount; k++) {
		if (piv_col == A->colsCount) {
			break;
		} else {
			//Find pivot
			size_t piv_row = k;

			while (true) {
				for (size_t i = k + 1; i < A->rowsCount; i++) {
					//Reference implementation of pivoting algorithm
					if (abs(a[i][piv_col]) > a[piv_row][piv_col]) {
						piv_row = i;
					}
				}

				if (iszero(a[piv_row][piv_col])) {
					A->isSingular = true;

					piv_col++;

					if (piv_col == A->colsCount) {
						goto rref_end;
					}
				} else {
					break;
				}
			}

			//Swap rows
			if (piv_row > k) {
				swapRows(A, piv_row, k);

				A->permutationSign *= -1; //-V127
			}

			if (!equals(a[k][piv_col], (entry_type) 1.0)) {
				entry_type d = a[k][piv_col];

				for (size_t j = 0; j < A->colsCount; j++) {
					a[k][j] /= d;
				}
			}

			for (size_t i = 0; i < A->rowsCount; i++) {
				if (i != k) {
					entry_type f = a[i][piv_col];

					for (size_t j = 0; j < A->colsCount; j++) {
						a[i][j] -= f * a[k][j];
					}
				}
			}

			piv_col++;
		}
	}

rref_end:
	Check$(!(A->isSingular), "Singular matrix.");

	return;
}
#pragma endregion "Transforming routines"


#pragma region "Solving routines"
/**
 \fn		Mat Solve_Gauss (Mat A, Mat B)

 \brief		Solves system of linear equations using Gauss elimination.

 \param	A	Coeffs matrix.
 \param	B	Right hand side as column-vector.

 \return	Solution as column-vector.
 */
Mat Solve_Gauss (Mat A, Mat B) {
	Assert$(A->rowsCount == B->rowsCount, "Number of equations must be equal to number of unknown variables.");
	Assert$(B->colsCount == 1, "Matrix B should be column-vector.");

	Mat AB = DeepCopy(A);

	joinColumns(AB, B);

	//Forward step (elimination with row pivoting)
	toRowEchelonForm(AB);

	Check$(!(AB->isSingular), "Cannot solve for singular matrix.");

	Mat X = AllocMat(B->rowsCount, 1);

	entry_type **ab = AB->data;
	entry_type **x = X->data;

	//Back-substitution step
	for (ssize_t i = AB->rowsCount - 1; i >= 0; i--) {
		x[i][0] = ab[i][AB->colsCount - 1];

		for (size_t j = (size_t) (i + 1); j < AB->rowsCount; j++) {
			x[i][0] -= ab[i][j] * x[j][0];
		}

		x[i][0] /= ab[i][i];
	}

	freeMat$(AB);

	return X;
}

/**
 \fn		Mat Solve_GaussJordan (Mat A, Mat B)

 \brief		Solves system of linear equations using Gauss-Jordan method.

 \param	A	Coefficients matrix.
 \param	B	Right hand side as column-vector.

 \return	Solution as column-vector.
 */
Mat Solve_GaussJordan (Mat A, Mat B) {
	Assert$(A->rowsCount == B->rowsCount, "Number of equations must be equal to the number of unknown variables.");

	Mat AB = DeepCopy(A);

	joinColumns(AB, B);

	toReducedRowEchelonForm(AB);

	Check$(!(AB->isSingular), "Cannot solve for singular matrix.");

	Mat X = AllocMat(B->rowsCount, B->colsCount);

	entry_type **ab = AB->data;
	entry_type **x = X->data;

	for (size_t i = 0; i < X->rowsCount; i++) {
		for (size_t j = A->colsCount; j < AB->colsCount; j++) {
			x[i][j - A->colsCount] = ab[i][j];
		}
	}

	freeMat$(AB);

	return X;
}

/**
 \fn		Mat Solve_Bareiss (Mat A, Mat B)

 \brief		Solves system of linear equations using Bareiss algorithm.

 \param	A	Coeffs matrix.
 \param	B	Right hand side as column-vector.

 \return	Solution as column-vector.
 */
Mat Solve_Bareiss (Mat A, Mat B) {
	Assert$(A->rowsCount == B->rowsCount, "Number of equations must be equal to number of unknown variables.");
	Assert$(B->colsCount == 1, "Matrix B should be column-vector.");

	Mat AB = DeepCopy(A);

	joinColumns(AB, B);

	//Forward step (elimination with row pivoting)
	toRowEchelonForm_Bareiss(AB);

	Check$(!(AB->isSingular), "Cannot solve for singular matrix.");

	Mat X = AllocMat(B->rowsCount, 1);

	entry_type **ab = AB->data;
	entry_type **x = X->data;

	//Back-substitution step
	for (ssize_t i = AB->rowsCount - 1; i >= 0; i--) {
		x[i][0] = ab[i][AB->colsCount - 1];

		for (size_t j = (size_t) (i + 1); j < AB->rowsCount; j++) {
			x[i][0] -= ab[i][j] * x[j][0];
		}

		x[i][0] /= ab[i][i];
	}

	freeMat$(AB);

	return X;
}

/**
 \fn		Mat Solve_Montante (Mat A, Mat B)

 \brief		Solves system of linear equations using Montante method (based on Bareiss algorithm).

 \param	A	Coefficients matrix.
 \param	B	Right hand side as column-vector.

 \return	Solution as column-vector.
 */
Mat Solve_Montante (Mat A, Mat B) {
	Assert$(A->rowsCount == B->rowsCount, "Number of equations must be equal to the number of unknown variables.");

	Mat AB = DeepCopy(A);

	joinColumns(AB, B);

	toReducedRowEchelonForm_Bareiss(AB);

	Check$(!(AB->isSingular), "Cannot solve for singular matrix.");

	Mat X = AllocMat(B->rowsCount, B->colsCount);

	entry_type **ab = AB->data;
	entry_type **x = X->data;

	for (size_t i = 0; i < X->rowsCount; i++) {
		for (size_t j = A->colsCount; j < AB->colsCount; j++) {
			x[i][j - A->colsCount] = ab[i][j] / ab[i][i];
		}
	}

	freeMat$(AB);

	return X;
}
#pragma endregion "Solving routines"
