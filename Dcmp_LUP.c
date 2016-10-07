#include <math.h>
#include <stdlib.h>

#include "Dcmp_LUP.h"
#include "Matrix.h"
#include "MatrixOps.h"
#include "Gauss.h"
#include "Maths.h"
#include "Extras.h"


#pragma region "Gauss"
/**
 \fn	Mat *Dcmp_LU_Gauss (Mat A)

 \brief	The Gaussian elimination algorithm (with partial (row) pivoting)
		for obtaining LU decomposition of Matrix A.
		NOTE: L and U can be stored in one matrix LU where diagonal
		(containing only 1.0's) of L is omitted.

 \param	A	The Mat to process.

 \return	The * to Matrices array containing L, U &amp; P.
 			\[0] is L, [1] is U, [2] is P. \.
 */
Mat *Dcmp_LU_Gauss (Mat A) {
	Mat LU = DeepCopy(A);
	entry_t **lu = LU->data;
	Mat P = Identity(A->rowsCount);
//	size_t *permutationVector = AllocVec_u(rows);

	//Pivotize and eliminate
	//TODO: one procedure to pivotize them all
	for (size_t k = 0; k < A->rowsCount; k++) {
		//Find pivot.
		size_t piv = k;

		for (size_t i = k + 1; i < A->rowsCount; i++) {
			//Reference implementation of pivoting algorithm
			if (abs(lu[i][k]) > abs(lu[piv][k])) {
				piv = i;
			}
		}

		if (isnotzero(lu[piv][k])) {
			//Swap rows
			if (piv > k) {
				swapRows(LU, piv, k);
				swapRows(P, piv, k);

				//_swap_i(permutationVector[piv], permutationVector[k]);
				P->permutationSign *= -1; //-V127
			}

			//Compute multipliers and eliminate k-th column.
			for (size_t i = k + 1; i < A->rowsCount; i++) {
				entry_t f = lu[i][k] / lu[k][k];

				for (size_t j = k + 1; j < A->colsCount; j++) {
					lu[i][j] -= f * lu[k][j];
				}

				lu[i][k] /= lu[k][k];
			}
		} else {
			A->isSingular = true;
			Check$(!(A->isSingular), "Singular matrix.");
		}
	}

	//Fill L
	Mat L = Identity(A->rowsCount);
	entry_t **l = L->data;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < i; j++) {
			l[i][j] = lu[i][j];
		}
	}

	//Fill U
	Mat U = Zeroes(A->rowsCount);
	entry_t **u = U->data;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = i; j < A->colsCount; j++) {
			u[i][j] = lu[i][j];
		}
	}

	Mat *result = (Mat*) malloc(3 * sizeof(*result));
	Assert$(result != NULL, "Cannot allocate.");
	result[0] = L;
	result[1] = U;
	result[2] = P;

	freeMat$(LU);
	//free$(permutationVector);

	return result;
}
#pragma endregion "Gauss"


#pragma region "Crout"
/**
 \fn	Mat Pivotize_LU (Mat A)

 \brief	Pivotize matrix for further using in LUP decomposition process.

 \param	[in,out] A	The Matrix to process.	Note that A will be modified too.

 \return	Permutation matrix.
 */
Mat Pivotize_LU (Mat A) {
	Mat P = Identity(A->rowsCount);
	entry_t **a = A->data;

	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t piv = k;

		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (abs(a[i][k]) > abs(a[piv][k])) {
				piv = i;
			}
		}

		if (piv > k) {
			swapRows(A, piv, k);
			swapRows(P, piv, k);

			P->permutationSign *= -1; //-V127
		}
	}

	return P;
}

/**
 \fn	Mat *Dcmp_LU_Crout (Mat A)

 \brief	LUP decomposition using Crout's method with partial pivoting, where L is lower triangular
		(has elements only on the diagonal and below) and U is upper triangular (has elements
		only on the diagonal and above). P is the permutation matrix of A produced by partial
		(row) pivoting method. Non-pivoted matrices can lead this algorithm to numerical
		instability (division by 0 or such shit).
		NOTE: L and U can be stored in one matrix LU
		where diagonal (containing only 1.0's) of L is omitted.

 \date	05-Jun-14

 \param	A	The Matrix to factorize.

 \return	The * to Matrices array containing L, U &amp; P. \[0] is L, [1] is U, [2] is P. \.
 */
Mat *Dcmp_LU_Crout (Mat A) {
	Mat L = Identity(A->rowsCount);
	Mat U = AllocMat(A->rowsCount, A->colsCount);
	fill_zeroes(U);
	Mat A_copy = DeepCopy(A);
	Mat P = Pivotize_LU(A_copy);

	entry_t **l = L->data;
	entry_t **u = U->data;
	entry_t **a = A_copy->data;

	for (size_t j = 0; j < A->rowsCount; j++) {
		for (size_t i = 0; i <= j; i++) {
			entry_t sum = 0.0;

			for (size_t k = 0; k < i; k++) {
				sum += u[k][j] * l[i][k];
			}

			u[i][j] = a[i][j] - sum;
		}

		for (size_t i = j; i < A->rowsCount; i++) {
			entry_t sum = 0.0;

			for (size_t k = 0; k < j; k++) {
				sum += u[k][j] * l[i][k];
			}

			l[i][j] = (a[i][j] - sum) / u[j][j];
		}
	}

	Mat *result = (Mat*) malloc(3 * sizeof(*result));
	Assert$(result != NULL, "Cannot allocate.");
	result[0] = L;
	result[1] = U;
	result[2] = P;

	freeMat$(A_copy);

	return result;
}
#pragma endregion "Crout"


/**
 \fn	Mat Solve_LUP (Mat *lup, Mat B)

 \brief	Solves system of linear equations using LUP decomposition. System can be solved
		directly by forward and backward substitution without using the Gaussian elimination
		process (however we do need this process or equivalent to compute the LU decomposition
		itself).

 \param [in] lup	* to Matrices array containing L, U &amp; P.
 \param	B		  	Right hand side.

 \return			Solution as column-vector.
 */
Mat Solve_LUP (Mat *LUP, Mat B) {
	Check$(!IsSingular_LUP(LUP), "Cannot solve for singular matrix.");
	Assert$(LUP[0]->rowsCount == B->rowsCount, "Rows count mismatch.");

	entry_t **l = LUP[0]->data;
	entry_t **u = LUP[1]->data;
	Mat PB = MatMul$(LUP[2], B);
	entry_t **pb = PB->data;

	Mat Y = AllocMat(B->rowsCount, B->colsCount);
	entry_t **y = Y->data;
	Mat X = AllocMat(B->rowsCount, B->colsCount);
	entry_t **x = X->data;

	for (size_t c = 0; c < B->colsCount; c++) {
		//forward solve L.y = b
		for (size_t i = 0; i < LUP[0]->rowsCount; i++) {
			y[i][c] = pb[i][c];

			for (size_t j = 0; j < i; j++) {
				y[i][c] -= l[i][j] * y[j][c];
			}

			//y[i][c] /= l[i][i];
		}

		//backward solve U.x = y
		for (ptrdiff_t i = LUP[1]->rowsCount - 1; i >= 0; i--) {
			x[i][c] = y[i][c];

			for (size_t j = (size_t) (i + 1); j < LUP[1]->colsCount; j++) {
				x[i][c] -= u[i][j] * x[j][c];
			}

			x[i][c] /= u[i][i];
		}
	}

	freeMats(Y, PB, NULL);

	return X;
}

/**
 \fn	double Det_LUP (Mat *lup)

 \brief	Calculates matrix determinant using LU decomposition.

 \date	12-Jun-14

 \param [in] lup	* to Matrices array containing L, U &amp; P.

 \return			Determinant value.
 */
double Det_LUP (Mat *LUP) {
	entry_t **u = LUP[1]->data;
	entry_t det = LUP[2]->permutationSign;

	for (size_t i = 0; i < LUP[1]->rowsCount; i++) {
		det *= u[i][i];
	}

	return det;
}

/**
 \fn	bool IsSingular_LUP (Mat *lup)

 \brief	Checks if matrix is singular.

 \param [in] LUP	* to Matrices array containing L, U &amp; P.

 \return			`true` if Matrix A singular, `false` if it is not.
 */
bool IsSingular_LUP (Mat *LUP) {
	for (size_t i = 0; i < LUP[1]->rowsCount; i++) {
		if (iszero(LUP[1]->data[i][i])) {
			return true;
		}
	}

	return false;
}
