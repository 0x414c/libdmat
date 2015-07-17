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
	entry_t **lu = LU->mat;
	Mat P = Identity(A->rowsCount);
	entry_t **p = P->mat;
	size_t cols = A->colsCount;
	size_t rows = A->rowsCount;
//	size_t *permutationVector = AllocVec_u(rows);
	int8_t permutationSign = 1;

	// Pivoting
	// TODO: one procedure to pivotize them all
	for (size_t k = 0; k < cols; k++) {
		// Find pivot.
		size_t pivot = k;
		for (size_t i = k + 1; i < rows; i++) {
			if (abs(lu[i][k]) > abs(lu[pivot][k])) {
				pivot = i;
			}
		}
		// Swap rows
		if (pivot != k) {
			for (size_t j = 0; j < cols; j++) {
				swap(lu[pivot][j], lu[k][j]);
				swap(p[pivot][j], p[k][j]);
			}
			//__swap_i(permutationVector[pivot], permutationVector[k]);
			permutationSign *= -1; //-V127
		}
		// Compute multipliers and eliminate k-th column.
		if (isnotzero(lu[k][k])) {
			for (size_t i = k + 1; i < rows; i++) {
				lu[i][k] /= lu[k][k];
				for (size_t j = k + 1; j < cols; j++) {
					lu[i][j] -= lu[i][k] * lu[k][j];
				}
			}
		}
	}
	P->permutationSign = permutationSign;

	// fill L
	Mat L = Identity(rows);
	entry_t **l = L->mat;
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < i; j++) {
			l[i][j] = lu[i][j];
		}
	}

	// fill U
	Mat U = AllocMat(rows, cols);
	fill_zeroes(U);
	entry_t **u = U->mat;
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = i; j < cols; j++) {
			u[i][j] = lu[i][j];
		}
	}

	Mat *result = (Mat*) malloc(3*sizeof(*result));
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

 \return	Pivoting matrix.
 */
Mat Pivotize_LU (Mat A) {
	Mat P = Identity(A->rowsCount);
	entry_t **a = A->mat;
	entry_t **p = P->mat;
	int8_t permutationSign = 1;

	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		for (size_t i = k; i < A->rowsCount; i++) {
			if (abs(a[i][k]) > abs(a[pivot][k])) {
				pivot = i;
			}
		}
		if (pivot != k) {
			permutationSign *= -1; //-V127
			for (size_t j = 0; j < A->colsCount; j++) {
				swap(p[k][j], p[pivot][j]);
				swap(a[k][j], a[pivot][j]);
			}
		}
	}
	P->permutationSign = permutationSign;

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
	Mat A_copy;
	Mat L, U, P;
	size_t n = A->rowsCount;

	L = Identity(A->rowsCount);
	U = AllocMat(A->rowsCount, A->colsCount);
	fill_zeroes(U);
	A_copy = DeepCopy(A);
	P = Pivotize_LU(A_copy);

	entry_t **l = L->mat;
	entry_t **u = U->mat;
	entry_t **a = A_copy->mat;

	for (size_t j = 0; j < n; j++) {
		for (size_t i = 0; i <= j; i++) {
			entry_t sum = 0.0;
			for (size_t k = 0; k < i; k++) {
				sum += u[k][j] * l[i][k];
			}
			u[i][j] = a[i][j] - sum;
		}
		for (size_t i = j; i < n; i++) {
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
	Assert$(LUP[0]->rowsCount == B->rowsCount, "Rows count mismatch.");
	Check$(isSingular_LUP(LUP) == false, "Cannot solve for singular matrix.");

	entry_t **l = LUP[0]->mat;
	entry_t **u = LUP[1]->mat;
	Mat PB = MatMul$(LUP[2], B);
	entry_t **b = PB->mat;

	Mat Y = AllocMat(B->rowsCount, B->colsCount);
	entry_t **y = Y->mat;
	Mat X = AllocMat(B->rowsCount, B->colsCount);
	entry_t **x = X->mat;

	for (size_t c = 0; c < B->colsCount; c++) {
		// forward solve Ly = b
		for (size_t i = 0; i < LUP[0]->rowsCount; i++) {
			y[i][c] = b[i][c];
			for (size_t j = 0; j < i; j++) {
				y[i][c] -= l[i][j] * y[j][c];
			}
			//y[i][c] /= l[i][i];
		}
		// backward solve Ux=y
		for (ptrdiff_t i = LUP[1]->rowsCount - 1; i >= 0; i--) {
			x[i][c] = y[i][c];
			for (size_t j = i + 1; j < LUP[1]->colsCount; j++) {
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
	entry_t **u = LUP[1]->mat;
	entry_t det = LUP[2]->permutationSign;

	for (size_t i = 0; i < LUP[1]->rowsCount; i++) {
		det *= u[i][i];
	}

	return det;
}

/**
 \fn	bool isSingular_LUP (Mat *lup)

 \brief	Checks if matrix is singular.

 \param [in] LUP	* to Matrices array containing L, U &amp; P.

 \return			`true` if Matrix A singular, `false` if it is not.
 */
bool isSingular_LUP (Mat *LUP) {
	for (size_t i = 0; i < LUP[1]->rowsCount; i++) {
		if (iszero(LUP[1]->mat[i][i])) {
			return true;
		}
	}

	return false;
}
