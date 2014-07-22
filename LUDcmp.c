#include <math.h>
#include <stdlib.h>

#include "LUDcmp.h"
#include "Matrix.h"
#include "MatrixOperations.h"
#include "Gauss.h"
#include "Extra.h"


Mat *LUDcmp_gauss (Mat A) {
	Mat LU = DeepCopy(A);
	double **lu = LU->a;
	Mat P = Identity(A->rowsCount);
	//double **p = P->a;
	size_t cols = A->colsCount;
	size_t rows = A->rowsCount;
	//size_t *permutationVector = uAllocVec(m);
	int permutationSign = 1;

	// Pivoting //TODO: generalize? to single procedure
	for (size_t k = 0; k < cols; k++) {
		// Find pivot.
		size_t pivot = k;
		for (size_t i = k + 1; i < rows; i++) {
			if (fabs(lu[i][k]) > fabs(lu[pivot][k])) {
				pivot = i;
			}
		}
		// Swap rows
		if (pivot != k) {
			for (size_t j = 0; j < cols; j++) {
				swap_d(lu[pivot][j], lu[k][j]);
				swap_d(P->a[pivot][j], P->a[k][j]);
			}
			//swap_i(permutationVector[pivot], permutationVector[k]);
			permutationSign *= -1; //-V127
		}
		// Compute multipliers and eliminate k-th column.
		if (fabs(lu[k][k]) > EPS) {
			for (size_t i = k + 1; i < rows; i++) {
				lu[i][k] /= lu[k][k];
				for (size_t j = k + 1; j < cols; j++) {
					lu[i][j] -= lu[i][k] * lu[k][j];
				}
			}
		}
	}
	P->permutationSign = permutationSign;

	Mat L = AllocMat(rows, cols);
	double **l = L->a;
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			if (i > j) {
				l[i][j] = lu[i][j];
			} else {
				l[i][j] = (double) (i == j);
			}
		}
	}

	Mat U = AllocMat(rows, cols);
	double **u = U->a;
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			if (i <= j) {
				u[i][j] = lu[i][j];
			}
		}
	}

	Mat *result = (Mat*) malloc(3*sizeof(*result));
	Assert(result != NULL, "Cannot allocate.");
	result[0] = L;
	result[1] = U;
	result[2] = P;

	FreeMat(LU);
	//free(permutationVector);	

	return result;
}

Mat LUPivotize (Mat A) {
	Mat P = Identity(A->rowsCount);
	double **a = A->a;
	double **p = P->a;
	int permutationSign = 1;

	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		for (size_t i = k; i < A->rowsCount; i++) {
			if (fabs(a[i][k]) > fabs(a[pivot][k])) {
				pivot = i;
			}
		}
		if (pivot != k) {
			permutationSign *= -1; //-V127
			for (size_t j = 0; j < A->colsCount; j++) {
				swap_d(p[k][j], p[pivot][j]);
				swap_d(A->a[k][j], A->a[pivot][j]);
			}
		}
	}
	P->permutationSign = permutationSign;

	return P;
}

/**
 \fn	void LUPDcmp (dMat A)
 \brief	LUP decomposition using Crout's method with partial pivoting,
 where L is lower triangular (has elements only on the diagonal and below) and U
 is upper triangular (has elements only on the diagonal and above).
 P is the permutation matrix of A produced by partial pivoting method.
 Non-pivoted matrices can lead this algorithm to numerical instability (division by 0 or such shit).
 TODO: L and U can be stored in one matrix LU where diagonal of L is omitted.	 
 \date	05-Jun-14																 
 \param	A	The dMat to process.
 */
Mat *LUDcmp_crout (Mat A) {
	Mat A_copy;
	Mat L, U, P;
	size_t n = A->rowsCount;

	L = Identity(A->rowsCount);
	U = AllocMat(A->rowsCount, A->colsCount);
	A_copy = DeepCopy(A);
	Assert(A_copy != NULL, "Cannot create copy...");

	P = LUPivotize(A_copy);

	double **l = L->a;
	double **u = U->a;
	double **a = A_copy->a;
	
	for (size_t j = 0; j < n; j++) {
		for (size_t i = 0; i <= j; i++) {
			long double sum = 0.0;
			for (size_t k = 0; k < i; k++) {
				sum += u[k][j] * l[i][k];
			}
			u[i][j] = a[i][j] - sum;
		}
		for (size_t i = j; i < n; i++) {
			double sum = 0.0;
			for (size_t k = 0; k < j; k++) {
				sum += u[k][j] * l[i][k];
			}
			l[i][j] = (a[i][j] - sum) / u[j][j];
		}
	}

	Mat *result = (Mat*) malloc(3 * sizeof(*result));
	Assert(result != NULL, "Cannot allocate.");
	result[0] = L;
	result[1] = U;
	result[2] = P;

	FreeMat(A_copy);
	
	return result;
}

// Ly=b	forward
// Ux=y	backward
// More precision loss than in QRSolve() 
Mat Solve_lup (Mat *lup, Mat B) {
	Assert(lup[0]->rowsCount == B->rowsCount, "Rows count mismatch.");
	Check(isSingular_lup(lup) == false, "Cannot solve for singular matrix.");

	double **l = lup[0]->a;
	double **u = lup[1]->a;
	Mat PB = MatMul(lup[2], B);
	double **b = PB->a;

	Mat Y = AllocMat(B->rowsCount, B->colsCount);
	double **y = Y->a;
	Mat X = AllocMat(B->rowsCount, B->colsCount);
	double **x = X->a;

	for (size_t c = 0; c < B->colsCount; c++) {
		for (size_t i = 0; i < lup[0]->rowsCount; i++) {
			y[i][c] = b[i][c];
			for (size_t j = 0; j < i; j++) {
				y[i][c] -= l[i][j] * y[j][c];
			}
			//y[i][c] /= l[i][i];
		}
		for (ptrdiff_t i = lup[1]->rowsCount - 1; i >= 0; i--) {
			x[i][c] = y[i][c];
			for (size_t j = i + 1; j < lup[1]->colsCount; j++) {
				x[i][c] -= u[i][j] * x[j][c];
			}
			x[i][c] /= u[i][i];
		}
	}

	freeMats(Y, PB, NULL);

	return X;
}

/**
 \fn	double LUDet (LUP lup)
 \brief	Calculate matrix determinant using LU decomposition.
 \date	12-Jun-14											
 \param	lup	The * to LUP struct to get data from.				
 \return	A double.
 */
double Det_lup (Mat *lup) {
	double** u = lup[1]->a;
	double det = lup[2]->permutationSign;

	for (size_t i = 0; i < lup[1]->rowsCount; i++) {
		det *= u[i][i];
	}

	return det;
}

bool isSingular_lup (Mat *lup) {
	for (size_t i = 0; i < lup[1]->rowsCount; i++) {
		if (fabs(lup[1]->a[i][i]) <= EPS) {
			return true;
		}
	}

	return false;
} 
