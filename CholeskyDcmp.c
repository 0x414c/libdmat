#include <stdbool.h>
#include <math.h>

#include "CholeskyDcmp.h"
#include "Matrix.h"
#include "Const.h"
#include "Extra.h"


/**
 \fn	Mat CholeskyDcmp (Mat A)

 \brief	Cholesky decomposition of symmetric positive-definite matrix. The Cholesky
		decomposition or Cholesky factorization is a decomposition of a positive-definite
		matrix into the product of a lower triangular matrix and its transpose, useful for
		efficient numerical solutions. It was discovered by Andre-Louis Cholesky for real
		matrices. When it is applicable, the Cholesky decomposition is roughly twice as
		efficient as the LU decomposition for solving systems of linear equations. The Cholesky
		decomposition is unique when A is positive definite.

 \date	23-Jul-14

 \param	A	The Mat to process.

 \return	A Lower Triangular form of A.
 */
Mat CholeskyDcmp (Mat A) {
	Assert(isSquare(A), "Matrix is not square.");

	Mat L = AllocMat(A->rowsCount, A->colsCount);
	bool isSPD = true;

	for (size_t j = 0; j < A->rowsCount; j++) {
		double d = 0.0;
		for (size_t k = 0; k < j; k++) {
			double s = 0.0;
			for (size_t i = 0; i < k; i++) {
				s += L->a[k][i] * L->a[j][i];
			}
			L->a[j][k] = s = (A->a[j][k] - s) / L->a[k][k];
			d += s*s;
			isSPD = isSPD && (equal_d(A->a[k][j], A->a[j][k]));
		}
		d = A->a[j][j] - d;
		isSPD = isSPD && (d > EPS);
		L->a[j][j] = sqrt(max(d, 0.0));
		for (size_t k = j + 1; k < A->rowsCount; k++) {
			L->a[j][k] = 0.0;
		}
	}
	L->isSPD = isSPD;

	return L;
}

/**
 \fn	Mat Solve_cholesky (Mat L, Mat B)

 \brief	Solves system of linear equations using Cholesky decomposition. The Cholesky
		decomposition is mainly used for the numerical solution of linear equations Ax = b. If A
		is symmetric and positive definite, then we can solve Ax = b by first computing the
		Cholesky decomposition A = LL*, then solving Ly = b for y by forward substitution, and
		finally solving L*x = y for x by back substitution.

 \param	L	The Matrix processed by CholeskyDcmp (Mat A).
 \param	B	The Matrix containing right hand side to solve over.

 \return	Matrix contatining solution as column-vector.
 */
Mat Solve_cholesky (Mat L, Mat B) {
	Assert(L->isSPD, "Matrix is not symmetric positive definite.");
	Assert(L->rowsCount == B->rowsCount, "Rows count mismatch.");

	Mat X = DeepCopy(B);

	// Solve L*Y = B;
	for (size_t k = 0; k < L->rowsCount; k++) {
		for (size_t j = 0; j < B->colsCount; j++) {
			for (size_t i = 0; i < k; i++) {
				X->a[k][j] -= X->a[i][j] * L->a[k][i];
			}
			X->a[k][j] /= L->a[k][k];
		}
	}

	// Solve L'*X = Y;
	for (ptrdiff_t k = L->rowsCount - 1; k >= 0; k--) {
		for (size_t j = 0; j < B->colsCount; j++) {
			for (size_t i = k + 1; i < L->rowsCount; i++) {
				X->a[k][j] -= X->a[i][j] * L->a[i][k];
			}
			X->a[k][j] /= L->a[k][k];
		}
	}

	return X;
}

/**
 \fn	double Det_cholesky (Mat L)

 \brief	Computes matrix determinant using Cholesky dcmp.

 \param	L	The Lower-triangular Matrix processed by CholeskyDcmp (Mat A).

 \return	Determinant value.
 */
double Det_cholesky (Mat L) {
	double det = 1.0;
	for (size_t i = 0; i < L->rowsCount; i++) {
		det *= L->a[i][i];
	}

	return det;
}
