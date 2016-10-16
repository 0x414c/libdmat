#include <stdbool.h>
#include <math.h>

#include "Dcmp_Cholesky.h"
#include "Extras.h"
#include "Maths.h"


/**
 \fn		Mat Dcmp_Cholesky (Mat A)

 \brief		Cholesky decomposition of symmetric positive-definite matrix. The Cholesky
			decomposition or Cholesky factorization is a decomposition of a positive-definite
			matrix into the product of a lower triangular matrix and its transpose, useful for
			efficient numerical solutions. It was discovered by Andre-Louis Cholesky for real
			matrices. When it is applicable, the Cholesky decomposition is roughly twice as
			efficient as the LU decomposition for solving systems of linear equations. The Cholesky
			decomposition is unique when A is positive definite.

 \date		23-Jul-14

 \param	A	The Matrix to process.

 \return	Matrix A in a Lower Triangular form.
 */
Mat Dcmp_Cholesky (Mat A) {
	Assert$(IsSquare$(A), "Matrix A must be square.");

	Mat L = AllocMat(A->rowsCount, A->colsCount);
	bool isSPD = true;

	for (size_t j = 0; j < A->rowsCount; j++) {
		entry_type d = 0.0;
		for (size_t k = 0; k < j; k++) {
			entry_type s = 0.0;

			for (size_t i = 0; i < k; i++) {
				s += L->data[k][i] * L->data[j][i];
			}

			L->data[j][k] = s = (A->data[j][k] - s) / L->data[k][k];
			d += s * s;

			isSPD = isSPD && (equals(A->data[k][j], A->data[j][k]));
		}

		d = A->data[j][j] - d;
		isSPD = isSPD && isnotzero(d);
		L->data[j][j] = sqrt(max(d, (entry_type) 0.0));

		for (size_t k = j + 1; k < A->rowsCount; k++) {
			L->data[j][k] = 0.0;
		}
	}

	L->isSPD = isSPD;

	return L;
}

/**
 \fn		Mat Solve_Cholesky (Mat L, Mat B)

 \brief		Solves system of linear equations using Cholesky decomposition. The Cholesky
			decomposition is mainly used for the numerical solution of linear equations Ax = b. If A
			is symmetric and positive definite, then we can solve Ax = b by first computing the
			Cholesky decomposition A = LL*, then solving Ly = b for y by forward substitution, and
			finally solving L*x = y for x by back substitution.

 \param	L	The Matrix processed by Dcmp_Cholesky (Mat A).
 \param	B	The Matrix containing right hand side to solve system over.

 \return	Matrix contatining system solution as column-vector.
 */
Mat Solve_Cholesky (Mat L, Mat B) {
	Assert$(L->isSPD, "Matrix A must be symmetric positive definite.");
	Assert$(L->rowsCount == B->rowsCount, "Rows count in A and B must be equal.");

	Mat X = DeepCopy(B);

	//Solve L*Y = B;
	for (size_t k = 0; k < L->rowsCount; k++) {
		for (size_t j = 0; j < B->colsCount; j++) {
			for (size_t i = 0; i < k; i++) {
				X->data[k][j] -= X->data[i][j] * L->data[k][i];
			}

			X->data[k][j] /= L->data[k][k];
		}
	}

	//Solve L'*X = Y;
	for (ssize_t k = L->rowsCount - 1; k >= 0; k--) {
		for (size_t j = 0; j < B->colsCount; j++) {
			for (size_t i = (size_t) (k + 1); i < L->rowsCount; i++) {
				X->data[k][j] -= X->data[i][j] * L->data[i][k];
			}

			X->data[k][j] /= L->data[k][k];
		}
	}

	return X;
}

/**
 \fn		double Det_cholesky (Mat L)

 \brief		Computes matrix determinant using Cholesky dcmp.

 \param	L	The Lower-triangular Matrix that is a result of Dcmp_Cholesky (Mat A).

 \return	Determinant value.
 */
entry_type Det_cholesky (Mat L) {
	entry_type det = 1.0;

	for (size_t i = 0; i < L->rowsCount; i++) {
		det *= L->data[i][i];
	}

	return det;
}
