#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "Dcmp_QR.h"
#include "Matrix.h"
#include "MatrixOps.h"
#include "Extras.h"
#include "Maths.h"


/*In linear algebra, a QR decomposition (also called a QR factorization) of a matrix
is a decomposition of a matrix A into a product A = QR of an orthogonal matrix Q
and an upper triangular matrix R.
QR decomposition is often used to solve the linear least squares problem,
and is the basis for a particular eigenvalue algorithm, the QR algorithm.

If A has n linearly independent columns, then the first n columns of Q form an orthonormal basis for the column space of A.
More specifically, the first k columns of Q form an orthonormal basis for the span of the first k columns of A for any 1 ≤ k ≤ n.
The fact that any column k of A only depends on the first k columns of Q is responsible for the triangular form of R.
*/
Mat *Dcmp_QR_Householder (Mat A) {
	size_t rows = A->rowsCount;
	size_t columns = A->colsCount;

	Mat QR = DeepCopy(A);
	Mat Q = AllocMat(rows, columns);
	Mat R = AllocMat(rows, columns);
	fill_zeroes(R);

	entry_t **qr = QR->data;

	for (size_t k = 0; k < columns; k++) {
		// Compute 2-norm of k-th column without under/overflow.
		entry_t norm = 0.0;

		for (size_t i = k; i < rows; i++) {
			norm = hypot(norm, qr[i][k]);
		}

		if (isnotzero(norm)) {
			// Form k-th Householder vector.
			if (qr[k][k] < (entry_t) 0.0) {
				norm = -norm;
			}

			for (size_t i = k; i < rows; i++) {
				qr[i][k] /= norm;
			}

			qr[k][k] += 1.0;

			// Apply transformation to remaining columns.
			for (size_t j = k + 1; j < columns; j++) {
				entry_t s = 0.0;

				for (size_t i = k; i < rows; i++) {
					s += qr[i][k] * qr[i][j];
				}

				s = -s / qr[k][k];

				for (size_t i = k; i < rows; i++) {
					qr[i][j] += s*(qr[i][k]);

					if (i < j) {
						R->data[i][j] = qr[i][j];
					}
				}
			}
		}

		R->data[k][k] = -norm;

		R->permutationSign *= -1;

		if (iszero(R->data[k][k])) {
			R->isRankDeficient = true;
		}
	}

	for (ssize_t k = columns - 1; k >= 0; k--) {
		for (size_t i = 0; i < rows; i++) {
			Q->data[i][k] = 0.0;
		}

		Q->data[k][k] = 1.0;

		for (size_t j = (size_t) k; j < columns; j++) {
			if (isnotzero(qr[k][k])) {
				entry_t s = 0.0;

				for (size_t i = (size_t) k; i < rows; i++) {
					s += qr[i][k] * Q->data[i][j];
				}

				s = -s / qr[k][k];

				for (size_t i = (size_t) k; i < rows; i++) {
					Q->data[i][j] += s * qr[i][k];
				}
			}
		}
	}

	freeMat$(QR);

	Mat *result = (Mat*) malloc(2 * sizeof(*result));
	Assert$(result != NULL, "Cannot allocate.");
	result[0] = Q;
	result[1] = R;

	return result;
}

entry_t Det_QR (Mat *QR) {
	entry_t det = (entry_t) QR[1]->permutationSign;
	entry_t **r = QR[1]->data;

	for (size_t i = 0; i < QR[1]->rowsCount; i++) {
		det *= r[i][i];
	}

	return det;
}

//Solve R*X=Q'*B
//Solve A*X = B, when A is represented by Q*R
Mat Solve_QR (Mat *QR, Mat B) {
	Assert$(QR[0]->rowsCount == B->rowsCount, "Rows count doesn't match.");
	Check$(QR[1]->isRankDeficient == false, "Rank deficient system.");

	Mat Qt = Transposed(QR[0]);
	matMul_inplace(&Qt, B);
	entry_t **qt = Qt->data;

	for (ssize_t k = QR[0]->colsCount - 1; k >= 0; k--) {
		for (size_t j = 0; j < B->colsCount; j++) {
			qt[k][j] /= QR[1]->data[k][k];
		}

		for (ssize_t i = 0; i < k; i++) {
			for (size_t j = 0; j < B->colsCount; j++) {
				qt[i][j] -= qt[k][j] * QR[1]->data[i][k];
			}
		}
	}

	return Qt;
}

//TODO: Solve X*A = B === A'*X' = B'
