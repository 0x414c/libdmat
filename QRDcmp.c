#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "QRDcmp.h"
#include "Matrix.h"
#include "MatrixOperations.h"
#include "Extra.h"


/*In linear algebra, a QR decomposition (also called a QR factorization) of a matrix
is a decomposition of a matrix A into a product A = QR of an orthogonal matrix Q 
and an upper triangular matrix R. 
QR decomposition is often used to solve the linear least squares problem,
and is the basis for a particular eigenvalue algorithm, the QR algorithm.

If A has n linearly independent columns, then the first n columns of Q form an orthonormal basis for the column space of A.
More specifically, the first k columns of Q form an orthonormal basis for the span of the first k columns of A for any 1 ≤ k ≤ n.
The fact that any column k of A only depends on the first k columns of Q is responsible for the triangular form of R.
*/
Mat *QRDcmp_householder (Mat A) {
	size_t rows = A->rowsCount;
	size_t columns = A->colsCount;

	Mat QRMat = DeepCopy(A);
	Mat Q = AllocMat(rows, columns);
	Mat R = AllocMat(rows, columns);

	double **qrPtr = QRMat->a;

	for (size_t k = 0; k < columns; k++) {
		// Compute 2-norm of k-th column without under/overflow.
		double norm = 0;
		for (size_t i = k; i < rows; i++) {
			norm = hypot(norm, qrPtr[i][k]);
		}

		if (fabs(norm) > EPS) {
			// Form k-th Householder vector.
			if (qrPtr[k][k] < 0) {
				norm = -norm;
			}
			for (size_t i = k; i < rows; i++) {
				qrPtr[i][k] /= norm;
			}
			qrPtr[k][k] += 1.0;

			// Apply transformation to remaining columns.
			for (size_t j = k + 1; j < columns; j++) {
				double s = 0.0;
				for (size_t i = k; i < rows; i++) {
					s += qrPtr[i][k] * qrPtr[i][j];
				}
				s = -s / qrPtr[k][k];
				for (size_t i = k; i < rows; i++) {
					qrPtr[i][j] += s*(qrPtr[i][k]);
					
					if (i < j) {
						R->a[i][j] = qrPtr[i][j];
					}
				}
			}
		}
		R->a[k][k] = -norm;
		if (fabs(R->a[k][k]) <= EPS) {
			R->isRankDeficient = true;
		}
	}

	for (ptrdiff_t k = columns - 1; k >= 0; k--) {
		for (size_t i = 0; i < rows; i++) {
			Q->a[i][k] = 0.0;
		}
		Q->a[k][k] = 1.0;
		for (size_t j = k; j < columns; j++) {
			if (fabs(qrPtr[k][k]) > EPS) {
				double s = 0.0;
				for (size_t i = k; i < rows; i++) {
					s += qrPtr[i][k] * Q->a[i][j];
				}
				s = -s / qrPtr[k][k];
				for (size_t i = k; i < rows; i++) {
					Q->a[i][j] += s*(qrPtr[i][k]);
				}
			}
		}
	}

	FreeMat(QRMat);

	Mat *result = (Mat*) malloc(2 * sizeof(*result));
	Assert(result != NULL, "Cannot allocate.");
	result[0] = Q;
	result[1] = R;

	return result;
}

double Det_qr (Mat *qr) {
	double det = 1.0;
	
	for (size_t i = 0; i < qr[1]->rowsCount; i++) {
		det *= qr[1]->a[i][i];
	}

	return det;
}

//Solve R*X=Q'*B
//Solve A*X = B, when A represented by Q*R
Mat Solve_qr (Mat *qr, Mat B) {
	Assert(qr[0]->rowsCount == B->rowsCount, "Rows count doesn't match.");
	Check(qr[1]->isRankDeficient == false, "Rank deficient system.");

	Mat Qt = Transposed(qr[0]);
	matMul(&Qt, B);

	for (ptrdiff_t k = qr[0]->colsCount - 1; k >= 0; k--) {
		for (size_t j = 0; j < B->colsCount; j++) {
			Qt->a[k][j] /= qr[1]->a[k][k];
		}
		for (ptrdiff_t i = 0; i < k; i++) { //warning C4018: '<' : signed/unsigned mismatch
			for (size_t j = 0; j < B->colsCount; j++) {
				Qt->a[i][j] -= Qt->a[k][j] * qr[1]->a[i][k];
			}
		}
	}

	return Qt;
}

//Solve X*A = B === A'*X' = B'
//Mat QRSolve_t (Mat A, Mat B) {
//}
