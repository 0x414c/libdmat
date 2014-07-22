#include <stdbool.h>
#include <math.h>

#include "CholeskyDcmp.h"
#include "Matrix.h"
#include "Const.h"
#include "Extra.h"


Mat CholeskyDcmp (Mat A) {
	Assert(A->rowsCount == A->colsCount, "");

	Mat L = AllocMat(A->rowsCount, A->colsCount);
	bool isSPD = (A->colsCount == A->rowsCount);

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
		isSPD = isSPD && (d > 0.0);
		L->a[j][j] = sqrt(max(d, 0.0));
		for (size_t k = j + 1; k < A->rowsCount; k++) {
			L->a[j][k] = 0.0;
		}
	}
	L->isSPD = isSPD;

	return L;
}

Mat CholeskySolve (Mat L, Mat B) {
	Assert(L->rowsCount == B->rowsCount, "Rows count mismatch.");
	Assert(L->isSPD, "Matrix is not symmetric positive definite.");

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

double CholeskyDet (Mat L) {
	double det = 1.0;
	for (size_t i = 0; i < L->rowsCount; i++) {
		det *= L->a[i][i];
	}

	return det;
}