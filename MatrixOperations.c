﻿#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#include "MatrixOperations.h"
#include "Matrix.h"
#include "Gauss.h"
#include "Const.h"
#include "Extra.h"
#include "SpinningIndicator.h"


#pragma region "Entrywise operations"
TElementWise_MatrixMatrix1$(add, +); TElementWise_MatrixMatrix1$(sub, -); TElementWise_MatrixMatrix1$(mul, *); TElementWise_MatrixMatrix1$(div, /);
TElementWise_MatrixScalar$(add, +); TElementWise_MatrixScalar$(sub, -); TElementWise_MatrixScalar$(mul, *); TElementWise_MatrixScalar$(div, / );
TElementWise_MatrixMatrix2$(add, +); TElementWise_MatrixMatrix2$(sub, -); TElementWise_MatrixMatrix2$(mul, *); TElementWise_MatrixMatrix2$(div, /);

Mat MatEntryWiseLerp(Mat A, Mat B, double t) {
	Mat L = AllocMat(A->rowsCount, A->colsCount);

	for (size_t i = 0; i < L->rowsCount; i++) {
		for (size_t j = 0; j < L->colsCount; j++) {
			L->a[i][j] = lerp(A->a[i][j], B->a[i][j], t);
		}
	}

	return L;
}
#pragma endregion "Entrywise operations"


#pragma region "Is?.."
/**
 \fn	bool IsEqual (Mat A, Mat B)

 \brief	Checks for element-wise equality of matrices A and B.

 \date	04-Jun-14

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	true if equal, else false.
 */
bool IsEqual (Mat A, Mat B) {
	double **a = A->a;
	double **b = B->a;

	if ((A->rowsCount != B->rowsCount) || (A->colsCount != B->colsCount)) {
		return false;
	}

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			if (!equals_d(a[i][j], b[i][j])) {
				return false;
			}
		}
	}

	return true;
}

/**
 \fn	bool IsIdentity (Mat A)

 \brief	Check$ if matrix A is an identity matrix.

 \date	17-May-14

 \param	A	The Mat to process.

 \return	Checking result (true or false).
 */
bool IsIdentity (Mat A) {
	double **a = A->a;
	
	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			if ((i != j) && (fabs(a[i][j]) > EPS)) {
				return false;
			} else {
				if ((i == j) && (equals_d(a[i][j], 1.0))) {
					return false;
				}
			}
		}
	}

	return true;
}

bool IsSingular (Mat A) {
	if (A->isSingular) {
		return true;
	} else {
		Mat T = DeepCopy(A);
		toRowEchelonForm(T);
		bool r = T->isSingular;
		freeMat$(T);
		A->isSingular = r;
		return r;
	}
}

/**
 \fn	bool IsSymmetric (Mat A)

 \brief	Checks if Matrix A is symmetric.

 \param	A	The Mat to process.

 \return	true or false.
 */
bool IsSymmetric (Mat A) {
	if (!isSquare$(A)) {
		return false;
	}

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < i; j++) {
			if (!equals_d(A->a[i][j], A->a[j][i])) {
				return false;
			}
		}
	}

	return true;
}

bool IsSkewSymmetric (Mat A) {
	if (!isSquare$(A)) {
		return false;
	}

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < i; j++) {
			if (!equals_d(A->a[i][j], -A->a[j][i])) {
				return false;
			}
		}
	}

	return true;
}
#pragma endregion "Is?.."


#pragma region "Transpose & Inverse"

/**
 \fn	void toTransposed_square (Mat A)

 \brief	In-place square matrix transposition.

 \date	04-Jun-14

 \param	A	The Mat to process.
 */
void toTransposed_square (Mat A) {
	Assert$(isSquare$(A), "Cannot transpose non-square matrix with this func.");

	for (size_t i = 0; i < A->rowsCount - 1; i++) {
		for (size_t j = i + 1; j < A->rowsCount; j++) {
			swap_d(A->a[i][j], A->a[j][i]);
		}
	}

	return;
}

/**
 \fn	Mat Transposed (Mat A)

 \brief	Returns A transposed.

 \date	04-Jun-14

 \param	A	The Mat to process.

 \return	A transposed.
 */
Mat Transposed (Mat A) {
	double **a = A->a;
	Mat T = AllocMat(A->colsCount, A->rowsCount);
	double **t = T->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			t[j][i] = a[i][j];
		}
	}

	return T;
}

/**
 \fn	void toTransposed (Mat *A)

 \brief	(In-place) Transposes the given Matrix A.

 \param	A	The Mat to process.
 */
void toTransposed (Mat *A) {
	if ((*A)->rowsCount == (*A)->colsCount) {
		toTransposed_square(*A);
	} else {
		Mat T = Transposed(*A);
		freeMat$(*A);
		*A = DeepCopy(T);
		freeMat$(T);
	}

	return;
}

/**
 \fn	Mat Inverse (Mat A)

 \brief	Returns the Inverse matrix of the given Matrix A.

 \param	A	The Mat to process.

 \return	A^(-1).
 */
Mat Inverse (Mat A) {
//	if ((A->isSingular) || (fabs(Det_gauss(A)) <= EPS) || (!isSquare$(A))) {
//		return NULL;
//	}

    if (!Check$(!((A->isSingular) || (fabs(Det_gauss(A)) <= EPS) || (!isSquare$(A))),
               "Cannot invert singular matrix.")) {
        return NULL;
    }

	double **a = A->a;
	Mat R = NULL, I = NULL;

	switch (A->rowsCount) {
		case 1:
			R = AllocMat(A->rowsCount, A->colsCount);
			R->a[0][0] = 1.0 / a[0][0];
			return R;
			break;
		case 2:
			R = AllocMat(A->rowsCount, A->colsCount);
			R->a[1][0] = a[1][0] * -1.0;
			R->a[0][1] = a[0][1] * -1.0;
			R->a[0][0] = a[1][1];
			R->a[1][1] = a[0][0];
			_ms_mul(R, 1.0 / A->det); //TODO: A->det is already computed at the very beginning of func when
			                        //it checks for matrix singularity. So if checks are disabled, 
			                        //you need to compute Det(A) manually)  
            return R;
            break;
        case 3:
			R = AllocMat(A->rowsCount, A->colsCount);
			R->a[0][0] = a[1][1] * a[2][2] - a[2][1] * a[1][2];
			R->a[0][1] = a[0][2] * a[2][1] - a[0][1] * a[2][2];
			R->a[0][2] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
			R->a[1][0] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
			R->a[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0];
			R->a[1][2] = a[1][0] * a[0][2] - a[0][0] * a[1][2];
			R->a[2][0] = a[1][0] * a[2][1] - a[2][0] * a[1][1];
			R->a[2][1] = a[2][0] * a[0][1] - a[0][0] * a[2][1];
			R->a[2][2] = a[0][0] * a[1][1] - a[1][0] * a[0][1];
			_ms_mul(R, 1.0 / A->det);
			return R;
			break;
		case 4:
			R = AllocMat(A->rowsCount, A->colsCount);
			R->a[0][0] = a[1][2] * a[2][3] * a[3][1] - a[1][3] * a[2][2] * a[3][1] + a[1][3] * a[2][1] * a[3][2] - a[1][1] * a[2][3] * a[3][2] - a[1][2] * a[2][1] * a[3][3] + a[1][1] * a[2][2] * a[3][3];
			R->a[0][1] = a[0][3] * a[2][2] * a[3][1] - a[0][2] * a[2][3] * a[3][1] - a[0][3] * a[2][1] * a[3][2] + a[0][1] * a[2][3] * a[3][2] + a[0][2] * a[2][1] * a[3][3] - a[0][1] * a[2][2] * a[3][3];
			R->a[0][2] = a[0][2] * a[1][3] * a[3][1] - a[0][3] * a[1][2] * a[3][1] + a[0][3] * a[1][1] * a[3][2] - a[0][1] * a[1][3] * a[3][2] - a[0][2] * a[1][1] * a[3][3] + a[0][1] * a[1][2] * a[3][3];
			R->a[0][3] = a[0][3] * a[1][2] * a[2][1] - a[0][2] * a[1][3] * a[2][1] - a[0][3] * a[1][1] * a[2][2] + a[0][1] * a[1][3] * a[2][2] + a[0][2] * a[1][1] * a[2][3] - a[0][1] * a[1][2] * a[2][3];
			R->a[1][0] = a[1][3] * a[2][2] * a[3][0] - a[1][2] * a[2][3] * a[3][0] - a[1][3] * a[2][0] * a[3][2] + a[1][0] * a[2][3] * a[3][2] + a[1][2] * a[2][0] * a[3][3] - a[1][0] * a[2][2] * a[3][3];
			R->a[1][1] = a[0][2] * a[2][3] * a[3][0] - a[0][3] * a[2][2] * a[3][0] + a[0][3] * a[2][0] * a[3][2] - a[0][0] * a[2][3] * a[3][2] - a[0][2] * a[2][0] * a[3][3] + a[0][0] * a[2][2] * a[3][3];
			R->a[1][2] = a[0][3] * a[1][2] * a[3][0] - a[0][2] * a[1][3] * a[3][0] - a[0][3] * a[1][0] * a[3][2] + a[0][0] * a[1][3] * a[3][2] + a[0][2] * a[1][0] * a[3][3] - a[0][0] * a[1][2] * a[3][3];
			R->a[1][3] = a[0][2] * a[1][3] * a[2][0] - a[0][3] * a[1][2] * a[2][0] + a[0][3] * a[1][0] * a[2][2] - a[0][0] * a[1][3] * a[2][2] - a[0][2] * a[1][0] * a[2][3] + a[0][0] * a[1][2] * a[2][3];
			R->a[2][0] = a[1][1] * a[2][3] * a[3][0] - a[1][3] * a[2][1] * a[3][0] + a[1][3] * a[2][0] * a[3][1] - a[1][0] * a[2][3] * a[3][1] - a[1][1] * a[2][0] * a[3][3] + a[1][0] * a[2][1] * a[3][3];
			R->a[2][1] = a[0][3] * a[2][1] * a[3][0] - a[0][1] * a[2][3] * a[3][0] - a[0][3] * a[2][0] * a[3][1] + a[0][0] * a[2][3] * a[3][1] + a[0][1] * a[2][0] * a[3][3] - a[0][0] * a[2][1] * a[3][3];
			R->a[2][2] = a[0][1] * a[1][3] * a[3][0] - a[0][3] * a[1][1] * a[3][0] + a[0][3] * a[1][0] * a[3][1] - a[0][0] * a[1][3] * a[3][1] - a[0][1] * a[1][0] * a[3][3] + a[0][0] * a[1][1] * a[3][3];
			R->a[2][3] = a[0][3] * a[1][1] * a[2][0] - a[0][1] * a[1][3] * a[2][0] - a[0][3] * a[1][0] * a[2][1] + a[0][0] * a[1][3] * a[2][1] + a[0][1] * a[1][0] * a[2][3] - a[0][0] * a[1][1] * a[2][3];
			R->a[3][0] = a[1][2] * a[2][1] * a[3][0] - a[1][1] * a[2][2] * a[3][0] - a[1][2] * a[2][0] * a[3][1] + a[1][0] * a[2][2] * a[3][1] + a[1][1] * a[2][0] * a[3][2] - a[1][0] * a[2][1] * a[3][2];
			R->a[3][1] = a[0][1] * a[2][2] * a[3][0] - a[0][2] * a[2][1] * a[3][0] + a[0][2] * a[2][0] * a[3][1] - a[0][0] * a[2][2] * a[3][1] - a[0][1] * a[2][0] * a[3][2] + a[0][0] * a[2][1] * a[3][2];
			R->a[3][2] = a[0][2] * a[1][1] * a[3][0] - a[0][1] * a[1][2] * a[3][0] - a[0][2] * a[1][0] * a[3][1] + a[0][0] * a[1][2] * a[3][1] + a[0][1] * a[1][0] * a[3][2] - a[0][0] * a[1][1] * a[3][2];
			R->a[3][3] = a[0][1] * a[1][2] * a[2][0] - a[0][2] * a[1][1] * a[2][0] + a[0][2] * a[1][0] * a[2][1] - a[0][0] * a[1][2] * a[2][1] - a[0][1] * a[1][0] * a[2][2] + a[0][0] * a[1][1] * a[2][2];
			_ms_mul(R, 1.0 / A->det);
			return R;
			break;
		default:
			R = DeepCopy(A);
			I = Identity(A->rowsCount);

			concat(R, I);
			toReducedRowEchelonForm(R);

			for (size_t i = 0; i < A->rowsCount; i++) {
				for (size_t j = 0; j < A->colsCount; j++) {
					I->a[i][j] = R->a[i][j + A->colsCount];
				}
			}

			freeMat$(R);
			return I;
			break;
	}
}

/**
 \fn	void toInverse (Mat *A)

 \brief	(In-place) Transforms square matrix A with size n into an inverse matrix A^(-1)
		using row operations (Gauss-Jordan method) for n>4.

 \date	24-May-14

 \param	A	The Mat to process.
 */
void toInverse (Mat *A) {	
	Mat I = Inverse(*A);
	freeMat$(*A);
	*A = I;

	return;
}
#pragma endregion "Transpose & Inverse"


#pragma region "Matrix multiplication"
/**
 \fn	Mat MatMul_naive (Mat A, Mat B)

 \brief	Matrix multiplication using naive, but cache-friendly method))
		ijk - 1.25 cache misses per iteration ikj - 0.5 cache misses (row-wise)

 \date	24-May-14

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	Matrix product of A & B.
 */
Mat MatMul_naive (Mat A, Mat B) {
	Assert$(A->colsCount == B->rowsCount, "Cannot multiply.");

	Mat C = AllocMat(A->rowsCount, B->colsCount);

	if ((isSquare$(A)) && (A->rowsCount == 1) && (isSquare$(B))) {
		C->a[0][0] = A->a[0][0] * B->a[0][0];
		return C;
	} else {
		if ((isSquare$(A)) && (A->rowsCount == 2) && (isSquare$(B))) {
			C->a[0][0] = A->a[0][0] * B->a[0][0] + A->a[0][1] * B->a[1][0];
			C->a[0][1] = A->a[0][0] * B->a[0][1] + A->a[0][1] * B->a[1][1];
			C->a[1][0] = A->a[1][0] * B->a[0][0] + A->a[1][1] * B->a[1][0];
			C->a[1][1] = A->a[1][0] * B->a[0][1] + A->a[1][1] * B->a[1][1];
			return C;
		} else {
			for (size_t i = 0; i < A->rowsCount; ++i) {
				for (size_t k = 0; k < A->colsCount; ++k) {
					double s = A->a[i][k];
					for (size_t j = 0; j < B->colsCount; ++j) {
						C->a[i][j] += s * B->a[k][j];
					}
				}
			}

			return C;
		}
	}
}

/**
 \fn	Mat MatMul_naive_recursive (Mat A, Mat B)

 \brief	Recursive implementation of naive matrix multiplication algorithm.
 Only for square Matrices with size=2^n. Not memory	efficient!
 A, B & product C will be split into 4 submatrices, then the product will be:
		C11 = A11B11 + A12B21
		C12 = A11B12 + A12B22
		C21 = A21B11 + A22B21
		C22 = A21B12 + A22B22,
		or:
		C11 ← A11×B11;  C11 ← C11 + A12×B21
		C12 ← A11×B12;  C12 ← C12 + A12×B22
		C21 ← A21×B11;  C21 ← C21 + A22×B21
		C22 ← A21×B12;  C22 ← C22 + A22×B22

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	Product of A & B.
 */
Mat MatMul_naive_recursive (Mat A, Mat B) {
	Mat C = NULL;

	if (A->rowsCount == 1) {
		C = AllocMat(1, 1);
		C->a[0][0] = A->a[0][0] * B->a[0][0];
		return C;
	} else {
		if (A->rowsCount < 3) {
			C = AllocMat(2, 2);
			C->a[0][0] = A->a[0][0] * B->a[0][0] + A->a[0][1] * B->a[1][0];
			C->a[0][1] = A->a[0][0] * B->a[0][1] + A->a[0][1] * B->a[1][1];
			C->a[1][0] = A->a[1][0] * B->a[0][0] + A->a[1][1] * B->a[1][0];
			C->a[1][1] = A->a[1][0] * B->a[0][1] + A->a[1][1] * B->a[1][1];
			return C;
		} else {
			if (A->rowsCount < MM_SIZE_THRESHOLD) {
				C = MatMul_naive(A, B);
				return C;
			} else {
				size_t Size = A->rowsCount;
				Mat tmp_1 = NULL;

				C = AllocMat(Size, Size);

				Size /= 2;

				Mat A11 = AllocMat(Size, Size);
				Mat A12 = AllocMat(Size, Size);
				Mat A21 = AllocMat(Size, Size);
				Mat A22 = AllocMat(Size, Size);

				Mat B11 = AllocMat(Size, Size);
				Mat B12 = AllocMat(Size, Size);
				Mat B21 = AllocMat(Size, Size);
				Mat B22 = AllocMat(Size, Size);

				for (size_t i = 0; i < Size; i++) {
					for (size_t j = 0; j < Size; j++) {
						A11->a[i][j] = A->a[i][j];
						A12->a[i][j] = A->a[i][j + Size];
						A21->a[i][j] = A->a[i + Size][j];
						A22->a[i][j] = A->a[i + Size][j + Size];

						B11->a[i][j] = B->a[i][j];
						B12->a[i][j] = B->a[i][j + Size];
						B21->a[i][j] = B->a[i + Size][j];
						B22->a[i][j] = B->a[i + Size][j + Size];
					}
				}

				//TODO: C22 can be used as temporary, so no tmp_1 needed!)
				tmp_1 = MatMul_naive_recursive(A12, B21);
				Mat C11 = MatMul_naive_recursive(A11, B11);
				_mm1_add(C11, tmp_1);
				freeMat$(tmp_1);

				tmp_1 = MatMul_naive_recursive(A12, B22);
				Mat C12 = MatMul_naive_recursive(A11, B12);
				_mm1_add(C12, tmp_1);
				freeMat$(tmp_1);

				tmp_1 = MatMul_naive_recursive(A22, B21);
				Mat C21 = MatMul_naive_recursive(A21, B11);
				_mm1_add(C21, tmp_1);
				freeMat$(tmp_1);

				tmp_1 = MatMul_naive_recursive(A22, B22);
				Mat C22 = MatMul_naive_recursive(A21, B12);
				_mm1_add(C22, tmp_1);
				freeMat$(tmp_1);

				for (size_t i = 0; i < Size; i++) {
					for (size_t j = 0; j < Size; j++) {
						C->a[i][j] = C11->a[i][j];
						C->a[i][j + Size] = C12->a[i][j];
						C->a[i + Size][j] = C21->a[i][j];
						C->a[i + Size][j + Size] = C22->a[i][j];
					}
				}

				freeMats(A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22, NULL);

				return C;
			}
		}
	}
}

/**
 \fn	void matMul(Mat *A, Mat B)

 \brief	In-place matrix multiplication (multiplicand is replaced by result).

 \date	24-May-14

 \param	A	The Mat A to process (will be replaced by product of A &amp; B).
 \param	B	The Mat B to process.
 */
void matMul (Mat *A, Mat B) {
	Assert$((*A)->colsCount == B->rowsCount,
		"Cannot multiply matrices\n(Number of columns of A must be equal to number of rows of B).");

	Mat P = MatMul(*A, B);
	freeMat$(*A);
	*A = DeepCopy(P);
	freeMat$(P);

	return;
}

/**
 \fn	Mat MatPow (Mat A, size_t pow)

 \brief	Raise matrix to power N.		
		TODO: use addition - chain exponentiation.

 \date	24-May-14
 
 \param	A  	The Mat to process.
 \param	pow	The power value to raise matrix to.

 \return	A^n.
 */
Mat MatPow (Mat A, size_t pow) {
	Assert$(pow >= 0, "Cannot raise power.");
	Assert$(isSquare$(A), "Cannot raise power of non-square matrix.");

	size_t c = 0;
	Mat R;

	switch (pow) {
		case 0:
			return Identity(A->rowsCount);
			break;
		case 1:
			R = DeepCopy(A);
			return R;
			break;
		case 2:
			R = DeepCopy(A);
			matMul(&R, A);

			return R;
			break;
		default:
			R = DeepCopy(A);
			do {
				spinActivityIndicator();
				matMul(&R, A);
				c++;
			} while (c < pow-1);
			clearActivityIndicator();

			return R;
			break;
	}
}
#pragma endregion "Matrix multiplication"


#pragma region "Norms, etc."

/**
 \fn	size_t Rank (Mat A)

 \brief	Calculates Rank of matrix A in RREF form.

 \date	17-May-14

 \param	A	The dMat to process.

 \return	Rank value.
 */
size_t Rank (Mat A) {
	size_t rank = 0;
	size_t i, j;
	double **r = A->a;

	for (i = 0; i < A->rowsCount; i++) {
		for (j = 0; j < A->colsCount; j++) {
			if (fabs(r[i][j]) > EPS) {
				rank++;
				break;
			}
		}
	}

	return rank;
}

/**
 \fn	double Trace (Mat A)

 \brief	Computes trace of A.

 \param	A	The Mat to process.

 \return	Trace value.
 */
double Trace (Mat A) {
	double **a = A->a;	
	double tr = 0.0;

	for (size_t i = 0; i < min(A->rowsCount, A->colsCount); i++) {
		tr += a[i][i];
	}
	A->trace = tr;

	return A->trace;
}

/**
 \fn	double OneNorm (Mat A)

 \brief	One-norm of A (maximum column sum).

 \param	A	The Mat to process.

 \return	1-norm.
 */
double OneNorm (Mat A) {
	double norm = 0.0;

	for (size_t i = 0; i < A->colsCount; i++) {
		double sum = 0.0;
		for (size_t j = 0; j < A->rowsCount; j++) {
			sum += fabs(A->a[j][i]);
		}
		norm = max(norm, sum);
	}

	return norm;
}

double TwoNorm (Mat A) {
	return EuclideanNorm(A); //TODO:
}

/**
 \fn	double InfinityNorm (Mat A)

 \brief	Infinity-norm of A (max row sum).

 \param	A	The Mat to process.

 \return	Inf-norm.
 */
double InfinityNorm (Mat A) {
	double norm = 0.0;

	for (size_t i = 0; i < A->rowsCount; i++) {
		double sum = 0.0;
		for (size_t j = 0; j < A->colsCount; j++) {
			sum += fabs(A->a[i][j]);
		}
		norm = max(norm, sum);
	}

	return norm;
}

/**
 \fn	double EuclideanNorm (Mat A)

 \brief	Euclidean norm of A.

 \param	A	The Mat to process.

 \return	Eucl. norm.
 */
double EuclideanNorm (Mat A) {
	double sum = 0.0;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			sum += square_d(A->a[i][j]);
		}
	}

	return sqrt(sum);
}

/**
 \fn	double ConditionNumber (Mat A)

 \brief	Condition number of operator A.

 \param	A	The Mat to process.

 \return	Condition number.
 */
double ConditionNumber (Mat A) {
	Mat Ai = DeepCopy(A);
	toInverse(&Ai);
	double r = InfinityNorm(A) * InfinityNorm(Ai);
	freeMat$(Ai);

	return r;
}
#pragma endregion "Norms, etc."


#pragma region "Kronecker"

/**
 \fn	Mat KroneckerProd (Mat A, Mat B)

 \brief	Kronecker product of A & B.

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	A(x)B.
 */
Mat KroneckerProd (Mat A, Mat B) {
	Mat K = AllocMat(A->rowsCount*B->rowsCount, A->colsCount*B->colsCount);
	double **a = A->a;
	double **b = B->a;
	double **k = K->a;

	for (size_t i = 0; i < K->rowsCount; i++) {
		for (size_t j = 0; j < K->colsCount; j++) {
			k[i][j] = (a[i/A->rowsCount][j/A->colsCount] *
				b[i%B->rowsCount][j%B->colsCount]);
		}
	}

	return K;
}

/**
 \fn	Mat KroneckerSum (Mat A, Mat B)

 \brief	Kronecker sum of A & B.

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	A(+)B.
 */
Mat KroneckerSum (Mat A, Mat B) {
	Assert$(isSquare$(A) && isSquare$(B), "A and B isn't square.");
	Mat Ib = Identity(B->rowsCount);
	Mat Ia = Identity(A->rowsCount);
	Mat AI = KroneckerProd(A, Ib);
	Mat IB = KroneckerProd(Ia, B);
	Mat KS = AllocMat(A->rowsCount*B->rowsCount, A->colsCount*B->colsCount);
	_mm2_add(AI, IB, KS);
	freeMat$(IB);
	freeMat$(Ia);
	freeMat$(AI);
	freeMat$(Ib);

	return KS;
}
#pragma endregion "Kronecker"


#pragma region "Strassen"
/**
 \fn	size_t _fixSize (size_t Size)

 \brief	Fix size.

 \param	Size	The size.

 \return	A size_t.
 */
size_t _fixSize(size_t Size) {
	if (!(Check$(ispowerof2_i(Size), "Matrix size is not a power of 2."))) {
		return (size_t) ((int64_t) (1) << (int64_t) (ceil(log2(Size)))); //-V113
	}
	return Size;
}

/**
 \fn	Mat MatMul_strassen (Mat A, Mat B)

 \brief	Finds product of A & B using Strassen algorithm.
 Source matrices A, B & its product C will be divided into 4 square blocks (submatrices),
 and this algorithm repeated recursively until blocks become numbers. 
 (or when size of blocks reaches some threshold value when naive algorithm will be used).
 Complexity: O(n^2.8)
 TODO: convert input matrices into 'recursion-friendly' form (e.g 'matrix if matrices')
		M1 := (A1,1 + A2,2)(B1,1 + B2,2)
		M2 := (A2,1 + A2,2)B1,1 M3 := A1,1(B1,2 − B2,2)
		M4 := A2,2(B2,1 − B1,1)
		M5 := (A1,1 + A1,2)B2, 2 M6 := (A2,1 − A1,1)(B1,1 + B1,2)
		M7 := (A1,2 − A2,2)(B2,1 + B2,2)
		C1,1 = M1 + M4 − M5 + M7
		C1,2 = M3 + M5
		C2,1 = M2 + M4
		C2,2 = M1 − M2 + M3 + M6.

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	Product of A & B.
 */
Mat MatMul_strassen (Mat A, Mat B) {
	Mat C = NULL;
	size_t Size = A->rowsCount;
	
	switch (Size) {
		case 1:
			C = AllocMat(Size, Size);
			C->a[0][0] = A->a[0][0] * B->a[0][0];
			return C;
			break;
		case 2:
			C = AllocMat(Size, Size);
			C->a[0][0] = A->a[0][0] * B->a[0][0] + A->a[0][1] * B->a[1][0];
			C->a[0][1] = A->a[0][0] * B->a[0][1] + A->a[0][1] * B->a[1][1];
			C->a[1][0] = A->a[1][0] * B->a[0][0] + A->a[1][1] * B->a[1][0];
			C->a[1][1] = A->a[1][0] * B->a[0][1] + A->a[1][1] * B->a[1][1];
            return C;
            break;
        case 3: case 4: case 5: case 6: case 7: case 8: //TODO: fuckin' MSVC does not support C99 case ranges upd(seems that it is GCC extension or smth)
			C = MatMul(A, B);							//TODO: use MM_SIZE_THRESHOLD
			return C;
			break;
		default:
			C = AllocMat(Size, Size);

			Size /= 2;

			Mat A11 = AllocMat(Size, Size);
			Mat A12 = AllocMat(Size, Size);
			Mat A21 = AllocMat(Size, Size);
			Mat A22 = AllocMat(Size, Size);

			Mat B11 = AllocMat(Size, Size);
			Mat B12 = AllocMat(Size, Size);
			Mat B21 = AllocMat(Size, Size);
			Mat B22 = AllocMat(Size, Size);

			Mat C11 = AllocMat(Size, Size);
			Mat C12 = AllocMat(Size, Size);
			Mat C21 = AllocMat(Size, Size);
			Mat C22 = AllocMat(Size, Size);

			Mat A_tmp = AllocMat(Size, Size);
			Mat B_tmp = AllocMat(Size, Size);

			// Split
			for (size_t i = 0; i < Size; i++) {
				for (size_t j = 0; j < Size; j++) {
					A11->a[i][j] = A->a[i][j];
					A12->a[i][j] = A->a[i][j + Size];
					A21->a[i][j] = A->a[i + Size][j];
					A22->a[i][j] = A->a[i + Size][j + Size];

					B11->a[i][j] = B->a[i][j];
					B12->a[i][j] = B->a[i][j + Size];
					B21->a[i][j] = B->a[i + Size][j];
					B22->a[i][j] = B->a[i + Size][j + Size];
				}
			}

            _mm2_add(A11, A22, A_tmp); // A11 + A22
            _mm2_add(B11, B22, B_tmp); // B11 + B22
			Mat p1 = MatMul_strassen(A_tmp, B_tmp); // p1 = (A11+A22) * (B11+B22)

			_mm2_add(A21, A22, A_tmp); // A21 + A22
			Mat p2 = MatMul_strassen(A_tmp, B11); // p2 = (A21+A22) * (B11)

            _mm2_sub(B12, B22, B_tmp); // B12 - B22
			Mat p3 = MatMul_strassen(A11, B_tmp); // p3 = (A11) * (B12 - B22)

            _mm2_sub(B21, B11, B_tmp); // B21 - B11
			Mat p4 = MatMul_strassen(A22, B_tmp); // p4 = (A22) * (B21 - B11)

            _mm2_add(A11, A12, A_tmp); // A11 + A12
			Mat p5 = MatMul_strassen(A_tmp, B22); // p5 = (A11+A12) * (B22)   

            _mm2_sub(A21, A11, A_tmp); // A21 - A11
            _mm2_add(B11, B12, B_tmp); // B11 + B12
			Mat p6 = MatMul_strassen(A_tmp, B_tmp); // p6 = (A21-A11) * (B11+B12)

            _mm2_sub(A12, A22, A_tmp); // A12 - A22
            _mm2_add(B21, B22, B_tmp); // B21 + B22
			Mat p7 = MatMul_strassen(A_tmp, B_tmp); // p7 = (A12-A22) * (B21+B22)

            _mm2_add(p3, p5, C12); // c12 = p3 + p5
            _mm2_add(p2, p4, C21); // c21 = p2 + p4

            _mm2_add(p1, p4, A_tmp); // p1 + p4
            _mm2_add(A_tmp, p7, B_tmp); // p1 + p4 + p7
            _mm2_sub(B_tmp, p5, C11); // c11 = p1 + p4 - p5 + p7

            _mm2_add(p1, p3, A_tmp); // p1 + p3
            _mm2_add(A_tmp, p6, B_tmp); // p1 + p3 + p6
			_mm2_sub(B_tmp, p2, C22); // c22 = p1 + p3 - p2 + p6

			// Join
			for (size_t i = 0; i < Size; i++) {
				for (size_t j = 0; j < Size; j++) {
					C->a[i][j] = C11->a[i][j];
					C->a[i][j + Size] = C12->a[i][j];
					C->a[i + Size][j] = C21->a[i][j];
					C->a[i + Size][j + Size] = C22->a[i][j];
				}
			}

			freeMats(A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22,
				A_tmp, B_tmp, p1, p2, p3, p4, p5, p6, p7, NULL);

			return C;
			break;
	}
}

/**
 \fn	Mat MatMul_strassen_optimized (Mat A, Mat B)

 \brief	P1 = (A01+A10)(B10+B01)
		P2 = (A10+A11)B10 P3 = A01(B11-B01)
		P4 = A10(B00-B10)
		P5 = (A01+A00)B01 P6 = (A11-A01)(B10+B11)
		P7 = (A00-A10)(B00+B01)
		C00 = P1+P4-P5+P7 C01 = P3+P5 C10 = P2+P4 C11 = P1-P2+P3+P6.

 \param	A	The Mat to process.
 \param	B	The Mat to process.

 \return	A Mat.
 */
Mat MatMul_strassen_optimized (Mat A, Mat B) {
	return MatMul_strassen(A, B); //TODO:
}
#pragma endregion "Strassen"
