#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#include "MatrixOps.h"
#include "Matrix.h"
#include "Const.h"
#include "Gauss.h"
#include "Extras.h"
#include "Maths.h"
#include "SpinningIndicator.h"


#pragma region "Entrywise operations"
EntryWiseFunc_MatrixMatrix_Inplace$(add, +)
EntryWiseFunc_MatrixMatrix_Inplace$(sub, -)
EntryWiseFunc_MatrixMatrix_Inplace$(mul, *)
EntryWiseFunc_MatrixMatrix_Inplace$(div, /)

EntryWiseFunc_MatrixMatrix$(add, +)
EntryWiseFunc_MatrixMatrix$(sub, -)
EntryWiseFunc_MatrixMatrix$(mul, *)
EntryWiseFunc_MatrixMatrix$(div, /)

EntryWiseFunc_MatrixScalar_Inplace$(add, +)
EntryWiseFunc_MatrixScalar_Inplace$(sub, -)
EntryWiseFunc_MatrixScalar_Inplace$(mul, *)
EntryWiseFunc_MatrixScalar_Inplace$(div, /)


Mat MatLerp_entrywise (Mat A, Mat B, entry_t t) {
	Mat L = AllocMat(A->rowsCount, A->colsCount);
	entry_t **a = A->data;
	entry_t **b = B->data;

	for (size_t i = 0; i < L->rowsCount; i++) {
		for (size_t j = 0; j < L->colsCount; j++) {
			L->data[i][j] = lerp(a[i][j], b[i][j], t);
		}
	}

	return L;
}
#pragma endregion "Entrywise operations"


#pragma region "Is?.."
/**
 \fn	bool IsEntriesEqual (Mat A, Mat B)

 \brief	Checks for element-wise equality of matrices A and B.

 \date	04-Jun-14

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	true if equal, else false.
 */
bool IsEntriesEqual (Mat A, Mat B) {
	entry_t **a = A->data;
	entry_t **b = B->data;

	if (!IsDimsEqual(A, B)) {
		return false;
	} else {
		for (size_t i = 0; i < A->rowsCount; i++) {
			for (size_t j = 0; j < A->colsCount; j++) {
				if (!(equals(a[i][j], b[i][j]))) {
					return false;
				}
			}
		}
	}

	return true;
}

inline bool IsDimsEqual (Mat A, Mat B) {
	return (A->rowsCount == B->rowsCount) && (A->colsCount == B->colsCount);
}

/**
 \fn	bool IsIdentity (Mat A)

 \brief	Checks if Matrix A is an identity matrix.

 \date	17-May-14

 \param	A	The Matrix to process.

 \return	Checking result (true or false).
 */
bool IsIdentity (Mat A) {
    Assert$(IsSquare$(A), "Matrix A should be square.");

	if (A->isIdentity) {
		return true;
	} else {
		entry_t **a = A->data;

		switch (A->rowsCount) {
			case 0:
				return false;
			case 1:
				return equals(a[0][0], 1.0);
			case 2:
				return equals(a[0][0], 1.0) && iszero(a[0][1]) && iszero(a[1][0]) && equals(a[1][1], 1.0);
			default:
				for (size_t i = 0; i < A->rowsCount; i++) {
					for (size_t j = 0; j < A->colsCount; j++) {
						if ((i != j) && (isnotzero(a[i][j]))) {
							return false;
						} else {
							if ((i == j) && (!(equals(a[i][j], 1.0)))) {
								return false;
							}
						}
					}
				}

				return true;
		}
	}
}

bool IsSingular (Mat A) {
	if (A->isSingular) {
		return true;
	} else {
		Mat T = DeepCopy(A);
		entry_t det = Det_Gauss(T);
		bool isSingular = T->isSingular;
		freeMat$(T);
		A->isSingular = isSingular;

		return isSingular;
	}
}

/**
 \fn	bool IsSymmetric (Mat A)

 \brief	Checks if Matrix A is symmetric.

 \param	A	The Mat to process.

 \return	true or false.
 */
bool IsSymmetric (Mat A) {
	if (!IsSquare$(A)) {
		return false;
	} else {
		entry_t **a = A->data;

		for (size_t i = 0; i < A->rowsCount; i++) {
			for (size_t j = 0; j < i; j++) {
				if (!equals(a[i][j], a[j][i])) {
					return false;
				}
			}
		}

		return true;
	}
}

bool IsSkewSymmetric (Mat A) {
	if (!IsSquare$(A)) {
		return false;
	} else {
		entry_t **a = A->data;

		for (size_t i = 0; i < A->rowsCount; i++) {
			for (size_t j = 0; j < i; j++) {
				if (!equals(a[i][j], -a[j][i])) {
					return false;
				}
			}
		}

		return true;
	}
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
	Assert$(IsSquare$(A), "Cannot transpose non-square matrix with this function. Use `toTransposed()' instead.");

	entry_t **a = A->data;

	for (size_t i = 0; i < A->rowsCount - 1; i++) {
		for (size_t j = i + 1; j < A->rowsCount; j++) {
			swap(a[i][j], a[j][i]);
		}
	}

	return;
}

/**
 \fn	Mat Transposed (Mat A)

 \brief	Returns A transposed.

 \date	04-Jun-14

 \param	A	The Mat to process.

 \return	Matrix A transposed.
 */
Mat Transposed (Mat A) {
	Mat T = AllocMat(A->colsCount, A->rowsCount);
	entry_t **t = T->data;
	entry_t **a = A->data;

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
	Assert$(IsSquare$(A), "Matrix A should be square.");

	entry_t det = Det_Gauss(A);
	if (!Check$(!(A->isSingular), "Cannot invert singular matrix.")) {
		return NULL;
	}

	Mat RI = NULL, I = NULL;
	entry_t **a = A->data;

	switch (A->rowsCount) {
		case 1:
			RI = AllocMat(A->rowsCount, A->colsCount);

			RI->data[0][0] = ((entry_t) 1.0) / A->det;

			return RI;
		case 2:
			RI = AllocMat(A->rowsCount, A->colsCount);

			RI->data[1][0] = -a[1][0];
			RI->data[0][1] = -a[0][1];
			RI->data[0][0] =  a[1][1];
			RI->data[1][1] =  a[0][0];

			_mul_m_s_inplace(RI, ((entry_t) 1.0) / A->det);

			return RI;
    	case 3:
			RI = AllocMat(A->rowsCount, A->colsCount);

			RI->data[0][0] = a[1][1] * a[2][2] - a[2][1] * a[1][2];
			RI->data[0][1] = a[0][2] * a[2][1] - a[0][1] * a[2][2];
			RI->data[0][2] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
			RI->data[1][0] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
			RI->data[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0];
			RI->data[1][2] = a[1][0] * a[0][2] - a[0][0] * a[1][2];
			RI->data[2][0] = a[1][0] * a[2][1] - a[2][0] * a[1][1];
			RI->data[2][1] = a[2][0] * a[0][1] - a[0][0] * a[2][1];
			RI->data[2][2] = a[0][0] * a[1][1] - a[1][0] * a[0][1];

			_mul_m_s_inplace(RI, ((entry_t) 1.0) / A->det);

			return RI;
		case 4:
			RI = AllocMat(A->rowsCount, A->colsCount);

			RI->data[0][0] = a[1][2] * a[2][3] * a[3][1] - a[1][3] * a[2][2] * a[3][1] + a[1][3] * a[2][1] * a[3][2] - a[1][1] * a[2][3] * a[3][2] - a[1][2] * a[2][1] * a[3][3] + a[1][1] * a[2][2] * a[3][3];
			RI->data[0][1] = a[0][3] * a[2][2] * a[3][1] - a[0][2] * a[2][3] * a[3][1] - a[0][3] * a[2][1] * a[3][2] + a[0][1] * a[2][3] * a[3][2] + a[0][2] * a[2][1] * a[3][3] - a[0][1] * a[2][2] * a[3][3];
			RI->data[0][2] = a[0][2] * a[1][3] * a[3][1] - a[0][3] * a[1][2] * a[3][1] + a[0][3] * a[1][1] * a[3][2] - a[0][1] * a[1][3] * a[3][2] - a[0][2] * a[1][1] * a[3][3] + a[0][1] * a[1][2] * a[3][3];
			RI->data[0][3] = a[0][3] * a[1][2] * a[2][1] - a[0][2] * a[1][3] * a[2][1] - a[0][3] * a[1][1] * a[2][2] + a[0][1] * a[1][3] * a[2][2] + a[0][2] * a[1][1] * a[2][3] - a[0][1] * a[1][2] * a[2][3];

			RI->data[1][0] = a[1][3] * a[2][2] * a[3][0] - a[1][2] * a[2][3] * a[3][0] - a[1][3] * a[2][0] * a[3][2] + a[1][0] * a[2][3] * a[3][2] + a[1][2] * a[2][0] * a[3][3] - a[1][0] * a[2][2] * a[3][3];
			RI->data[1][1] = a[0][2] * a[2][3] * a[3][0] - a[0][3] * a[2][2] * a[3][0] + a[0][3] * a[2][0] * a[3][2] - a[0][0] * a[2][3] * a[3][2] - a[0][2] * a[2][0] * a[3][3] + a[0][0] * a[2][2] * a[3][3];
			RI->data[1][2] = a[0][3] * a[1][2] * a[3][0] - a[0][2] * a[1][3] * a[3][0] - a[0][3] * a[1][0] * a[3][2] + a[0][0] * a[1][3] * a[3][2] + a[0][2] * a[1][0] * a[3][3] - a[0][0] * a[1][2] * a[3][3];
			RI->data[1][3] = a[0][2] * a[1][3] * a[2][0] - a[0][3] * a[1][2] * a[2][0] + a[0][3] * a[1][0] * a[2][2] - a[0][0] * a[1][3] * a[2][2] - a[0][2] * a[1][0] * a[2][3] + a[0][0] * a[1][2] * a[2][3];

			RI->data[2][0] = a[1][1] * a[2][3] * a[3][0] - a[1][3] * a[2][1] * a[3][0] + a[1][3] * a[2][0] * a[3][1] - a[1][0] * a[2][3] * a[3][1] - a[1][1] * a[2][0] * a[3][3] + a[1][0] * a[2][1] * a[3][3];
			RI->data[2][1] = a[0][3] * a[2][1] * a[3][0] - a[0][1] * a[2][3] * a[3][0] - a[0][3] * a[2][0] * a[3][1] + a[0][0] * a[2][3] * a[3][1] + a[0][1] * a[2][0] * a[3][3] - a[0][0] * a[2][1] * a[3][3];
			RI->data[2][2] = a[0][1] * a[1][3] * a[3][0] - a[0][3] * a[1][1] * a[3][0] + a[0][3] * a[1][0] * a[3][1] - a[0][0] * a[1][3] * a[3][1] - a[0][1] * a[1][0] * a[3][3] + a[0][0] * a[1][1] * a[3][3];
			RI->data[2][3] = a[0][3] * a[1][1] * a[2][0] - a[0][1] * a[1][3] * a[2][0] - a[0][3] * a[1][0] * a[2][1] + a[0][0] * a[1][3] * a[2][1] + a[0][1] * a[1][0] * a[2][3] - a[0][0] * a[1][1] * a[2][3];

			RI->data[3][0] = a[1][2] * a[2][1] * a[3][0] - a[1][1] * a[2][2] * a[3][0] - a[1][2] * a[2][0] * a[3][1] + a[1][0] * a[2][2] * a[3][1] + a[1][1] * a[2][0] * a[3][2] - a[1][0] * a[2][1] * a[3][2];
			RI->data[3][1] = a[0][1] * a[2][2] * a[3][0] - a[0][2] * a[2][1] * a[3][0] + a[0][2] * a[2][0] * a[3][1] - a[0][0] * a[2][2] * a[3][1] - a[0][1] * a[2][0] * a[3][2] + a[0][0] * a[2][1] * a[3][2];
			RI->data[3][2] = a[0][2] * a[1][1] * a[3][0] - a[0][1] * a[1][2] * a[3][0] - a[0][2] * a[1][0] * a[3][1] + a[0][0] * a[1][2] * a[3][1] + a[0][1] * a[1][0] * a[3][2] - a[0][0] * a[1][1] * a[3][2];
			RI->data[3][3] = a[0][1] * a[1][2] * a[2][0] - a[0][2] * a[1][1] * a[2][0] + a[0][2] * a[1][0] * a[2][1] - a[0][0] * a[1][2] * a[2][1] - a[0][1] * a[1][0] * a[2][2] + a[0][0] * a[1][1] * a[2][2];

			_mul_m_s_inplace(RI, ((entry_t) 1.0) / A->det);

			return RI;
		default:
			RI = DeepCopy(A);
			I = Identity(A->rowsCount);

			join(RI, I);

			toReducedRowEchelonForm(RI);

			for (size_t i = 0; i < A->rowsCount; i++) {
				for (size_t j = 0; j < A->colsCount; j++) {
					I->data[i][j] = RI->data[i][j + A->colsCount];
				}
			}

			freeMat$(RI);

			return I;
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

 \brief	Matrix multiplication using naive, but cache-friendly method)
		ijk - 1.25 cache misses per iteration ikj - 0.5 cache misses (row-wise)

 \date	24-May-14

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	Matrix product of A & B.
 */
Mat MatMul_naive (Mat A, Mat B) {
	Assert$(A->colsCount == B->rowsCount, "Number of columns in A must be equal to number of rows in B.");

	Mat C = AllocMat(A->rowsCount, B->colsCount);
	entry_t **a = A->data;
	entry_t **b = B->data;
	entry_t **c = C->data;

	if (IsSquare$(A) && IsDimsEqual(A, B)) {
		switch (A->rowsCount) {
			case 1:
				c[0][0] = a[0][0] * b[0][0];

				return C;
			case 2:
				c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
				c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];

				c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
				c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];

				return C;
			case 3:
				c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
				c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
				c[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];

				c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
				c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
				c[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];

				c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
				c[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
				c[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];

				return C;
			default:
				goto generalized;
		}
	} else {
		goto generalized;
	}

generalized:
	fill_zeroes(C);
	for (size_t i = 0; i < A->rowsCount; ++i) {
		for (size_t k = 0; k < A->colsCount; ++k) {
			entry_t s = a[i][k];
			for (size_t j = 0; j < B->colsCount; ++j) {
				c[i][j] += s * b[k][j];
			}
		}
	}

	return C;
}

/**
 \fn	Mat MatMul_naive_recursive (Mat A, Mat B)

 \brief	Recursive implementation of naive matrix multiplication algorithm.
 		Only for square matrices with size=2^n. Not memory efficient!
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
	Assert$(IsSquare$(A) && IsDimsEqual(A, B), "Matrices A and B should be square and equal by sizes.");

	Mat C = NULL;
	size_t size = A->rowsCount;
	entry_t **a = A->data;
	entry_t **b = B->data;

	switch (size) {
		case 1:
			C = AllocMat(1, 1);
			C->data[0][0] = a[0][0] * b[0][0];

			return C;
		case 2:
			C = AllocMat(2, 2);
			C->data[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
			C->data[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];
			C->data[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
			C->data[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];

			return C;
		case 3:
			C = AllocMat(3, 3);

			C->data[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
			C->data[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
			C->data[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];

			C->data[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
			C->data[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
			C->data[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];

			C->data[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
			C->data[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
			C->data[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];

			return C;
		default:
			if (A->rowsCount < RECURSION_LIMIT) {
				C = MatMul_naive(A, B);

				return C;
			} else {
				C = AllocMat(size, size);

				Mat tmp_1 = NULL;

				size /= 2;

				Mat A11 = AllocMat(size, size);
				Mat A12 = AllocMat(size, size);
				Mat A21 = AllocMat(size, size);
				Mat A22 = AllocMat(size, size);

				Mat B11 = AllocMat(size, size);
				Mat B12 = AllocMat(size, size);
				Mat B21 = AllocMat(size, size);
				Mat B22 = AllocMat(size, size);

				// Split
				for (size_t i = 0; i < size; i++) {
					for (size_t j = 0; j < size; j++) {
						A11->data[i][j] = a[i][j];
						A12->data[i][j] = a[i][j + size];
						A21->data[i][j] = a[i + size][j];
						A22->data[i][j] = a[i + size][j + size];

						B11->data[i][j] = b[i][j];
						B12->data[i][j] = b[i][j + size];
						B21->data[i][j] = b[i + size][j];
						B22->data[i][j] = b[i + size][j + size];
					}
				}

				//TODO: C22 can be used as temporary, so no tmp_1 needed!)
				tmp_1 = MatMul_naive_recursive(A12, B21);
				Mat C11 = MatMul_naive_recursive(A11, B11);
				_add_m_m_inplace(C11, tmp_1);
				freeMat$(tmp_1);

				tmp_1 = MatMul_naive_recursive(A12, B22);
				Mat C12 = MatMul_naive_recursive(A11, B12);
				_add_m_m_inplace(C12, tmp_1);
				freeMat$(tmp_1);

				tmp_1 = MatMul_naive_recursive(A22, B21);
				Mat C21 = MatMul_naive_recursive(A21, B11);
				_add_m_m_inplace(C21, tmp_1);
				freeMat$(tmp_1);

				tmp_1 = MatMul_naive_recursive(A22, B22);
				Mat C22 = MatMul_naive_recursive(A21, B12);
				_add_m_m_inplace(C22, tmp_1);
				freeMat$(tmp_1);

				// Join
				for (size_t i = 0; i < size; i++) {
					for (size_t j = 0; j < size; j++) {
						C->data[i][j] = C11->data[i][j];
						C->data[i][j + size] = C12->data[i][j];
						C->data[i + size][j] = C21->data[i][j];
						C->data[i + size][j + size] = C22->data[i][j];
					}
				}

				freeMats(A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22, NULL);

				return C;
			}
	}
}

/**
 \fn	void matMul_inplace(Mat *A, Mat B)

 \brief	In-place matrix multiplication (multiplicand is replaced by result).

 \date	24-May-14

 \param	A	The Mat A to process (will be replaced by product of A &amp; B).
 \param	B	The Mat B to process.
 */
void matMul_inplace (Mat *A, Mat B) {
	Assert$((*A)->colsCount == B->rowsCount, "Number of columns of A must be equal to number of rows of B.");

	Mat P = MatMul$(*A, B);
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
	Assert$(IsSquare$(A), "Cannot raise power of non-square matrix.");

	size_t c = 0;
	Mat R;

	switch (pow) {
		case 0:
			return Identity(A->rowsCount);
		case 1:
			R = DeepCopy(A);

			return R;
		case 2:
			R = DeepCopy(A);

			matMul_inplace(&R, A);

			return R;
		default:
			R = DeepCopy(A);

			do {
				matMul_inplace(&R, A);
				c++;
			} while (c < pow - 1);

			return R;
	}
}
#pragma endregion "Matrix multiplication"


#pragma region "Norms, etc."
/**
 \fn	size_t Rank (Mat A)

 \brief	Computes Rank of Matrix A in RREF form.

 \date	17-May-14

 \param	A	The Mat to process.

 \return	Rank value.
 */
size_t Rank (Mat RREF) {
	entry_t **a = RREF->data;
	size_t rank = 0;

	for (size_t i = 0; i < RREF->rowsCount; i++) {
		for (size_t j = 0; j < RREF->colsCount; j++) {
			if (isnotzero(a[i][j])) {
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
entry_t Trace (Mat A) {
	entry_t **a = A->data;
	entry_t tr = 0.0;

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
entry_t OneNorm (Mat A) {
	entry_t norm = 0.0;

	for (size_t i = 0; i < A->colsCount; i++) {
        entry_t sum = 0.0;

		for (size_t j = 0; j < A->rowsCount; j++) {
			sum += abs(A->data[j][i]);
		}

		norm = max(norm, sum);
	}

	return norm;
}

entry_t TwoNorm (Mat A) {
	Assert$(false, "Not Yet Implemented."); //TODO:

    return 0.0;
}

/**
 \fn	double InfinityNorm (Mat A)

 \brief	Infinity-norm of A (max row sum).

 \param	A	The Mat to process.

 \return	Inf-norm.
 */
entry_t InfinityNorm (Mat A) {
    entry_t norm = 0.0;

	for (size_t i = 0; i < A->rowsCount; i++) {
        entry_t sum = 0.0;

		for (size_t j = 0; j < A->colsCount; j++) {
			sum += abs(A->data[i][j]);
		}

		norm = max(norm, sum);
	}

	return norm;
}

/**
 \fn	double EuclideanNorm (Mat A)

 \brief	Euclidean norm of A.
 		Aka Frobenius norm.

 \param	A	The Mat to process.

 \return	Euclidean norm.
 */
entry_t EuclideanNorm (Mat A) {
    entry_t sum = 0.0;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			sum += square(A->data[i][j]); //TODO: check for overflow
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
entry_t ConditionNumber (Mat A) {
	Mat Ai = DeepCopy(A);

	toInverse(&Ai);

	if (Ai != NULL) {
		entry_t c = InfinityNorm(A) * InfinityNorm(Ai);

		freeMat$(Ai);

		return c;
	} else {
		return NAN;
	}
}

entry_t DiagProd (Mat A) {
	entry_t prod = 1.0;

	for (size_t i = 0; i < A->rowsCount; ++i) {
		prod *= A->data[i][i];
	}

	return prod;
}
#pragma endregion "Norms, etc."


#pragma region "Kronecker"
/**
 \fn	Mat KroneckerProd (Mat A, Mat B)

 \brief	Kronecker product of A and B.

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	A (⨯) B.
 */
Mat KroneckerProd (Mat A, Mat B) {
	Mat K = AllocMat(A->rowsCount * B->rowsCount, A->colsCount * B->colsCount);
	entry_t **a = A->data;
	entry_t **b = B->data;
	entry_t **k = K->data;

	for (size_t i = 0; i < K->rowsCount; i++) {
		for (size_t j = 0; j < K->colsCount; j++) {
			k[i][j] = (a[i / A->rowsCount][j / A->colsCount] * b[i % B->rowsCount][j % B->colsCount]);
		}
	}

	return K;
}

/**
 \fn	Mat KroneckerSum (Mat A, Mat B)

 \brief	Kronecker sum of A and B.

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	A (+) B.
 */
Mat KroneckerSum (Mat A, Mat B) {
	Assert$(IsSquare$(A) && IsSquare$(B), "Matrices A and B should be square.");

	Mat Ib = Identity(B->rowsCount);
	Mat Ia = Identity(A->rowsCount);
	Mat AI = KroneckerProd(A, Ib);
	Mat IB = KroneckerProd(Ia, B);
	Mat KS = AllocMat(A->rowsCount * B->rowsCount, A->colsCount * B->colsCount);

	_add_m_m_m(AI, IB, KS);

	freeMat$(IB);
	freeMat$(Ia);
	freeMat$(AI);
	freeMat$(Ib);

	return KS;
}
#pragma endregion "Kronecker"


#pragma region "Strassen"
/**
 \fn	size_t _adjustSize (size_t Size)

 \brief	Fix size.

 \param	Size	The size.

 \return		Fixed Size value.
 */
size_t _adjustSize (size_t size) {
	if (!(Check$(ispowerof2_i(size), "Matrix size is not a power of 2."))) {
		return (size_t) ((int64_t) (1) << (int64_t) (ceil(log2(size)))); //-V113 //TODO: get rid of FP operations here
	} else {
		return size;
	}
}

/**
 \fn	Mat MatMul_Strassen (Mat A, Mat B)

 \brief	Finds product of A & B using Strassen's algorithm.
 Source matrices A, B & its product C will be divided into 4 square blocks (submatrices),
 and this algorithm repeated recursively until blocks become numbers.
 (or when size of blocks reaches some threshold value when naive algorithm will be used).
 Complexity: O(n^2.8)
 TODO: convert input matrices into 'recursion-friendly' form (e.g 'matrix of matrices')
		M1 := (A11 + A22)(B11 + B22)
		M2 := (A21 + A22)B11
		M3 := A11(B12 − B22)
		M4 := A22(B21 − B11)
		M5 := (A11 + A12)B22
		M6 := (A21 − A11)(B11 + B12)
		M7 := (A12 − A22)(B21 + B22)
		C11 = M1 + M4 − M5 + M7
		C12 = M3 + M5
		C21 = M2 + M4
		C22 = M1 − M2 + M3 + M6.

 \param	A	The Mat A to process.
 \param	B	The Mat B to process.

 \return	Product of A & B.
 */
Mat MatMul_Strassen (Mat A, Mat B) {
	Mat C = NULL;
	size_t size = A->rowsCount;
	entry_t **a = A->data;
	entry_t **b = B->data;

	switch (size) {
		case 1:
			C = AllocMat(size, size);

			C->data[0][0] = a[0][0] * B->data[0][0];

			return C;
		case 2:
			C = AllocMat(size, size);

			C->data[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
			C->data[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];
			C->data[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
			C->data[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];

            return C;
		case 3:
			C = AllocMat(size, size);

			C->data[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
			C->data[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
			C->data[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];

			C->data[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
			C->data[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
			C->data[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];

			C->data[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
			C->data[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
			C->data[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];

			return C;
		default:
			if (size < RECURSION_LIMIT) {
				C = MatMul_naive(A, B);

				return C;
			} else {
				C = AllocMat(size, size);

				size /= 2;

				Mat A11 = AllocMat(size, size);
				Mat A12 = AllocMat(size, size);
				Mat A21 = AllocMat(size, size);
				Mat A22 = AllocMat(size, size);

				Mat B11 = AllocMat(size, size);
				Mat B12 = AllocMat(size, size);
				Mat B21 = AllocMat(size, size);
				Mat B22 = AllocMat(size, size);

				// Split
				for (size_t i = 0; i < size; i++) {
					for (size_t j = 0; j < size; j++) {
						A11->data[i][j] = a[i][j];
						A12->data[i][j] = a[i][j + size];
						A21->data[i][j] = a[i + size][j];
						A22->data[i][j] = a[i + size][j + size];

						B11->data[i][j] = b[i][j];
						B12->data[i][j] = b[i][j + size];
						B21->data[i][j] = b[i + size][j];
						B22->data[i][j] = b[i + size][j + size];
					}
				}

				Mat A_tmp = AllocMat(size, size);
				Mat B_tmp = AllocMat(size, size);

				_add_m_m_m(A11, A22, A_tmp); // A11 + A22
				_add_m_m_m(B11, B22, B_tmp); // B11 + B22
				Mat P1 = MatMul_Strassen(A_tmp, B_tmp); // P1 = (A11+A22) * (B11+B22)

				_add_m_m_m(A21, A22, A_tmp); // A21 + A22
				Mat P2 = MatMul_Strassen(A_tmp, B11); // P2 = (A21+A22) * (B11)

				_sub_m_m_m(B12, B22, B_tmp); // B12 - B22
				Mat P3 = MatMul_Strassen(A11, B_tmp); // P3 = (A11) * (B12 - B22)

				_sub_m_m_m(B21, B11, B_tmp); // B21 - B11
				Mat P4 = MatMul_Strassen(A22, B_tmp); // P4 = (A22) * (B21 - B11)

				_add_m_m_m(A11, A12, A_tmp); // A11 + A12
				Mat P5 = MatMul_Strassen(A_tmp, B22); // P5 = (A11 + A12) * (B22)

				_sub_m_m_m(A21, A11, A_tmp); // A21 - A11
				_add_m_m_m(B11, B12, B_tmp); // B11 + B12
				Mat P6 = MatMul_Strassen(A_tmp, B_tmp); // P6 = (A21 - A11) * (B11 + B12)

				_sub_m_m_m(A12, A22, A_tmp); // A12 - A22
				_add_m_m_m(B21, B22, B_tmp); // B21 + B22
				Mat P7 = MatMul_Strassen(A_tmp, B_tmp); // P7 = (A12 - A22) * (B21 + B22)

				Mat C11 = AllocMat(size, size);
				Mat C12 = AllocMat(size, size);
				Mat C21 = AllocMat(size, size);
				Mat C22 = AllocMat(size, size);

				_add_m_m_m(P3, P5, C12); // C12 = P3 + P5
				_add_m_m_m(P2, P4, C21); // C21 = P2 + P4

				_add_m_m_m(P1, P4, A_tmp); // P1 + P4
				_add_m_m_m(A_tmp, P7, B_tmp); // P1 + P4 + P7
				_sub_m_m_m(B_tmp, P5, C11); // C11 = P1 + P4 - P5 + P7

				_add_m_m_m(P1, P3, A_tmp); // P1 + P3
				_add_m_m_m(A_tmp, P6, B_tmp); // P1 + P3 + P6
				_sub_m_m_m(B_tmp, P2, C22); // C22 = P1 + P3 - P2 + P6

				// Join
				for (size_t i = 0; i < size; i++) {
					for (size_t j = 0; j < size; j++) {
						C->data[i][j] = C11->data[i][j];
						C->data[i][j + size] = C12->data[i][j];
						C->data[i + size][j] = C21->data[i][j];
						C->data[i + size][j + size] = C22->data[i][j];
					}
				}

				freeMats(
						A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22,
						A_tmp, B_tmp, P1, P2, P3, P4, P5, P6, P7, NULL
				);

				return C;
			}
	}
}

/**
 \fn	Mat MatMul_Strassen_optimized (Mat A, Mat B)

 \brief	P1 = (A01+A10)(B10+B01)
		P2 = (A10+A11)B10
		P3 = A01(B11-B01)
		P4 = A10(B00-B10)
		P5 = (A01+A00)B01
		P6 = (A11-A01)(B10+B11)
		P7 = (A00-A10)(B00+B01)
		C00 = P1+P4-P5+P7
		C01 = P3+P5
		C10 = P2+P4
		C11 = P1-P2+P3+P6.

 \param	A	The Mat to process.
 \param	B	The Mat to process.

 \return	A . B.
 */
Mat MatMul_Strassen_optimized (Mat A, Mat B) {
	return MatMul_Strassen(A, B); //TODO:
}
#pragma endregion "Strassen"
