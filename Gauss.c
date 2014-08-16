#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <limits.h>

#include "Matrix.h"
#include "MatrixOperations.h"
#include "Extra.h"
#include "Gauss.h"


#pragma region "Determinant computation"

/**
 \fn	double Det_gauss (Mat A)

 \brief	Calculates matrix determinant by Gauss' method for n&gt;3 and with rule of Sarrus for
		n&lt;=3. 
		TODO: long double?

 \date	13-May-14

 \param	A	The Mat to process.

 \return	Determinant value.
 */
double Det_gauss (Mat A) {
	double **a = A->a;
	Mat T = NULL;
	double det = 1.0;

	switch (A->rowsCount) {
		case 1:
			det = (a[0][0]);
			if (fabs(det) <= EPS) { A->isSingular = true; }
			break;
		case 2:
			det = (((a[0][0] * a[1][1]) - (a[0][1] * a[1][0])));
			if (fabs(det) <= EPS) { A->isSingular = true; }
			break;
		case 3:
			det = (
				(a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])) -
				(a[0][1] * (a[2][2] * a[1][0] - a[1][2] * a[2][0])) +
				(a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
			);
			if (fabs(det) <= EPS) { A->isSingular = true; }
			break;
		case 4:
			det = (
				a[0][0] * (
					a[1][1] * (a[2][2] * a[3][3] - a[2][3] * a[3][2]) -
					a[2][1] * (a[1][2] * a[3][3] - a[1][3] * a[3][2]) +
					a[3][1] * (a[1][2] * a[2][3] - a[1][3] * a[2][2])
				) -
				a[1][0] * (
					a[0][1] * (a[2][2] * a[3][3] - a[2][3] * a[3][2]) -
					a[2][1] * (a[0][2] * a[3][3] - a[0][3] * a[3][2]) +
					a[3][1] * (a[0][2] * a[2][3] - a[0][3] * a[2][2])
				) +
				a[2][0] * (
					a[0][1] * (a[1][2] * a[3][3] - a[1][3] * a[3][2]) -
					a[1][1] * (a[0][2] * a[3][3] - a[0][3] * a[3][2]) +
					a[3][1] * (a[0][2] * a[1][3] - a[0][3] * a[1][2])
				) -
				a[3][0] * (
					a[0][1] * (a[1][2] * a[2][3] - a[1][3] * a[2][2]) -
					a[1][1] * (a[0][2] * a[2][3] - a[0][3] * a[2][2]) +
					a[2][1] * (a[0][2] * a[1][3] - a[0][3] * a[1][2])
				)
			);
			if (fabs(det) <= EPS) { A->isSingular = true; }
			break;
		default:
			T = DeepCopy(A);

			toRowEchelonForm(T);
			
			if (!(T->isSingular)) {
				for (size_t i = 0; i < A->rowsCount; i++) {
					det *= T->a[i][i];
				}
				det *= T->permutationSign;
			} else {
				det = 0.0;
				A->isSingular = true;
			}			
			FreeMat(T);
			break;
	}
	//if (fabs(det) <= EPS) { A->isSingular = true; }
	A->det = det;

	return det;
}

/**
 \fn	double Det_bareiss (Mat A)

 \brief	Computes Matrix determinant by Baeriss' algorithm.
		The Bareiss Algorithm is fraction-free method for determinant
		computation.
		However, it can also be thought of as a sophisticated form of row reduction.
		Note that the divisions computed at any step are exact; thus Bareiss’ Algorithm is
		indeed fraction-free. Entry a[n][n] is the determinant of A (after Bareiss and pivoting steps).

 \param	A	The Mat to process.

 \return	Determinant.
 */
double Det_bareiss (Mat A) {
	Mat T = DeepCopy(A);

	// Pivotize
	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		for (size_t i = k; i < A->rowsCount; i++) {
			if (fabs(T->a[i][k]) > fabs(T->a[pivot][k])) {
				pivot = i;
			}
		}
		if (pivot != k) {
			T->permutationSign *= -1; //-V127
			for (size_t j = 0; j < A->colsCount; j++) {
				swap_d(T->a[k][j], T->a[pivot][j]);
			}
		}
	}
	
	// Bareiss algorithm main step
	for (size_t i = 0; i < T->rowsCount - 1; i++) {
		//Assert(T->a[i][i] > EPS, "Singularity...");
		for (size_t j = i + 1; j < T->rowsCount; j++)
			for (size_t k = i + 1; k < T->rowsCount; k++) {
			T->a[j][k] = T->a[j][k]*T->a[i][i] - T->a[j][i]*T->a[i][k];
			if (i) {
				T->a[j][k] /= T->a[i-1][i-1];
			}
		}
	}

	double det = T->a[T->rowsCount - 1][T->rowsCount - 1];
	FreeMat(T);

	return det;
}
#pragma endregion "Determinant computation"


#pragma region "Transforming routines"

/**
 \fn	void toRowEchelonForm (Mat A)

 \brief	Transforms matrix A into a row echelon form.

 \param	A	The Mat to process.
 */
void toRowEchelonForm (Mat A) {
	double **a = A->a;

	//optimized? with immediate rows swapping
	for (size_t k = 0; k < A->rowsCount; k++) {
		// Pivotize
		for (size_t i = k + 1; i < A->rowsCount; i++) {
				if (fabs(a[i][k]) > EPS) {
					for (size_t j = 0; j < A->colsCount; j++) {
						swap_d(a[k][j], a[i][j]);
					}
					A->permutationSign *= -1; //-V127
				}
		}
		// Eliminate
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (fabs(a[k][k]) > EPS) {
				double factor = a[i][k] / a[k][k];
				for (size_t j = k + 0; j < A->colsCount; j++) { //no need to set lower<diag to 0
					a[i][j] -= factor * a[k][j];
				}
			} else {
				A->isSingular = true;
				Check(0, "toRowEchelonForm: singular mat.");
				return;
			}
		}
	}

	return;
}

/**
 \fn	void toRowEchelonForm_r (Mat A)

 \brief	Transforms A to a row echelon form ('reference' implementation).

 \param	A	The Mat to process.
 */
void toRowEchelonForm_r (Mat A) {
	double **a = A->a;

	//reference implementation of pivoting algorithm
	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		// Find pivot
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (fabs(a[i][k]) > fabs(a[pivot][k])) {
				pivot = i;
			}
		}
		if (fabs(a[pivot][k]) <= EPS) {
			A->isSingular = true;
			Check(0, "toRowEchelonForm_r: Singular matrix.");
			return;
		}
		// Swap rows
		if (pivot != k) {
			for (size_t i = 0; i < A->colsCount; i++) {
				swap_d(a[k][i], a[pivot][i]);
			}
			A->permutationSign *= -1; //-V127
		}
		// Eliminate
		//if (fabs(a[k][k]) > EPS) {
			for (size_t i = k + 1; i < A->rowsCount; i++) {
				//a[i][k] /= a[k][k];
				double f = a[i][k] / a[k][k];
				for (size_t j = k + 0; j < A->colsCount; j++) {
					//a[i][j] -= a[i][k] * a[k][j];
					a[i][j] -= f * a[k][j];
					//a[i][j] -= a[k][j] * (a[i][k] / a[k][k]);
				}
				//a[i][k] = 0.0;
			}
		//}
	}

	return;
}

/**
 \fn	void toReducedRowEchelonForm (Mat A)

 \brief	(In-place) Transforms matrix into reduced row echelon form (aka row canonical form).
		The reduced row echelon form of A is unique, the pivot positions are
		uniquely determined and do not depend on whether or not row interchanges
		are performed in the reduction process.

 \date	13-May-14

 \param	A	The double-valued matrix to process.
 */
void toReducedRowEchelonForm (Mat A) {
	double **a = A->a;
	size_t pivotCol = 0;
	size_t rowsCount = A->rowsCount, colsCount = A->colsCount;
	size_t i, j, k, r;

	for (r = 0; r < rowsCount; r++) {
		if (colsCount < pivotCol + 1) {
			break;
		}
		i = r;
		while (fabs(a[i][pivotCol]) <= EPS) {
			i++;
			if (i == rowsCount) {
				i = r;
				pivotCol++;
				if (colsCount == pivotCol) {
					pivotCol--;
					break;
				}
			}
		}
		for (j = 0; j < colsCount; j++) {
			swap_d(a[i][j], a[r][j]);
		}
		if (fabs(a[r][pivotCol]) > EPS) {
			double divisor = a[r][pivotCol];
			for (j = 0; j < colsCount; j++) {
				a[r][j] /= divisor;
			}
		}
		for (j = 0; j < rowsCount; j++) {
			if (j != r)	{
				double sub = a[j][pivotCol];
				for (k = 0; k < colsCount; k++) {
					a[j][k] -= (sub * a[r][k]);
				}
			}
		}
		pivotCol++;
	}

	return;
}
#pragma endregion "Transforming routines"


#pragma region "Solving routines"

/**
 \fn	Mat Solve_gaussjordan (Mat A, Mat B)

 \brief	Solves system of linear equations using Gauss-Jordan method.

 \param	A	Coeffs matrix.
 \param	B	Right hand side.

 \return	Solution as column-vector.
 */
Mat Solve_gaussjordan (Mat A, Mat B) {
	Assert(B->colsCount == 1, "");
	Assert(A->rowsCount == B->rowsCount, "Number of equations doesn't equal to number of unknown variables.");

	Mat X = AllocMat(A->rowsCount, 1);
	Mat AU = DeepCopy(A);

	concat(AU, B);
	double **au = AU->a;
	double **x = X->a;

	toReducedRowEchelonForm(AU);

	for (size_t i = 0; i < X->rowsCount; i++) {
		x[i][0] = au[i][AU->colsCount - 1];
	}

	FreeMat(AU);

	return X;
}

/**
 \fn	Mat Solve_gauss (Mat A, Mat B)

 \brief	Solves system of linear equations using Gauss elimination.

 \param	A	Coeffs matrix.
 \param	B	Right hand side.

 \return	Solution as column-vector.
 */
Mat Solve_gauss (Mat A, Mat B) {
	Assert(B->colsCount == 1, "");
	Assert(A->rowsCount == B->rowsCount, "Number of equations doesn't equal to number of unknown variables.");

	Mat X = AllocMat(B->rowsCount, B->colsCount);
	Mat AU = DeepCopy(A);
	
	concat(AU, B);

	//Forward step (elimination with row pivoting)
	toRowEchelonForm_r(AU);

	//Back-substitution
	for (ptrdiff_t i = AU->rowsCount - 1; i >= 0; i--) {
		X->a[i][0] = AU->a[i][AU->colsCount-1];
		for (size_t j = i + 1; j < AU->rowsCount; j++) {
			X->a[i][0] -= AU->a[i][j] * X->a[j][0];
		}
		X->a[i][0] /= AU->a[i][i];
	}

	FreeMat(AU);

	return X;
}
#pragma endregion "Solving routines"


//-----------------------------somewhat outdated, but working...----------------
#pragma region old
/**
 \fn	void simpleSolver (double **a, size_t Size, double *x)
 \brief	Simple back-substitution routine for undeterminedSolver.
 \date	22-May-14												
 \param [in]	a		If non-null, the double ** to matrix in rref form.
 \param	Size		 	The matrix size.
 \param [out]	x		* to array to write the solution to.
 */
void simpleSolver_h (double **a, size_t Size, double *x) {
	ptrdiff_t i;
	size_t j;

	for (i = Size - 1; i >= 0; i--) {
		for (j = i + 1; j < Size; j++) {
			x[i] -= a[i][j];
		}
		if ((fabs(x[i]) > EPS) && (fabs(a[i][i]) > EPS)) {
			x[i] /= a[i][i];
		}
	}

	return;
}

/**
 \fn	void undeterminedSolver (dMat RREF, dMat A, dMat R)
 \brief	Solves undetermined system of linear equations.
 Gives fundamental system of solution vectors as rows of matrix.
 \date	22-May-14												
 \param	[in]	RREF	The matrix A in rref form.
 \param	[in]	A		The matrix A.
 \param	[out]	R   	The matrix with solution vectors.
 */
void undeterminedSolver_h (Mat RREF, Mat A, Mat R) {
	double **rref = RREF->a;
	size_t *f = uAllocVec(A->rowsCount);
	size_t i, j, c = 0, row = 0;

	for (i = 0; i < RREF->rowsCount; i++) {
		for (j = 0; j < RREF->colsCount; j++) {
			if (fabs(rref[i][j]) > EPS) {
				f[c++] = i;
				break;
			}
		}
	}
	c = 0;

	for (i = 0; i < RREF->rowsCount; i++) {
		if (!(exists_u(f, RREF->rowsCount, 0, i))) {
			f[c++] = i;
		}
	}
	j = RREF->rowsCount-1;
	
	do {
		Mat Copy = DeepCopy(A);

		double **copy = Copy->a;
		for (i = 0; i < c-1; i++) {
			copy[(f[i])][(f[i])] += 1;
		}
		toReducedRowEchelonForm(Copy);
		simpleSolver_h(copy, Copy->rowsCount, (R->a[row]));
		for (i = RREF->rowsCount - c; i < RREF->rowsCount; i++) {
			*(R->a[row] + i) = 0.0;
		}
		*(R->a[row] + (j--)) = 1.0;
		FreeMat(Copy);
		row++;
	} while (nextCombination(f, c-1, RREF->rowsCount-1));

	free(f); f = NULL;

	return;
}

/**
\fn	double *GaussianSolve (dMat A)
\brief	Solves homogeneous system of linear equations defined by NxN matrix.
Theorem: Every homogeneous system has either exactly one solution or infinitely many solutions.
If a homogeneous system has more unknowns than equations, then it has infinitely many solutions.
\date	15-May-14								   
\param	[in,out]	A  	The double-valued matrix to process. 
\return	Matrix with solution vectors written as rows.
*/
Mat GaussianSolve_h (Mat A)	{
	double **a = A->a;
	Mat copyA = DeepCopy(A);

	toReducedRowEchelonForm(A);
	A->rank = Rank(A);

	if (A->rank < A->rowsCount) {
		Mat RES = AllocMat(A->rowsCount - A->rank, A->rowsCount);
		undeterminedSolver_h(A, copyA, RES);					  
		FreeMat(copyA);
		FreeMat(A);
		
		return RES;
	} else {
		Mat RES = AllocMat(1, A->rowsCount);
		simpleSolver_h(a, A->rowsCount, RES->a[0]);
		FreeMat(A);
		
		return RES;
	}
}
#pragma endregion old