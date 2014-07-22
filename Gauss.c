#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <limits.h>

#include "Matrix.h"
#include "MatrixOperations.h"
#include "Extra.h"
#include "Gauss.h"


/**
 \fn	double GaussianDeterminant (dMat A)	 
 \brief	Calculates matrix determinant by Gauss' method for n>3
 and with rule of Sarrus for n<=3.			 
 \date	13-May-14							 
 \param	A   	The double-valued matrix to process.
 \return		Determinant value.
 */
//TODO: long double?
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
		default:
			T = DeepCopy(A);
			Assert(T != NULL, "Cannot create copy...");

			toRowEchelonForm(T);
			
			if (!(T->isSingular)) {
				double **t = T->a;
				for (size_t i = 0; i < A->rowsCount; i++) {
					det *= t[i][i];
				}
				det *= T->permutationSign;
			} else {
				det = 0.0;
				A->isSingular = true;
			}			
			FreeMat(T);
			break;
	}
	A->det = det;

	return det;
}

void toRowEchelonForm (Mat A) {
	double **a = A->a;

	//optimized? with immediate rows swapping
	for (size_t k = 0; k < A->rowsCount; k++) {
		for (size_t i = k + 1; i < A->rowsCount; i++) {
				if (fabs(a[i][k]) > EPS) {
					for (size_t j = 0; j < A->colsCount; j++) {
						swap_d(a[k][j], a[i][j]);
					}
					A->permutationSign *= -1; //-V127
					break;
				} 
		}
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			if (fabs(a[k][k]) > EPS) {
				double factor = a[i][k] / a[k][k];
				for (size_t j = k + 0; j < A->colsCount; j++) { //no need to set lower<diag to 0
					a[i][j] -= factor * a[k][j];
				}
			} else {
				A->isSingular = true;
				//Check(0, "toRowEchelonForm: singular mat.");
				return;
			}
		}
	}

	return;
}
void toRowEchelonForm_r (Mat A) {
	double **a = A->a;

	//reference implementation of pivoting algorithm
	for (size_t k = 0; k < A->rowsCount; k++) {
		size_t pivot = k;
		//double Max = a[k][k];
		for (size_t i = k + 1; i < A->rowsCount; i++) {
			//double abs = fabs(a[i][k]);
			//if (abs > Max)
			if (fabs(a[i][k]) > fabs(a[pivot][k])) {
				pivot = i;
				//Max = abs;
			}
		}
		if (fabs(a[pivot][k]) <= EPS) {
			A->isSingular = true;
			//Check(0, "toRowEchelonForm_r: Singular matrix.");
			return;
		}
		if (pivot != k) {
			for (size_t i = 0; i < A->colsCount; i++) {
				swap_d(a[k][i], a[pivot][i]);
			}
			A->permutationSign *= -1; //-V127
		}
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
 \fn	void toReducedRowEchelonForm (dMat A)
 \brief	(In-place) Transforms matrix into reduced row echelon form
 (aka row canonical form).
 \date	13-May-14		  
 \param	A   	The double-valued matrix to process.
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

Mat Solve_gaussjordan (Mat A, Mat B) {
	Assert(B->colsCount == 1, "Use LUSolver for this.");
	Assert(A->rowsCount == B->rowsCount, "Number of equations doesn't equal to number of unknown variables.");

	Mat X = AllocMat(A->rowsCount, 1);
	Mat AU = DeepCopy(A);
	Assert(AU != NULL, "Cannot create copy...");

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

Mat Solve_gauss (Mat A, Mat B) {
	Assert(B->colsCount == 1, "");
	Assert(A->rowsCount == B->rowsCount, "Number of equations doesn't equal to number of unknown variables.");

	Mat X = AllocMat(B->rowsCount, B->colsCount);
	Mat AU = DeepCopy(A);
	Assert(AU != NULL, "Cannot create copy...");

	
	concat(AU, B);
	toRowEchelonForm_r(AU);

	//print(AU);

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


//------------------------------------------------------------------------------

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
		Assert(Copy != NULL, "Cannot create copy...");

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
	Assert(copyA != NULL, "Cannot create copy...");	

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
