#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

#include "EigenDcmp_n.h"
#include "Matrix.h"
#include "CharPoly.h"
#include "Gauss.h"
#include "Extra.h"
#include "SpinningIndicator.h"
#include "MatrixOperations.h"


/**
 \fn	void Eigendecomposition (dMat A, size_t n)			  
 \brief	'Naive' implementation of eigenvalue algorithm.
 Constructs eigendecomposition of matrix A,
 using operator' characteristic polynomial roots as its eigenvalues.
 It is not efficient and/or stable method, but it's very
 simple & easy to understand
 So, use it only for reference, not real computations)
 \date	16-May-14											  
 \param	A	The dMat to process.
 */
Mat *EigenDcmp_n (Mat A) {
	size_t i = 0;
	Mat *result = NULL;

	Assert$(A->rowsCount == A->colsCount, "This method is intended only for square matrices.");
	Check$(A->rowsCount <= 18, "This method is very slow on sizes > 18.");

	int64_t *polynomialCoeffs = GetCharPolyCoeffs(A);
	puts(">>>Characteristic equation:");
	printCharacteristicEquation(polynomialCoeffs, A->rowsCount, stdout);
	
	//TODO: error handling
	for (i = 0; i < A->rowsCount; i++) {
		Assert$(LLONG_MAX - llabs(polynomialCoeffs[i]) >= 0, "Overflow.");
	}

	size_t rootsCount = 0;
	size_t *pRootsCount = &rootsCount;
	long double *polynomialRoots = GetPolynomialRoots(polynomialCoeffs, A->rowsCount, pRootsCount);

	if (rootsCount != 0) {
		size_t *k = AllocVec_u(A->rowsCount);
		Mat dCopyA = DeepCopy(A);
		size_t LIVectorsCounter = 0, valuesCount = 0;
		Mat C = AllocMat(A->rowsCount, A->rowsCount);

		for (i = 0; i < rootsCount; i++) {
//			printf(">>>Eigenvalue v%Iu=%.2Lf\n", i, polynomialRoots[i]);
			printf(">>>Eigenvalue v%zu=%.2Lf\n", i, polynomialRoots[i]);

			Mat eigenvectors = GetEigenvectors(dCopyA, polynomialRoots[i]);
			printf(" >>Corresponding eigenvectors for v=%.2Lf (written as rows):\n", polynomialRoots[i]);
			printMat$(eigenvectors);

			LIVectorsCounter += eigenvectors->rowsCount;
			valuesCount++;
			k[i] = eigenvectors->rowsCount;
			fillEigenvectorMatrix(eigenvectors, C, 0, LIVectorsCounter - eigenvectors->rowsCount); //TODO: use concat

			freeMat$(eigenvectors);
		}
		//----------------------eigendecomposition--------------------------------------

		if (LIVectorsCounter != A->rowsCount) {
			puts("Matrix A is not diagonalizable (by this method) so eigendecomposition cannot be performed.\n(too few eigenvectors).");
			return NULL;
		}

		Mat *results = (Mat*) malloc(3 * sizeof(*results));
		Assert$(results != NULL, "Cannot allocate.");

		Mat D = AllocMat(A->rowsCount, A->rowsCount);
		for (i = 0; i < valuesCount; i++) {
			fillEigenvalueMatrix(D, polynomialRoots, i, k[i]);
		}
		
		Mat Ci = DeepCopy(C);
		toInverse(&Ci);

		results[0] = C;
		results[1] = D;
		results[2] = Ci;

		freeMat$(dCopyA);
		free$(k);

		return results;
	} else {
		printf("Matrix A is not diagonalizable so eigendecomposition cannot be performed.\n(no real eigenvalues found).\n");
		return NULL;
	}
	//----------------------free up resources---------------------------------------
	free$(polynomialCoeffs);
	free$(polynomialRoots);

	return result;
}

Mat Spectrum(Mat *evd) {
	Mat Sp = Diag(evd[1]);

	return Sp;
}

double SpectralRadius (Mat Sp) {
	double rad = 0.0;
	for (size_t i = 0; i < Sp->colsCount; i++) {
		rad = max(rad, fabs(Sp->a[0][i]));
	}

	return rad;
}

//void PrintEigendecomposition (EVD res) {
//	Assert$(res != NULL, "Cannot print.");
//
//	printf(" >>Matrix C:\n");
//	printMat$(res->C);
//
//	printf(" >>Matrix D:\n");
//	printMat$(res->D);
//
//	printf(" >>Matrix C^(-1):\n");
//	printMat$(res->Ci);
//
//	return;
//}
//
//void FreeEigendecomposition (EVD res) {
//	Assert$(res != NULL, "Cannot free.");
//
//	freeMat$(res->C);
//	freeMat$(res->D);
//	freeMat$(res->Ci);
//	free(res); res = NULL;
//
//	return;
//}

/**
\fn	dMat EDMatPow (EVD res, size_t n)	
\brief Raise matrix A to power n with eigenvalue decomposition.		
\date	04-Jun-14			   
\param	res		Struct with decomposition results.
\param	n		Power to raise to.
*/
Mat EDMatPow (Mat *evd, size_t n) {
	Mat DP = MatPow(evd[1], n);
	Mat CDP = MatMul(evd[0], DP);
	Mat CDPC_ = MatMul(CDP, evd[2]);

	freeMat$(DP);
	freeMat$(CDP);

	return CDPC_;
}

// A^n = Q*L^n*Q^-1
// A = Q*L*Q^-1
// An n × n matrix A is diagonalizable if and only if the sum of the dimensions of the eigenspaces is n.
void fillEigenvectorMatrix (Mat EV, Mat OUT, size_t s, size_t end) {
	double **e = EV->a;
	double **o = OUT->a;

	for (size_t i = s; i < EV->rowsCount; i++) {
		for (size_t j = s; j < EV->colsCount; j++) {
			o[j][i + end] = e[i][j];
		}
	}

	return;
}

void fillEigenvalueMatrix (Mat OUT, long double *eValues, size_t c, size_t k) {
	double **o = OUT->a;

	for (size_t i = c; i < k+1; i++) {
        o[i][i] = eValues[c]; //TODO: Values of type 'long double' may not fit into the receiver type 'double'
	}

	return;
}

/**
\fn	double *GetEigenvector (dMat A, int eigenValue)	 
\brief	Computes eigenvector of matrix A associated with given eigenvalue.	
\date	15-May-14											  
\param	A		  	The double-valued matrix to process.
\param	eigenvalue	The eigenvalue.				   
\return	null if it fails, else the matrix with eigenvectors.
*/
Mat GetEigenvectors (Mat A, long double eigenvalue) {
	Mat R = DeepCopy(A);
	double **r = R->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		r[i][i] -= eigenvalue;
	}

	return GaussianSolve_h(R);
}

/**
\fn	void StrikeOut (dMat A, size_t d)	
\brief	Strikes out dth row & column (to make principal minor).
TODO: rewrite?			
\date	15-May-14	 
\param	A   	The double-valued matrix to process.
\param	Size	The matrix size.
\param	d   	The D value.
*/
void strikeOut (Mat A, size_t d) {
	double **a = A->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		a[i][d] = 0.0;
		a[d][i] = 0.0;
	}
	a[d][d] = 1.0;

	return;
}

/**
\fn	int GetPrincipalMinorsSum (dMat A, size_t order)  
\brief	Calculates sum of all the principal minors of order K in matrix A.	 
\date	15-May-14															 
\param	A	 	The integer-valued matrix to process.
\param	Size 	The matrix size.
\param	order	The minors order.											 
\return	The principal minors sum.
*/
int64_t GetPrincipalMinorsSum (Mat A, size_t order) {
	int64_t sum = 0;
	size_t *index = AllocVec_u(A->rowsCount);
	fillSequential_u(index, A->rowsCount, 1);

	do {
		Mat M = DeepCopy(A);
		spinActivityIndicator(); 
		for (size_t d = 0; d < A->rowsCount - order; d++) {
			strikeOut(M, index[d] - 1);
		}
		sum += round(Det_gauss(M));
		freeMat$(M);
	} while (nextCombination(index, A->rowsCount - order, A->rowsCount));
	clearActivityIndicator();
	
	free(index); index = NULL;

	return sum;
}
