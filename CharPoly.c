#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "CharPoly.h"
#include "Matrix.h"
#include "Extra.h"
#include "EigenDcmp_n.h"
#include "Const.h"
#include "Config.h"
#include "SpinningIndicator.h"


/**
 \fn	int64_t *GetCharPolyCoeffs (Mat A)

 \brief	Computes operator' characteristic polynomial coeffs by the Method of Direct Expansion.

 \date	14-May-14

 \param	A	The matrix to process.

 \return	null if it fails, else the * to array of polynomial coeffs.
 */
int64_t *GetCharPolyCoeffs (Mat A) {
	int64_t *coeffs = AllocVec_i(A->rowsCount);

	for (size_t i = 0; i < A->rowsCount; i++) {
		coeffs[i] = GetPrincipalMinorsSum(A, i + 1);
		if (!(i % 2)) {
			coeffs[i] *= -1; //-V127
		}
	}

	return coeffs;
}

/**
 \fn	void printCharacteristicEquation (int64_t *c, size_t s, FILE *file)

 \brief	Prints characteristic equation.

 \date	14-May-14

 \param [in]	c		If non-null, the * to array of polynomial coeffs to read from.
 \param	s				The array size.
 \param [out]	file	If non-null, the file.
 */
void printCharacteristicEquation (int64_t *c, size_t s, FILE *file) {
	Assert$(c != NULL, "");
	Assert$(file != NULL, "");

	for (size_t i = 0; i < s; i++) {
		fprintf(file, "x^" FMT_SIZET, s - i);
		fprintf(file, " + (%lld)", c[i]);
	}
	fprintf(file, " = 0\n");

	return;
}


#pragma region "Evaluation & root finding"

/**
 \fn	long double EvalPolyAt (long double x, int64_t *c, size_t s)

 \brief	Evaluates polynomial value at point X using Horner's method.  

 \date	14-May-14

 \param	x		 	The x coordinate.
 \param [in]	c	If non-null, the * to array of coeffs.
 \param	s		 	The array size.

 \return	Evaluated value.
 */
long double EvalPolyAt (long double x, int64_t *c, size_t s) {
	long double res = 0.0;

	res = 1.0 + res*x;
	for (size_t i = 0; i < s; i++)	{
		res = c[i] + ((res) * (x));
	}

	return res;
}

/**
 \fn	long double FindRoot (int64_t *c, size_t s, long double a, long double b, size_t maxIters)

 \brief	Searches for the equation's first root in range [a, b] using Regula Falsi method.

 \date	15-May-14

 \param [in]	c	If non-null, the * to array of polynomial coefficients.
 \param	s		 	The coeffs array size.
 \param	a		 	Start position.
 \param	b		 	End position.
 \param	maxIters 	Number of maximum iterations.

 \return	The found root value.
 */
long double FindRoot (int64_t *c, size_t s, long double a, long double b, size_t maxIters) {
	long double falseRoot = 0.0;
    long double root = 0.0;
    int x = 0;
	long double fA = EvalPolyAt(a, c, s);
	long double fB = EvalPolyAt(b, c, s);

	for (size_t i = 0; i < maxIters; i++) {
		root = (fA*b - fB*a) / (fA - fB);
		if (equals_ld(a, b)) {
			//printf("# %f %f \n", a, b);
			break;
		}
		falseRoot = EvalPolyAt(root, c, s);

		if (falseRoot * fB > 0) {
			b = root;
			fB = falseRoot;
			if (x == -1) {
				fA /= 2;
			}
			x = -1;
		} else {
			if (fA * falseRoot > 0) {
				a = root;
				fA = falseRoot;
				if (x == +1) {
					fB /= 2;
				}
				x = +1;
			} else {
				break;
			}
		}
	}
	return root;
}

/**
 \fn	long double *GetPolynomialRoots (int64_t *c, size_t s, size_t *rootsCount)

 \brief	Searches for all the roots of given equation.

 \date	15-May-14

 \param [in]	c		  	If non-null, the * to coeffs array.
 \param	s				  	The array size to process.
 \param [out]	rootsCount	If non-null, number of roots.

 \return	null if it fails, else the * to array of polynomial roots.
 */
long double *GetPolynomialRoots (int64_t *c, size_t size, size_t *rootsCount) {
	long double a = START, b = a + DELTA, e = END;
	long double root = 0.0;
	long double *roots = AllocVec_ld(size);
	bool isValid = false;
	size_t i = 0;

	for (i = 0; i < size; i++) {
		if (c[i] != 0) {
			isValid = true;
			break;
		}
	}
	if (!isValid) {
		*rootsCount = 0;
		roots[0] = 0;
		printf("No roots exists. Only the trivial solution (x=0) is available.\n");
		return roots;
	}

	while (b < e) {
		spinActivityIndicator();
		root = round(FindRoot(c, size, a, b, MAXITERS)); //To round or not to round?

		if (fabsl(EvalPolyAt(root, c, size)) <= EPS) {
			if (!(exists_d(roots, i+1, root))) {
				roots[i++] = (root);
			}
		}
		a+=DELTA; b+=DELTA;
	}
	clearActivityIndicator();

	*rootsCount = i;
	if (i == 0) {
		printf("No eigenvalues found,\nbecause function failed to converge because of no integer roots exists\nor invalid roots search interval (you can change it in <Const.h>).\n");
		free(roots); roots = NULL;
		return NULL;
	}

	return roots;
}
#pragma endregion "Evaluation & root finding"
