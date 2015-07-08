#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "Matrix.h"
#include "Config.h"
#include "Const.h"
#include "Extras.h"
#include "Maths.h"
#include "SpinningIndicator.h"


#pragma region "Alloc & dealloc"
/**
 \fn	Mat AllocMat (size_t rowsCount, size_t columnsCount)

 \brief	Allocate memory for Matrix w/ specified number of rows and columns.

 \date	17-May-14

 \param	rowsCount   	The rows count.
 \param	columnsCount	The columns count.

 \return	Pointer to __Mat_struct if allocation was successful, NULL otherwise.
 */
Mat AllocMat (size_t rowsCount, size_t columnsCount) {
	Mat A = NULL;

	Assert$(((rowsCount != 0) && (columnsCount != 0)), "Rows and/or columns count can't be set to zero");

	A = (Mat) malloc(sizeof(*A));
	Assert$(A != NULL, "Cannot allocate memory for __Mat_struct");

	A->a = NULL;
//	A->a = (entry_t**) calloc(rowsCount, sizeof(entry_t*));
	A->a = (entry_t**) malloc(rowsCount * sizeof(entry_t*));
	Assert$(A->a != NULL, "Cannot allocate memory for row pointers");

	for (size_t i = 0; i < rowsCount; i++) {
		A->a[i] = NULL;
//		A->a[i] = (entry_t*) calloc(columnsCount, sizeof(entry_t));
		A->a[i] = (entry_t*) malloc(columnsCount * sizeof(entry_t));
		Assert$(A->a[i] != NULL, "Cannot allocate memory for row");
	}

	A->rowsCount = rowsCount;
	A->colsCount = columnsCount;
	A->rank = A->rowsCount;
	A->trace = 0.0;
	A->det = 0.0;
	A->permutationSign = 1;
	A->isSingular = false;
	A->isSPD = false;
	A->isRankDeficient = false;

	return A;
}

/**
 \fn	Mat freeMat (Mat A)

 \brief	Free memory consumed by Matrix.

 \date	15-May-14

 \param	A	The Matrix to process.
 */
void freeMat (Mat *A) {
	Assert$(*A != NULL, "Pointer cannot be NULL.");
	Assert$(((*A)->rowsCount != 0) && ((*A)->colsCount != 0), "Size must not be equal to 0.");
	Assert$((*A)->a != NULL, "Seems that there exists no contents to free...");

	for (size_t i = 0; i < (*A)->rowsCount; i++) {
		free((*A)->a[i]);
	}
	free$((*A)->a);
	free$(*A);

	return;
}

/**
 \fn	size_t freeMats (Mat A, ...)

 \brief	Deallocates many Matrices.

 \param	A	The Mat to process.

 \return	Deallocated matrices count.
 */
size_t freeMats (Mat A, ...) {
	Assert$(A != NULL, "Cannot free.");
	size_t n = 0;
	va_list vl;
	va_start(vl, A);
	for (Mat T = A; T; T = va_arg(vl, Mat), n++) {
		freeMat$(T);
	}
	va_end(vl);

	return n;
}
#pragma endregion "Alloc & dealloc"


#pragma region "Resizing"

/**
 \fn	void resize (Mat A, size_t newRows, size_t newCols)

 \brief	Resizes matrix.

 \param	A	   	The Mat to process.
 \param	newRows	The new rows count.
 \param	newCols	The new cols count.
 */
void resize (Mat A, size_t newRows, size_t newCols) {
	//C6385	Read overrun	Reading invalid data from 'A->a':  the readable size is 'newRows*sizeof(double *)' bytes, but '8' bytes may be read.	5	matrix.c	132
	//	'a' is an Input to 'realloc' (declared at c : \program files(x86)\microsoft visual studio 12.0\vc\include\stdlib.h:644)			127
	//	Enter this branch, (assume 'A!=(((void *)0))')			124
	//	Enter this branch, (assume 'newRows!=0&&newCols!=0')			125
	//	'a' is an Input to 'realloc' (declared at c : \program files(x86)\microsoft visual studio 12.0\vc\include\stdlib.h:644)			127
	//	Enter this branch, (assume 'newA!=(((void *)0))')			128
	//	Enter this loop, (assume '<branch condition>')			131
	//	Enter this branch, (assume 'newAi!=(((void *)0))')			133
	//	'i' may equal 1			131
	//	Continue this loop, (assume '<branch condition>')			131
	//	'i' is an Input to 'realloc' (declared at c : \program files(x86)\microsoft visual studio 12.0\vc\include\stdlib.h:644)			132
	//	Invalid read from 'A->a', (outside its readable range)			132


	Assert$(A != NULL, "Cannot resize");
	Assert$((newRows != 0 && newCols != 0), "Rows and/or columns count can't be set to 0");
	if ((A->rowsCount == newRows) && (A->colsCount == newCols)) {
		Check$(0, "Resizing has no effect because size difference is 0");
		return;
	}

	entry_t **newA = (entry_t**) realloc(A->a, newRows*sizeof(entry_t*));
	Assert$(newA != NULL, "Reallocating space for row pointers failed");
	A->a = newA;

	for (size_t i = 0; i < ((newRows < A->rowsCount) ? newRows : A->rowsCount); i++) {
		entry_t *newAi = (entry_t*) realloc(A->a[i], newCols*sizeof(entry_t));
		Assert$(newAi != NULL, "Reallocating space for rows failed");
		A->a[i] = newAi;
#ifdef __STDC_IEC_559__
		if (newCols > A->colsCount)	{
			//NOTE: memsetting double array w/ 0 will be UB if double isn't IEEE754-compliant.
			memset(A->a[i] + A->colsCount, 0, sizeof(entry_t) * (newCols - A->colsCount));
		}
#endif // __STDC_IEC_559__
	}

	if (newRows > A->rowsCount)	{
		for (size_t i = A->rowsCount; i < newRows; i++) {
//			A->a[i] = (entry_t*) calloc(newCols, sizeof(entry_t));
			A->a[i] = (entry_t*) malloc(newCols * sizeof(entry_t));
			Assert$(A->a[i] != NULL, "Allocating space for rows failed.");
		}
	}

	A->rowsCount = newRows;
	A->colsCount = newCols;

	return;
}

/**
 \fn	void concat (Mat A, Mat B)

 \brief	Merges matrices A & B into one, result will be in A.

 \param	A	The source Mat A to process.
 \param	B	The source Mat B to process.
 */
void concat (Mat A, Mat B) {
	Assert$((A != NULL) && (B != NULL), "Cannot concatenate matrices.");

	resize(A, max(A->rowsCount, B->rowsCount), A->colsCount + B->colsCount);

	entry_t **a = A->a;
	entry_t **b = B->a;
	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = A->colsCount - B->colsCount; j < A->colsCount; j++) {
			a[i][j] = b[i][j - (A->colsCount - B->colsCount)];
		}
	}

	return;
}
#pragma endregion "Resizing"


#pragma region "Printing"
/**
 \fn	void printMatToFile (Mat A, FILE *file, char *format)

 \brief	Prints matrix to file.

 \date	17-May-14

 \param	A			  	The Mat to process.
 \param [out]	file  	If non-null, the file to write to.
 \param [in]	format	If non-null, the format string.
 */
void printMatrixToFile (Mat A, FILE *file, char *format) {
	Assert$(file != NULL, "File reading error");
	Assert$(A != NULL, "Cannot print. Pointer is NULL.");

	entry_t **a = A->a;
#ifdef PRETTYOUTPUT
	static char buf[PRINTBUFSZ];
#endif // PRETTYOUTPUT

	for (size_t i = 0; i < A->rowsCount; i++) {
		fprintf(file, "[");
		for (size_t j = 0; j < A->colsCount; j++) {
#ifdef PRETTYOUTPUT
			snprintf(buf, PRINTBUFSZ-1, format, a[i][j]);
            buf[PRINTBUFSZ-1] = '\0';
			_trimTrailingZeroes(buf);
			fprintf(file, "%s|", buf);
#else
			fprintf(file, format, a[i][j]);
#endif
		}
		fprintf(file, "\b]\n");
	}

	return;
}

/**
 \fn	size_t printMatricesToFile (Mat A, ...)

 \brief	Print many matrices to file.

 \param	A	The Mat to process.

 \return	Printed matrices count.
 */
size_t printMatricesToFile (Mat A, ...) {
	Assert$(A != NULL, "Cannot print because of A is NULL.");
	size_t n = 0;
	va_list vl;
	va_start(vl, A);
	for (Mat T = A; T; T = va_arg(vl, Mat), n++) {
        printf(" >>Mat["FMT_SIZET"]\n", n);
		printMat$(T);
	}
	va_end(vl);

	return n;
}

/**
 \fn	void toString (Mat A, FILE *file, char *format)

 \brief	Prints matrix A in form of Mathematica-compatible string.

 \param	A			  	The Mat to process.
 \param [out]	file  	If non-null, the file to write to.
 \param [in]	format	If non-null, the format string.
 */
void toString (Mat A, FILE *file, char *format) {
	Assert$(file != NULL, "File access error.");
	Assert$(A != NULL, "Cannot print.");

	entry_t **a = A->a;
//#ifdef PRETTYOUTPUT
//	char buf[PRINTBUFSZ];
//#endif // PRETTYOUTPUT

	fprintf(file, "{");
	for (size_t i = 0; i < A->rowsCount; i++) {
		fprintf (file, "\n{");
		for (size_t j = 0; j < A->colsCount; j++) {
//#ifdef PRETTYOUTPUT
//			snprintf(buf, PRINTBUFSZ-1, format, a[i][j]);
//			_trimTrailingZeroes(buf); //TODO: buffer overflow
//			fprintf(file, "%s,", buf);
//#else
			fprintf(file, format, a[i][j]);
			fprintf(file, ",");
//#endif // PRETTYOUTPUT
		}
		fprintf(file, "\b},");
	}
	fprintf(file, "\b \n}\n");

	return;
}

#ifdef PRETTYOUTPUT
/**
 \fn	void _trimTrailingZeroes (char *str)

 \brief	Trims trailing zeroes in string containing decimal floating-point number.

 \param [in,out] str	If non-null, the string to process.
 */
void _trimTrailingZeroes (char *str) {
	Assert$(str != NULL, "");

	char *point = strchr(str, '.');
	if (point == NULL) {
		return;
	}
	char *lastZero = strrchr(point, '\0') - 1;
	if (lastZero != NULL) {
		while (*lastZero != '.') {
			if (*lastZero == '0') {
				memcpy(lastZero--, " ", 1);
			} else {
				break;
			}
		}
		if (*lastZero == '.') {
			memcpy(lastZero, " ", 1);
		}
	}
	char *s = strstr(str, "-0 ");
	if (s != NULL) {
		memcpy(s, " ", 1);
	}

	return;
}
//#else
//#define _trimTrailingZeroes
#endif // PRETTYOUTPUT
#pragma endregion "Printing"


#pragma region "Constructors"

/**
\fn	Mat DeepCopy (Mat A)

\brief	Constructs a deep copy of matrix A.

\date	22-May-14

\param	A	The source Mat to process.

\return	Deep copy of source Matrix.
*/
Mat DeepCopy (Mat A) {
	Assert$(A != NULL, "Cannot copy NULL...");
	Mat B = AllocMat(A->rowsCount, A->colsCount);
	Assert$(B != NULL, "Allocating memory for copy failed.");
	entry_t **a = A->a;
	entry_t **b = B->a;

	for (size_t i = 0; i < B->rowsCount; i++) {
		for (size_t j = 0; j < B->colsCount; j++) {
			b[i][j] = a[i][j];
		}
	}

	B->rank = A->rank;
	B->trace = A->trace;
	B->det = A->det;
	B->isSingular = A->isSingular;
	B->isSPD = A->isSPD;
	B->isRankDeficient = A->isRankDeficient;
	B->permutationSign = A->permutationSign;

	return B;
}

/**
\fn	Mat Copy (Mat A)

\brief	Constructs an shallow copy of matrix A.

\date	01-Jul-15

\param	A	The source Mat to process.

\return	Shallow copy of source Matrix.
*/
Mat Copy (Mat A) {
	Assert$(A != NULL, "Cannot copy NULL...");
	Mat B = AllocMat(A->rowsCount, A->colsCount);
	Assert$(B != NULL, "Allocating memory for copy failed.");
	entry_t **a = A->a;
	entry_t **b = B->a;

	for (size_t i = 0; i < B->rowsCount; i++) {
		for (size_t j = 0; j < B->colsCount; j++) {
			b[i][j] = a[i][j];
		}
	}

	return B;
}

/**
 \fn	Mat Diag (Mat A)

 \brief	Returns main diagonal of Matrix A in form of column-vector.

 \param	A	The Mat to process.

 \return	Main diagonal of Matrix A.
 */
Mat Diag (Mat A) {
	Mat D = AllocMat(1, A->rowsCount);
	entry_t **a = A->a;
	entry_t **d = D->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		d[0][i] = a[i][i];
	}

	return D;
}

/**
 \fn	Mat Sub (Mat A, size_t row, size_t col)

 \brief	Returns submatrix of A.

 \param	A  	The Mat to process.
 \param	row	The row.
 \param	col	The col.

 \return	A Mat.
 */
Mat SubMatrix (Mat A, size_t row, size_t col) {
	Assert$(false, "Not implemented"); //TODO:

	return NULL;
}

/**
 \fn	Mat Identity (size_t Size)

 \brief	Constructs an identity matrix with the given size.

 \date	30-May-14

 \param	Size	The size.

 \return	Identity Matrix.
 */
Mat Identity (size_t size) {
	Mat E = AllocMat(size, size);
	entry_t **e = E->a;

	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			if (i != j) {
				e[i][j] = 0.0;
			} else {
				e[i][j] = 1.0;
			}
		}
	}

	return E;
}

Mat Zeroes (size_t size) {
	Mat Z = AllocMat(size, size);
	fill_zeroes(Z);

	return Z;
}

/**
 \fn	Mat Minor (Mat A, size_t d)

 \brief	Returns leading Minor of A.

 \param	A	The Mat to process.
 \param	d	The order.

 \return	Minor.
 */
Mat Minor (Mat A, size_t d) { //TODO:
	Assert$(A->rowsCount <= d, "Rows count of A must be <= minor degree.");
	Mat M = AllocMat(A->rowsCount, A->colsCount);

	for (size_t i = 0; i < d; i++)	{
		M->a[i][i] = 1.0; //TODO: Array access results in a null pointer dereference
	}
	for (size_t i = d; i < A->rowsCount; i++) {
		for (size_t j = d; j < A->colsCount; j++) {
			M->a[i][j] = A->a[i][j];
		}
	}

	return M;
}
#pragma endregion "Constructors"


#pragma region "Filling routines"
/**
 \fn	size_t fill_fromFile (Mat A, FILE *file)

 \brief	Fills matrix from file.

 \date	17-May-14

 \param [out]	A   	The Matrix to fill to.
 \param [in]	file	If non-null, the file to read from.

 \return	Number of read items.
 */
size_t fill_fromFile (Mat A, FILE *file) {
	Assert$(A != NULL, "Cannot fill. It is NULL.");
	entry_t **a = A->a;
	entry_t tmp = 0.0;
	size_t r = 0;

	Assert$(A->rowsCount != 0 && A->colsCount != 0, "Invalid size.");
	Assert$(file != NULL, "File reading error");

	while (!(feof(file)) && !(ferror(file))) {
		spinActivityIndicator();
		for (size_t i = 0; i < A->rowsCount; i++) {
			for (size_t j = 0; j < A->colsCount; j++) {
				r += fscanf(file, FMT_FLT_INPUT, &tmp);
				a[i][j] = tmp;
			}
		}
	}
	clearActivityIndicator();

	return r;
}

/**
 \fn	void fill_random (Mat A)

 \brief	Fills matrix with random integer values.

 \date	30-May-14

 \param	A	The Mat to process.
 */
void fill_random (Mat A) {
	Assert$(A != NULL, "Cannot fill. It is NULL.");

	entry_t **a = A->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = round(((double) rand()) / 1000); //-V636
		}
	}

	return;
}

/**
 \fn	void fill_zeroes (Mat A)

 \brief	Fills double-valued matrix with zeroes.

 \date	04-Jun-14

 \param	A	The Mat to process.
 */
void fill_zeroes (Mat A) {
	Assert$(A != NULL, "Cannot fill. It is NULL.");

    entry_t **a = A->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = 0.0;
		}
	}

	return;
}

void fill_ones (Mat A) {
	Assert$(A != NULL, "Cannot fill. It is NULL.");

	entry_t **a = A->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = 1.0;
		}
	}

	return;
}

/**
 \fn	void fill_sequential (Mat A, int64_t start, int64_t inc)

 \brief	Fill Matrix A with sequential numbers starting from some value.

 \param	A	 	The Mat to process.
 \param	start	The start value.
 \param inc     The increment.
 */
void fill_sequential (Mat A, int64_t start, int64_t inc) {
	Assert$(A != NULL, "Cannot fill. It is NULL.");

    entry_t **a = A->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = (entry_t) start; //-V113
            start += inc;
		}
	}

	return;
}

/**
 \fn	void fill_spiral (Mat A, int64_t start)

 \brief	Spiral-like fill.

 \param	A	 	The Mat to process.
 \param	start	The start value.
 */
void fill_spiral (Mat A, int64_t start) {
	Assert$(A != NULL, "Cannot fill. It is NULL.");

	size_t sideLen = A->rowsCount;
	size_t numConcentricSquares = (size_t) (ceil(((double) (sideLen)) / 2.0));

	for (size_t i = 0; i < numConcentricSquares; i++) {
		// do top side
		for (size_t j = 0; j < sideLen; j++) {
			A->a[i][i + j] = (entry_t) start++;
		}

		// do right side
		for (size_t j = 1; j < sideLen; j++) {
			A->a[i + j][A->rowsCount - 1 - i] = (entry_t) start++;
		}

		// do bottom side
		for (ptrdiff_t j = sideLen - 2; j > -1; j--) {
			A->a[A->rowsCount - 1 - i][i + j] = (entry_t) start++;
		}

		// do left side
		for (ptrdiff_t j = sideLen - 2; j > 0; j--) {
			A->a[i + j][i] = (entry_t) start++;
		}

		sideLen -= 2;
	}
	return;
}

/**
 \fn	void fill_zigZag (Mat A, int64_t start)

 \brief	Fill 'zig-zag'-like.

 \param	A	 	The Mat to process.
 \param	start	The start value.
 */
void fill_zigZag (Mat A, int64_t start) {
	Assert$(A != NULL, "Cannot fill. It is NULL.");

	size_t lastValue = A->rowsCount * A->rowsCount - 1;
	size_t currDiag = 0;
	size_t loopFrom, loopTo;
	size_t row, col;
	do {
		if (currDiag < A->rowsCount) {// if doing the upper-left triangular half
			loopFrom = 0;
			loopTo = currDiag;
		} else {// doing the bottom-right triangular half
			loopFrom = currDiag - A->rowsCount + 1;
			loopTo = A->rowsCount - 1;
		}

		for (size_t i = loopFrom; i <= loopTo; i++)			{
			if (currDiag % 2 == 0) {// want to fill upwards
				row = loopTo - i + loopFrom;
				col = i;
			} else {// want to fill downwards
				row = i;
				col = loopTo - i + loopFrom;
			}

			A->a[row][col] = (entry_t) start++;
		}

		currDiag++;
	} while (currDiag <= lastValue);

	return;
}

void fill_tabulate (Mat A, Function2_u_u_e func) {
    Assert$(A != NULL, "Cannot fill. It is NULL.");

    for (size_t i = 0; i < A->rowsCount; ++i) {
        for (size_t j = 0; j < A->colsCount; ++j) {
            A->a[i][j] = (*func) (i, j);
        }
    }

    return;
}
#pragma endregion "Filling routines"
