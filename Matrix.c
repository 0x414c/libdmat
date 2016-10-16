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

 \return	Pointer to `_mat_s' if allocation was successful, NULL otherwise.
 */
Mat AllocMat (size_t rowsCount, size_t columnsCount) {
	Assert$((rowsCount > 0) && (columnsCount > 0), "Size should be greater than 0");

	Mat A = NULL;

	A = (Mat) malloc(sizeof(struct _mat_s));
	Assert$(A != NULL, "Cannot allocate memory for `_mat_s'.");

	A->data = NULL;
//	A->a = (entry_type**) calloc(rowsCount, sizeof(entry_type*));
	A->data = (entry_type**) malloc(rowsCount * sizeof(entry_type*));
	Assert$(A->data != NULL, "Cannot allocate memory for row pointers.");

	for (size_t i = 0; i < rowsCount; i++) {
		A->data[i] = NULL;
//		A->a[i] = (entry_type*) calloc(columnsCount, sizeof(entry_type));
		A->data[i] = (entry_type*) malloc(columnsCount * sizeof(entry_type));
		Assert$(A->data[i] != NULL, "Cannot allocate memory for row.");
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

 \brief	Free memory allocated for Matrix.

 \date	15-May-14

 \param	A	The Matrix to process.
 */
void freeMat (Mat *A) {
	Assert$(A != NULL, "*A should not be NULL.");
	Assert$((*A) != NULL, "A should not be NULL.");
	Assert$(((*A)->rowsCount > 0) && ((*A)->colsCount > 0), "Size should be greater than 0.");
	Assert$((*A)->data != NULL, "A contains no data.");

	for (size_t i = 0; i < (*A)->rowsCount; i++) {
		free((*A)->data[i]);
	}
	free$((*A)->data);
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
	Assert$(A != NULL, "A should not be NULL.");

	size_t n = 0;
	va_list vl;
	va_start(vl, A);

	for (Mat T = A; T != NULL; T = va_arg(vl, Mat), n++) {
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


	Assert$(A != NULL, "A should not be NULL.");
	Assert$((newRows > 0) && (newCols > 0), "Rows and/or columns count cannot be set to 0.");
	if ((A->rowsCount == newRows) && (A->colsCount == newCols)) {
		Check$(false, "Resizing has no effect because size difference is 0.");

		return;
	}

	entry_type **newA = (entry_type**) realloc(A->data, newRows*sizeof(entry_type*));
	Assert$(newA != NULL, "Reallocating memory for row pointers failed.");
	A->data = newA;

	for (size_t i = 0; i < ((newRows < A->rowsCount) ? newRows : A->rowsCount); i++) {
		entry_type *newAi = (entry_type*) realloc(A->data[i], newCols * sizeof(entry_type));
		Assert$(newAi != NULL, "Reallocating memory for rows failed.");
		A->data[i] = newAi;
#ifdef __STDC_IEC_559__
		if (newCols > A->colsCount)	{
			//NOTE: memsetting double array w/ 0 will be UB if double isn't IEEE754-compliant.
			memset(A->a[i] + A->colsCount, 0, sizeof(entry_type) * (newCols - A->colsCount));
		}
#endif //__STDC_IEC_559__
	}

	if (newRows > A->rowsCount)	{
		for (size_t i = A->rowsCount; i < newRows; i++) {
//			A->a[i] = (entry_type*) calloc(newCols, sizeof(entry_type));
			A->data[i] = (entry_type*) malloc(newCols * sizeof(entry_type));
			Assert$(A->data[i] != NULL, "Allocating memoty for rows failed.");
		}
	}

	A->rowsCount = newRows;
	A->colsCount = newCols;

	return;
}

/**
 \fn	void joinColumns (Mat A, Mat B)

 \brief	Merges matrices A & B into one, result will be in A.

 \param	A	The source Mat A to process.
 \param	B	The source Mat B to process.
 */
void joinColumns (Mat A, Mat B) {
	Assert$((A != NULL) && (B != NULL), "A and B should not be NULL.");

	resize(A, max(A->rowsCount, B->rowsCount), A->colsCount + B->colsCount);

	entry_type **a = A->data;
	entry_type **b = B->data;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = A->colsCount - B->colsCount; j < A->colsCount; j++) {
			a[i][j] = b[i][j - (A->colsCount - B->colsCount)];
		}
	}

	return;
}
#pragma endregion "Resizing"


#pragma region "Swap*"
void swapRows (Mat A, size_t i, size_t j) {
	Assert$(A != NULL, "A should not be NULL.");

	if (Check$(i != j, "Swap has no effect because source and target rows are equal.")) {
		entry_type **a = A->data;

		entry_type *a_i = a[i];
		a[i] = a[j];
		a[j] = a_i;

//		for (size_t k = 0; k < A->colsCount; k++) {
//			swap(a[i][k], a[j][k]); //TODO: Swap row pointers
//		}

//		A->permutationSign *= -1; //-V127 //TODO:
	}

	return;
}

void swapCols (Mat A, size_t i, size_t j) {
	Assert$(A != NULL, "A should not be NULL.");

	if (Check$(i != j, "Swap has no effect because source and target columns are equal.")) {
		entry_type **a = A->data;

		for (size_t k = 0; k < A->rowsCount; k++) {
			swap(a[k][i], a[k][j]);
		}
	}

	return;
}
#pragma endregion "Swap*"


#pragma region "Printing"
/**
 \fn	void printMatToFile (Mat A, FILE *file, char *format)

 \brief	Prints matrix to file.

 \date	17-May-14

 \param	A			  	The Mat to process.
 \param [out]	file  	If non-null, the file to write to.
 \param [in]	format	If non-null, the format string.
 */
void printMatToFile (Mat A, FILE *file, char *format) {
	Assert$(A != NULL, "A should not be NULL.");
	Assert$(file != NULL, "file should not be NULL.");
	Assert$(format != NULL, "format should not be NULL.");

	entry_type **a = A->data;
#ifdef WITH_PRETTYPRINT
	static char buf[PRINTBUFSZ];
#endif //WITH_PRETTYPRINT

	for (size_t i = 0; i < A->rowsCount; i++) {
		fprintf(file, "[");
		fflush(file);
		for (size_t j = 0; j < A->colsCount; j++) {
#ifdef WITH_PRETTYPRINT
			snprintf(buf, PRINTBUFSZ - 1, format, a[i][j]);
            buf[PRINTBUFSZ-1] = '\0';
			_trimTrailingZeroes(buf);
			fprintf(file, "%s|", buf);
			fflush(file);
#else //WITH_PRETTYPRINT
			fprintf(file, format, a[i][j]);
			fflush(file);
#endif //WITH_PRETTYPRINT
		}
		fprintf(file, "\b]\n");
		fflush(file);
	}

	return;
}

/**
 \fn	size_t printMatsToFile (Mat A, ...)

 \brief	Print many matrices to file.

 \param	A	The Mat to process.

 \return	Printed matrices count.
 */
size_t printMatsToFile (Mat A, ...) {
	Assert$(A != NULL, "A should not be NULL.");

	size_t n = 0;
	va_list vl;
	va_start(vl, A);
	for (Mat T = A; T != NULL; T = va_arg(vl, Mat), n++) {
		printf(" > Mat[" FMT_SIZET "]:\n", n);
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
	Assert$(A != NULL, "A should not be NULL.");
	Assert$(file != NULL, "file should not be NULL.");
	Assert$(format != NULL, "format should not be NULL.");

	entry_type **a = A->data;
//#ifdef WITH_PRETTYPRINT
//	char buf[PRINTBUFSZ];
//#endif //WITH_PRETTYPRINT

	fprintf(file, "{");
	for (size_t i = 0; i < A->rowsCount; i++) {
		fprintf (file, "\n  {");
		for (size_t j = 0; j < A->colsCount; j++) {
//#ifdef WITH_PRETTYPRINT
//			snprintf(buf, PRINTBUFSZ-1, format, a[i][j]);
//			_trimTrailingZeroes(buf); //TODO: buffer overflow
//			fprintf(file, "%s,", buf);
//#else //WITH_PRETTYPRINT
			fprintf(file, format, a[i][j]);
			fprintf(file, ",");
//#endif //WITH_PRETTYPRINT
		}
		fprintf(file, "\b},");
	}
	fprintf(file, "\b \n}\n");

	return;
}

#ifdef WITH_PRETTYPRINT
/**
 \fn	void _trimTrailingZeroes (char *str)

 \brief	Trims trailing zeroes in string containing decimal floating-point number.

 \param [in,out] str	If non-null, the string to process.
 */
void _trimTrailingZeroes (char *str) {
	Assert$(str != NULL, "str should not be NULL.");

	char *point = strchr(str, '.');
	if (point == NULL) {
		return;
	} else {
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
	}

	return;
}
#endif //WITH_PRETTYPRINT
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
	Assert$(A != NULL, "A should not be NULL.");

	entry_type **a = A->data;

	Mat A1 = AllocMat(A->rowsCount, A->colsCount);
	Assert$(A1 != NULL, "Allocating memory for copy failed.");

	entry_type **b = A1->data;

	for (size_t i = 0; i < A1->rowsCount; i++) {
		for (size_t j = 0; j < A1->colsCount; j++) {
			b[i][j] = a[i][j];
		}
	}

	A1->rank = A->rank;
	A1->trace = A->trace;
	A1->det = A->det;

	A1->isSingular = A->isSingular;
	A1->isIdentity = A->isIdentity;
	A1->isSPD = A->isSPD;
	A1->isRankDeficient = A->isRankDeficient;
	A1->permutationSign = A->permutationSign;

	return A1;
}

/**
\fn	Mat Copy (Mat A)

\brief	Constructs an shallow copy of matrix A.

\date	01-Jul-15

\param	A	The source Mat to process.

\return	Shallow copy of source Matrix.
*/
Mat Copy (Mat A) {
	Assert$(A != NULL, "A should not be NULL.");

	entry_type **a = A->data;

	Mat A1 = AllocMat(A->rowsCount, A->colsCount);
	Assert$(A1 != NULL, "Allocating memory for copy failed.");
	entry_type **b = A1->data;

	for (size_t i = 0; i < A1->rowsCount; i++) {
		for (size_t j = 0; j < A1->colsCount; j++) {
			b[i][j] = a[i][j];
		}
	}

	return A1;
}

/**
 \fn	Mat Diag (Mat A)

 \brief	Returns main diagonal of Matrix A in form of column-vector.

 \param	A	The Mat to process.

 \return	Main diagonal of Matrix A.
 */
Mat Diag (Mat A) {
	Assert$(A != NULL, "A should not be NULL.");

	entry_type **a = A->data;
	Mat D = AllocMat(1, A->rowsCount);
	entry_type **d = D->data;

	for (size_t i = 0; i < A->rowsCount; i++) {
		d[0][i] = a[i][i];
	}

	return D;
}

/**
 \fn	Mat SubMat (Mat A, size_t row, size_t col)

 \brief	Returns submatrix of A.

 \param	A  	The Mat to process.
 \param	row	The row.
 \param	col	The col.

 \return	A Mat.
 */
Mat SubMat (Mat A, size_t row, size_t col) {
	Assert$(false, "Not Yet Implemented."); //TODO:

	return NULL;
}

/**
 \fn	Mat Identity (size_t size)

 \brief	Constructs an identity matrix with the given size.

 \date	30-May-14

 \param	Size	The size.

 \return	Identity Matrix.
 */
Mat Identity (size_t size) {
	Assert$(size > 0, "size should be greater than 0.");

	Mat E = AllocMat(size, size);
	entry_type **e = E->data;

	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			if (i != j) {
				e[i][j] = 0.0;
			} else {
				e[i][j] = 1.0;
			}
		}
	}

	E->isIdentity = true;

	return E;
}

Mat Zeroes (size_t size) {
	Mat Z = AllocMat(size, size);
	fill_zeroes(Z);

	return Z;
}

Mat Ones (size_t size) {
	Mat O = AllocMat(size, size);
	fill_ones(O);

	return O;
}

/**
 \fn	Mat Minor (Mat A, size_t d)

 \brief	Returns leading Minor of A.

 \param	A	The Mat to process.
 \param	d	The order.

 \return	Minor.
 */
Mat Minor (Mat A, size_t d) { //TODO:
	Assert$(A != NULL, "A should not be NULL.");
	Assert$(A->rowsCount <= d, "Rows count of A should be less than minor degree.");

	Mat M = AllocMat(A->rowsCount, A->colsCount);

	for (size_t i = 0; i < d; i++)	{
		M->data[i][i] = 1.0; //TODO: Array access results in a null pointer dereference
	}
	for (size_t i = d; i < A->rowsCount; i++) {
		for (size_t j = d; j < A->colsCount; j++) {
			M->data[i][j] = A->data[i][j];
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
	Assert$(A != NULL, "A should not be NULL.");
	Assert$(A->rowsCount != 0 && A->colsCount != 0, "Size of matrix A should be greater than 0.");
	Assert$(file != NULL, "file should not be NULL.");

	entry_type **a = A->data;
	entry_type tmp = 0.0;
	size_t r = 0;

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
	Assert$(A != NULL, "A should not be NULL.");

	entry_type **a = A->data;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = round(((entry_type) rand()) / 1000); //-V636
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
	Assert$(A != NULL, "A should not be NULL.");

	entry_type **a = A->data;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = 0.0;
		}
	}

	return;
}

void fill_ones (Mat A) {
	Assert$(A != NULL, "A should not be NULL.");

	entry_type **a = A->data;

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
	Assert$(A != NULL, "A should not be NULL.");

	entry_type **a = A->data;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = (entry_type) start; //-V113
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
	Assert$(A != NULL, "A should not be NULL.");
	Assert$(IsSquare$(A), "Matrix A should be square."); //TODO:

	size_t sideLen = A->rowsCount;
	size_t numConcentricSquares = (size_t) (ceil(((double) (sideLen)) / 2.0));

	for (size_t i = 0; i < numConcentricSquares; i++) {
		//do top side
		for (size_t j = 0; j < sideLen; j++) {
			A->data[i][i + j] = (entry_type) start++;
		}

		//do right side
		for (size_t j = 1; j < sideLen; j++) {
			A->data[i + j][A->rowsCount - 1 - i] = (entry_type) start++;
		}

		//do bottom side
		for (ptrdiff_t j = sideLen - 2; j > -1; j--) {
			A->data[A->rowsCount - 1 - i][i + j] = (entry_type) start++;
		}

		//do left side
		for (ptrdiff_t j = sideLen - 2; j > 0; j--) {
			A->data[i + j][i] = (entry_type) start++;
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
	Assert$(A != NULL, "A should not be NULL.");
	Assert$(IsSquare$(A), "Matrix A should be square."); //TODO:

	size_t lastValue = A->rowsCount * A->rowsCount - 1;
	size_t currDiag = 0;
	size_t loopFrom, loopTo;
	size_t row, col;

	do {
		if (currDiag < A->rowsCount) { //if doing the upper-left triangular half
			loopFrom = 0;
			loopTo = currDiag;
		} else { //doing the bottom-right triangular half
			loopFrom = currDiag - A->rowsCount + 1;
			loopTo = A->rowsCount - 1;
		}

		for (size_t i = loopFrom; i <= loopTo; i++) {
			if (currDiag % 2 == 0) { //want to fill upwards
				row = loopTo - i + loopFrom;
				col = i;
			} else { //want to fill downwards
				row = i;
				col = loopTo - i + loopFrom;
			}

			A->data[row][col] = (entry_type) start++;
		}

		currDiag++;
	} while (currDiag <= lastValue);

	return;
}

void fill_tabulate (Mat A, function_s_s_e_t func) {
	Assert$(A != NULL, "A should not be NULL.");
	Assert$(func != NULL, "func should not be NULL.");

	for (size_t i = 0; i < A->rowsCount; ++i) {
		for (size_t j = 0; j < A->colsCount; ++j) {
			A->data[i][j] = (*func) (i, j);
		}
	}

	return;
}
#pragma endregion "Filling routines"
