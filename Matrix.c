#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define snprintf _snprintf
#endif

#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Matrix.h"
#include "Const.h"
#include "Extra.h"
#include "SpinningIndicator.h"


/**
 \fn	dMat dAllocMat (size_t RowsCount, size_t ColumnsCount)
 \brief	Allocate double-valued MxN matrix.
 \date	17-May-14						  
 \param	RowsCount		The rows count.
 \param	ColumnsCount	The columns count.
 \return	Pointer to _dMat struct.
 */
Mat AllocMat (size_t RowsCount, size_t ColumnsCount) {
	Mat A = NULL;

	Assert(((RowsCount != 0) && (ColumnsCount != 0)), "Matrix size can't be set to zero.");

	A = (Mat) malloc(sizeof(*A));
	Assert(A != NULL, "No space for matrix");
	
	A->a = NULL;
	A->a = (double**) calloc(RowsCount, sizeof(double*));
	Assert(A->a != NULL, "No space for row pointers");
	
	for (size_t i = 0; i < RowsCount; i++) {
		A->a[i] = NULL;
		A->a[i] = (double*) calloc(ColumnsCount, sizeof(double));
		Assert(A->a[i] != NULL, "No space for rows");
	}

	A->rowsCount = RowsCount;
	A->colsCount = ColumnsCount;
	A->rank = A->rowsCount;
	A->trace = 0.0;
	A->det = 0.0;
	A->isSingular = false;
	A->isSPD = false;
	A->isRankDeficient = false;
	A->permutationSign = 1;

	return A;
}

/**
 \fn	void FreeMat (dMat A)	
 \brief	Deallocate matrix.		
 \date		15-May-14			
 \param	A	The dMat to process.
 */
Mat freeMat (Mat A) {
	Assert(A != NULL, "Cannot free Mat."); //TODO: fucking 0xFEEEFEEE; upd: just do not free already freed mem chunks 
	Assert((A->rowsCount != 0) && (A->colsCount != 0), "Invalid size.");
	//Assert(A->a != NULL, "Seems that there exists no contents to free...");
	
	for (size_t i = 0; i < A->rowsCount; i++) {
		free(A->a[i]); A->a[i] = NULL;
	}
	free(A->a); A->a = NULL;
	free(A); A = NULL;

	return NULL; //TODO: ? 
}

size_t freeMats (Mat A, ...) {
	Assert(A != NULL, "Cannot free.");
	size_t n = 0;
	va_list vl;
	va_start(vl, A);
	for (Mat T = A; T; T = va_arg(vl, Mat), n++) {
		FreeMat(T);
	}
	va_end(vl);

	return n;
}

/**
 \fn	dMat DeepCopy (dMat A)
 \brief	Constructs a new deep copy of matrix A.	
 \date	22-May-14	 
 \param	A   	The dMat to process.		
 \return	A dMat.
 */
Mat DeepCopy (Mat A) {
	Assert(A != NULL, "Cannot copy NULL...");
	double **a = A->a;
	Mat B = AllocMat(A->rowsCount, A->colsCount);
	Assert(B != NULL, "Cannot create copy...");	//TODO: remove asserts like this?
	double **b = B->a;

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


	Assert(A != NULL, "Cannot resize.");
	Assert((newRows != 0 && newCols != 0), "Rows or columns count can't be set to 0.");
	if ((A->rowsCount == newRows) && (A->colsCount == newCols)) {
		Check(0, "Resizing has no effect.");
		return;
	}

	double **newA = (double**) realloc(A->a, newRows*sizeof(double*));
	Assert(newA != NULL, "Reallocating space for row pointers failed.");
	A->a = newA;

	for (size_t i = 0; i < ((newRows < A->rowsCount)? newRows: A->rowsCount); i++) {
		double *newAi = (double*) realloc(A->a[i], newCols*sizeof(double));
		Assert(newAi != NULL, "Reallocating space for rows failed.");
		A->a[i] = newAi;
		if (newCols > A->colsCount)	{
			memset(A->a[i] + A->colsCount, 0, sizeof(double) * (newCols - A->colsCount));
		}
	}

	if (newRows > A->rowsCount)	{
		for (size_t i = A->rowsCount; i < newRows; i++) {
			A->a[i] = (double*) calloc(newCols, sizeof(double));
			Assert(A->a[i] != NULL, "Allocating space for rows failed.");
		}
	}

	A->rowsCount = newRows;
	A->colsCount = newCols;

	return;
}

void concat (Mat A, Mat B) {
	Assert((A != NULL) && (B != NULL), "Cannot join matrices.");

	resize(A, max(A->rowsCount, B->rowsCount), A->colsCount + B->colsCount);

	double **a = A->a;
	double **b = B->a;
	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = A->colsCount - B->colsCount; j < A->colsCount; j++) {
			a[i][j] = b[i][j - (A->colsCount - B->colsCount)];
		}
	}

	return;
}

/**
 \fn	void PrintMatrixToFile (dMat A, FILE *file, char *format) 
 \brief	Prints matrix to file.								  
 \date	17-May-14											  
 \param	A							The dMat to process.
 \param [out]		file			If non-null, the file to write to.
 \param [in]		formatString	If non-null, the format string.
 */
void printMatrixToFile (Mat A, FILE *file, char *format) {
	Assert(file != NULL, "File reading error");
	Assert(A != NULL, "Cannot print.");

	double **a = A->a;
	char buf[PRINTBUFSZ];

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {			 
			snprintf(buf, PRINTBUFSZ-1, format, a[i][j]);
			_cleanTrailingZeroes(buf);
			fprintf(file, "%s", buf);
		}
		fprintf(file, "\n");
	}

	return;
}

size_t printMatricesToFile (Mat A, ...) {
	Assert(A != NULL, "Cannot print.");
	size_t n = 0;
	va_list vl;
	va_start(vl, A);
	for (Mat T = A; T; T = va_arg(vl, Mat), n++) {
		printf(" >>Mat[%Iu]\n", n);
		printMat(T);
	}
	va_end(vl);

	return n;
}


/**
\fn	void toString (dMat A, FILE *file, char *format)

\brief	Prints matrix A in form of Mathematica-compatible string.

\param	A						The dMat to process.
\param [out]	file			If non-null, the file to write to.
\param [in]		formatString	If non-null, the format string.
*/
void toString (Mat A, FILE *file, char *format) {
	Assert(file != NULL, "File reading error");
	Assert(A != NULL, "Cannot print.");

	double **a = A->a;
	char buf[PRINTBUFSZ];
	
	fprintf(file, "{");
	for (size_t i = 0; i < A->rowsCount; i++) {
		fprintf (file, "\n{");
		for (size_t j = 0; j < A->colsCount; j++) {
			snprintf(buf, PRINTBUFSZ-1, format, a[i][j]);
			//_cleanTrailingZeroes(buf);
			fprintf(file, "%s,", buf);
		}
		fprintf(file, "\b},");
	}
	fprintf(file, "\b \n}\n");

	return;
}

void _cleanTrailingZeroes (char *str) {
	char *s;
	char *start = strchr(str, '.');
	s = strrchr(start, '0');
	if ((s) && (*(s + 1) == '\0')) {
		while ((*s != '.')) {
			if (*s == '0') {
				memcpy(s--, " ", 1);
			} else {
				break;
			}
		}
	}
	s = strstr(str, "-0. ");
	if (s) {
		memcpy(s, " ", 1);
	}

	return;
}

Mat Diag (Mat A) {
	double **a = A->a;
	Mat D = AllocMat(1, A->colsCount);
	double **d = D->a;
	for (size_t i = 0; i < A->rowsCount; i++) {
		d[0][i] = a[i][i];
	}

	return D;
}

Mat Sub (Mat A, size_t row, size_t col) {
	Mat S = AllocMat(A->rowsCount / 2, A->colsCount / 2);
	return S; //TODO:
}

/**
 \fn	size_t fillMatrixFromFile (dMat A, FILE *file)	
 \brief	Fills matrix from file.	
 \date		17-May-14		
 \param	[out]		A		The dMat to fill to.
 \param [in]		file	If non-null, the file to read from.	 
 \return			Number of read items.
 */
size_t fillFromFile (Mat A, FILE *file) {
	Assert(A != NULL, "Cannot fill. It is NULL.");
	double **a = A->a;
	double tmp = 0.0;
	size_t r = 0;

	Assert(A->rowsCount != 0, "Invalid size.");
	Assert(file != NULL, "File reading error");

	while (!(feof(file))) {
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
 \fn	void randomFill (dMat A)   
 \brief	Fills matrix with random integer values.
 \date	30-May-14						 
 \param	A	The dMat to process.
 */
void fillRandom (Mat A) {
	Assert(A != NULL, "Cannot fill. It is NULL.");

	double **a = A->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = round((double) rand()/1000); //-V636
		}
	}

	return;
}

/**
 \fn	void zeroFill (dMat A) 
 \brief	Fills double-valued matrix with zeroes.	
 \date	04-Jun-14							  
 \param	A	The dMat to process.
 */						
void fillZero (Mat A) {
	Assert(A != NULL, "Cannot fill. It is NULL.");

	double **a = A->a;

	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = 0.0;
		}
	}

	return;
}

void fillNumbers (Mat A, int64_t start) {
	Assert(A != NULL, "Cannot fill. It is NULL.");

	double **a = A->a;
	
	for (size_t i = 0; i < A->rowsCount; i++) {
		for (size_t j = 0; j < A->colsCount; j++) {
			a[i][j] = (double) start++; //-V113
		}
	}

	return;
}

void fillSpiral (Mat A, int64_t start) {
	Assert(A != NULL, "Cannot fill. It is NULL.");

	size_t sideLen = A->rowsCount;
	size_t numConcentricSquares = (size_t) (ceil((double)(sideLen) / 2.0));

	for (size_t i = 0; i < numConcentricSquares; i++) {
		// do top side
		for (size_t j = 0; j < sideLen; j++) {
			A->a[i][i + j] = (double) start++;
		}

		// do right side
		for (size_t j = 1; j < sideLen; j++) {
			A->a[i + j][A->rowsCount - 1 - i] = (double) start++;
		}

		// do bottom side
		for (ptrdiff_t j = sideLen - 2; j > -1; j--) {
			A->a[A->rowsCount - 1 - i][i + j] = (double) start++;
		}

		// do left side
		for (ptrdiff_t j = sideLen - 2; j > 0; j--) {
			A->a[i + j][i] = (double) start++;
		}

		sideLen -= 2;
	}
	return;
}

void fillZigZag (Mat A, int64_t start) {
	Assert(A != NULL, "Cannot fill. It is NULL.");

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

			A->a[row][col] = (double) start++;
		}

		currDiag++;
	} while (currDiag <= lastValue);

	return;
}

/**
 \fn	dMat Identity (size_t Size)	
 \brief	Constructs an identity matrix with the given size. 
 \date	30-May-14		   
 \param	Size	The size.	 
 \return	A dMat.
 */
Mat Identity (size_t Size) {
	Mat I = AllocMat(Size, Size);

	for (size_t i = 0; i < Size; i++) {
		*(I->a[i] + i) = 1.0;
	}
	I->det = 1.0;

	return I;
}

Mat Minor (Mat A, size_t d) { //TODO: 
	Assert(A->rowsCount <= d, "");
	Mat M = AllocMat(A->rowsCount, A->colsCount);

	for (size_t i = 0; i < d; i++)	{
		M->a[i][i] = 1.0;
	}
	for (size_t i = d; i < A->rowsCount; i++) {
		for (size_t j = d; j < A->colsCount; j++) {
			M->a[i][j] = A->a[i][j];
		}
	}
	
	return M;
}
