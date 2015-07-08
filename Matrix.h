#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>

#include "Config.h"
#include "Const.h"
#include "Datatypes.h"
#include "Functional.h"


Mat AllocMat (size_t rowsCount, size_t columnsCount);
void freeMat (Mat *A);
#define freeMat$(p) ( (p!=NULL)?(freeMat(&(p)),(p)=NULL,true):(printf("Cannot free memory."),false) )
size_t freeMats (Mat A, ...);

void resize (Mat A, size_t newRows, size_t newCols);
void concat (Mat A, Mat B);

Mat Copy (Mat A);
Mat DeepCopy (Mat A);
Mat Identity (size_t size);
Mat Zeroes (size_t size);
Mat Diag (Mat A);
Mat SubMatrix (Mat A, size_t row, size_t col);
Mat Minor (Mat A, size_t d);

size_t fill_fromFile (Mat A, FILE *file);
void fill_random (Mat A);
void fill_zeroes (Mat A);
void fill_ones (Mat A);
void fill_sequential (Mat A, int64_t start, int64_t inc);
void fill_tabulate (Mat A, entry_t (*func) (size_t, size_t));
void fill_spiral (Mat A, int64_t start);
void fill_zigZag (Mat A, int64_t start);

#define IsSquare$(A) (((A)->rowsCount) == ((A)->colsCount))

void printMatrixToFile (Mat A, FILE *file, char *format);
size_t printMatricesToFile (Mat A, ...);
void toString (Mat A, FILE *file, char *format);

#ifdef PRETTYOUTPUT
extern void _trimTrailingZeroes (char *str);
#else
#define _cleanTrailingZeroes(str)
#endif // PRETTYOUTPUT

#ifdef PRINTING_ENABLED
#define printMat$(a)       printMatrixToFile(a, OUTFILE, FMT_FLT)
#define printMats$(a, ...) printMatricesToFile(a, __VA_ARGS__)
#define printAsStr$(a)     toString(a, OUTFILE, FMT_FLT_STR)
#else
#define printMat$(a)
#define printMats$(a, ...)
#define printAsStr$(a)
#endif
