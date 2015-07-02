#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>

#include "Config.h"
#include "Const.h"
#include "Datatypes.h"
#include "Functional.h"


Mat AllocMat (size_t SizeR, size_t SizeC);
void freeMat (Mat *A);
#define freeMat$(p) ( (p) ? (freeMat(&(p)), (p) = NULL) : (printf("Cannot free..."), NULL) )
size_t freeMats (Mat A, ...);

void resize (Mat A, size_t newRows, size_t newCols);
void concat (Mat A, Mat B);
#define IsSquare$(A) (((A)->rowsCount) == ((A)->colsCount))

void printMatrixToFile (Mat A, FILE *file, char *format);
size_t printMatricesToFile (Mat A, ...);
void toString (Mat A, FILE *file, char *format);

#ifdef PRETTYOUTPUT
extern void _cleanTrailingZeroes (char *str);
#else
#define _cleanTrailingZeroes(str)
#endif // PRETTYOUTPUT

#ifndef PRINTING_ENABLED
#define printMat$(a)
#define printMats$(a, ...)
#define printStr$(a)
#else
#define printMat$(a)       printMatrixToFile(a, OUTFILE, FMT_FLT)
#define printMats$(a, ...) printMatricesToFile(a, __VA_ARGS__)
#define printStr$(a)       toString(a, OUTFILE, FMT_FLT_STR)
#endif

size_t fill_fromFile (Mat A, FILE *file);
void fill_random (Mat A);
void fill_zeroes (Mat A);
void fill_sequential (Mat A, int64_t start, int64_t inc);
void fill_spiral (Mat A, int64_t start);
void fill_zigZag (Mat A, int64_t start);
void fill_tabulate (Mat A, entry_t (*func)(size_t, size_t));

Mat DeepCopy (Mat A);
Mat Copy (Mat A);
Mat Identity (size_t Size);
Mat Diag (Mat A);
Mat Minor (Mat A, size_t d);
