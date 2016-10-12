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
#define freeMat$(p) ( ((p) != (NULL)) ? (freeMat(&(p)), (p) = (NULL), true) : (fprintf(stderr, "Cannot deallocate memory."), false) )
size_t freeMats (Mat A, ...);

void resize (Mat A, size_t newRows, size_t newCols);
void join (Mat A, Mat B);

void swapRows (Mat A, size_t i, size_t j);
void swapCols (Mat A, size_t i, size_t j);

Mat Copy (Mat A);
Mat DeepCopy (Mat A);
Mat Identity (size_t size);
Mat Zeroes (size_t size);
Mat Diag (Mat A);
Mat SubMat (Mat A, size_t row, size_t col);
Mat Minor (Mat A, size_t d);

size_t fill_fromFile (Mat A, FILE *file);
void fill_random (Mat A);
void fill_zeroes (Mat A);
void fill_ones (Mat A);
void fill_sequential (Mat A, int64_t start, int64_t inc);
void fill_tabulate (Mat A, function_s_s_e_t func);
void fill_spiral (Mat A, int64_t start);
void fill_zigZag (Mat A, int64_t start);

#define IsSquare$(A) ( ((A)->rowsCount) == ((A)->colsCount) )

void printMatToFile (Mat A, FILE *file, char *format);
size_t printMatsToFile (Mat A, ...);
void toString (Mat A, FILE *file, char *format);

#ifdef WITH_PRETTYPRINT
extern void _trimTrailingZeroes (char *str);
#else //WITH_PRETTYPRINT
#define _trimTrailingZeroes(str)
#endif //WITH_PRETTYPRINT

#ifdef WITH_PRINTING
#define printMat$(a)       do { printMatToFile((a), OUT, FMT_FLT); } while (false)
#define printMats$(a, ...) do { printMatsToFile((a), __VA_ARGS__); } while (false)
#define printAsStr$(a)     do { toString((a), OUT, FMT_FLT_STR); } while (false)
#else
#define printMat$(a)
#define printMats$(a, ...)
#define printAsStr$(a)
#endif //WITH_PRINTING
