#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>

#include "Types.h"
#include "Const.h"


/**
 \def	DOPRINT();

 \brief	Will otutput be produced?
 */
#define DOPRINT

/**
 \def	OUTFILE();

 \brief	Output file used in printing.
 */
#define OUTFILE	   stdout

/**
 \def	PRETTYOUTPUT();

 \brief	Will output be in readable form?
 */
#define PRETTYOUTPUT

#ifdef PRETTYOUTPUT
/**
 \def	PRINTBUFSZ();

 \brief	Size of printing buffer.
 Note that buffer must be capable to hold string defined by FMT_FLT
 */
#define PRINTBUFSZ (64)
#endif // PRETTYOUTPUT


/**
 \def	FMT_INT();

 \brief	Formatting string for ints.
 */
#define FMT_INT	"%5d"

/**
 \def	FMT_FLT();

 \brief	Formatting string for floats.
 */
#define FMT_FLT	"%12.3f"

/**
 \def	FMT_FLT_STR();

 \brief	Format string for floats (used in printStr()).
 */
#define FMT_FLT_STR		"% .50g"

/**
 \def	FMT_FLT_INPUT();

 \brief	Format string for filling from file containing floats (fillFromFile()).
 */
#define FMT_FLT_INPUT   "%lf"


#ifdef _MSC_VER
/**
 \def	FMT_SIZET();

 \brief	Format string for size_t on Windows.
 */
#define FMT_SIZET "Iu"

/**
 \def	FMT_PTRDIFFT();

 \brief	Format string for ptrdiff_t on Windows. */
#define FMT_PTRDIFFT "Id"
#else
#ifdef __GNUC__
/**
\def	FMT_SIZET();

\brief	Format string for size_t on Linux.
*/
#define FMT_SIZET "zu"

/**
\def	FMT_PTRDIFFT();

\brief	Format string for ptrdiff_t on Linux. */
#define FMT_PTRDIFFT "zd"
#endif // __GNUC__
#endif // _MSC_VER //TODO:


Mat AllocMat (size_t SizeR, size_t SizeC);
void freeMat (Mat *A);
#define FreeMat(p) ( (p)?(freeMat(&p), p=NULL):(printf("Cannot free..."), NULL) )
size_t freeMats (Mat A, ...);

void resize (Mat A, size_t newRows, size_t newCols);
void concat (Mat A, Mat B);
#define square(A) (((A)->rowsCount)==((A)->colsCount))

void printMatrixToFile (Mat A, FILE *file, char *format);
size_t printMatricesToFile (Mat A, ...); 
void toString (Mat A, FILE *file, char *format);

#ifdef PRETTYOUTPUT
static void _cleanTrailingZeroes(char *str);
#else
#define _cleanTrailingZeroes()
#endif // PRETTYOUTPUT

#ifndef DOPRINT
#define printMat(a)
#define printMats(a, ...)
#define printStr(a)
#else
#define printMat(a)       printMatrixToFile(a, OUTFILE, FMT_FLT)
#define printMats(a, ...) printMatricesToFile(a, __VA_ARGS__)
#define printStr(a)       toString(a, OUTFILE, FMT_FLT_STR) 
#endif

size_t fillFromFile (Mat A, FILE *file);
void fillRandom (Mat A);
void fillZero (Mat A);
void fillNumbers (Mat A, int64_t start);
void fillSpiral (Mat A, int64_t start);
void fillZigZag (Mat A, int64_t start);

Mat DeepCopy (Mat A);
Mat Identity (size_t Size);
Mat Diag (Mat A);
Mat Minor (Mat A, size_t d);
