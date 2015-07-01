#pragma once

/**
 \def	DOUBLE_PRECISION

 \brief	Enable double precision?
        Experimental. Do not touch!
 */
#define DOUBLE_PRECISION

/**
 \def	CHECKS_ENABLED

 \brief	Enable various checks?
 */
#define CHECKS_ENABLED

/**
 \def	ASSERTS_ENABLED

 \brief	Enable asserts?
 */
#define ASSERTS_ENABLED

/**
 \def	PRINTING_ENABLED();

 \brief	Will otutput be produced?
 */
#define PRINTING_ENABLED

/**
 \def	OUTFILE();

 \brief	Output file used in printing.
 */
#define OUTFILE stdout

/**
 \def	PRETTYOUTPUT();

 \brief	Will output be formatted to easy-readable form?
 */
//#define PRETTYOUTPUT

#ifdef PRETTYOUTPUT
/**
 \def	PRINTBUFSZ();

 \brief	Size of printing buffer.
        Note that buffer must be capable to hold string defined by FMT_FLT
 */
#define PRINTBUFSZ 1280
#endif // PRETTYOUTPUT

/**
 \def	FMT_INT();

 \brief	Formatting string for ints.
 */
#define FMT_INT	"%5d"

/**
 \def	FMT_FLT();

 \brief	Format string for floats.
 */
#define FMT_FLT	"%20.3f"

/**
 \def	FMT_FLT_STR();

 \brief	Format string for floats (used in printStr$()).
 */
#define FMT_FLT_STR		"% .e"

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
#define FMT_SIZET "%Iu"

/**
 \def	FMT_PTRDIFFT();

 \brief	Format string for ptrdiff_t on Windows. */
#define FMT_PTRDIFFT "%Id"
#else
#ifdef __GNUC__
/**
\def	FMT_SIZET();

\brief	Format string for size_t on Linux.
*/
#define FMT_SIZET "%zu"

/**
\def	FMT_PTRDIFFT();

\brief	Format string for ptrdiff_t on Linux. */
#define FMT_PTRDIFFT "%zd"

#ifdef __MINGW32__
#define __USE_MINGW_ANSI_STDIO 1 // TODO: doesn't work :C
#endif // __MINGW32__
#endif // __GNUC__
#endif // _MSC_VER //TODO:

#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER

#ifdef __MINGW32__
#include <stdio.h>
#define printf __mingw_printf
#endif // __MINGW32__

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif // _MSC_VER
