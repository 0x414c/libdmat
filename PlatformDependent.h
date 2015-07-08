#pragma once

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
#endif // __GNUC__
#endif // _MSC_VER //TODO:

#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER

//#ifdef __MINGW32__
//#define __USE_MINGW_ANSI_STDIO 1 // TODO: doesn't work :C
//#endif // __MINGW32__

#ifdef __MINGW32__
#include <stdio.h>
#define printf __mingw_printf
#define fprintf __mingw_fprintf
#endif // __MINGW32__

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif // _MSC_VER

#ifndef __GNUC__
#define __asm__ asm
#define __inline__ inline
#endif // __GNUC__
