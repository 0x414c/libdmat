#pragma once

#include "PlatformDependent.h"
#include "Const.h"


/**
 \def	DOUBLE_PRECISION

 \brief	Enable double precision storage?
        If defined, Matrix entry type `entry_t` will be defined as `double_t`,
        `float_t` otherwise.
        Computations will follow rules defined by `FLT_EVAL_METHOD`:
            -x  (not including `-1`) Implementation-defined.
            -1  Indeterminate.
             0  Evaluate all operations and constants just to the range and
                precision of the corresponding types.
                `float_t` and `double_t` are equivalent to `float` and `double` respectively.
             1  Evaluate operations and constants of type `float` and `double`
                to the range and precision of the `double` type,
                evaluate `long double` operations and constants to the range and
                precision of the `long double` type.
                Both `float_t` and `double_t` are equivalent to `double`.
             2  Evaluate all operations and constants to the range and
                precision of the `long double` type.
                Both `float_t` and `double_t` are equivalent to `long double`.
 */
//#define DOUBLE_PRECISION


/**
 \def	MM_SIZE_THRESHOLD

 \brief
 */
#define MM_SIZE_THRESHOLD (16)


/**
 \def	CHECKS_ENABLED

 \brief	Enable checks?
 */
#define CHECKS_ENABLED


/**
 \def	ASSERTS_ENABLED

 \brief	Enable asserts?
 */
#define ASSERTS_ENABLED


/**
 \def	PRINTING_ENABLED();

 \brief	Will output be produced?
 */
#define PRINTING_ENABLED


/**
 \def	OUTFILE();

 \brief	Output file used in printing functions.
 */
#define OUTFILE stdout


/**
 \def	PRETTYOUTPUT();

 \brief	Will output be formatted to easy-readable form?
 */
#define PRETTYOUTPUT


#ifdef PRETTYOUTPUT
/**
 \def	PRINTBUFSZ();

 \brief	Size of printing buffer.
        Note that buffer must be capable to hold string formatted w/ FMT_FLT
 */
#define PRINTBUFSZ (320)
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
#define FMT_FLT	"%14.3f"


/**
 \def	FMT_FLT_STR();

 \brief	Format string for floats (used in printAsStr$()).
 */
#define FMT_FLT_STR		"% .10f"


/**
 \def	FMT_FLT_INPUT();

 \brief	Format string for filling from file containing floats (used in fill_fromFile()).
 */
#ifdef DOUBLE_PRECISION
#define FMT_FLT_INPUT   "%lf"
#else
#define FMT_FLT_INPUT   "%f"
#endif // DOUBLE_PRECISION
