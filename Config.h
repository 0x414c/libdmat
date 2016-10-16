#pragma once

#include "PlatformDependent.h"


/**
 \def	WITH_DOUBLE

 \brief	Enable double precision storage?
        If defined, Matrix entry type `entry_type` will be defined as `double_t`,
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
//#define WITH_DOUBLE


/**
 \def	RECURSION_LIMIT();

 \brief Matrix size limit for recursive procedures.
 */
#define RECURSION_LIMIT ( 8 )


/**
 \def	WITH_CHECKS();

 \brief	Enable checks?
 */
#define WITH_CHECKS


/**
 \def	WITH_ASSERTS();

 \brief	Enable asserts?
 */
#define WITH_ASSERTS


/**
 \def	WITH_PRINTING();

 \brief	Will output be produced?
 */
#define WITH_PRINTING


/**
 \def	OUT();

 \brief	Output file used in printing functions.
 */
#define OUT ( stdout )


/**
 \def	ERROUT();

 \brief	Output file used in printing functions.
 */
#define ERROUT ( stdout )


/**
 \def	WITH_PRETTYPRINT();

 \brief	Will output be formatted to easy-readable form?
 */
#define WITH_PRETTYPRINT


#ifdef WITH_PRETTYPRINT
/**
 \def	PRINTBUFSZ();

 \brief	Size of printing buffer.
        Note that buffer must be capable to hold string formatted using `FMT_FLT'
 */
#define PRINTBUFSZ ( 320 )
#endif //WITH_PRETTYPRINT


/**
 \def	FMT_INT();

 \brief	Formatting string for integers.
 */
#define FMT_INT "%5d"


/**
 \def	FMT_FLT();

 \brief	Format string for floats.
 */
#define FMT_FLT "%14.3f"


/**
 \def	FMT_FLT_STR();

 \brief	Format string for floats (used in `printAsStr$()').
 */
#define FMT_FLT_STR "% .10f"


/**
 \def	FMT_FLT_INPUT();

 \brief	Format string for filling from file containing floats (used in `fill_fromFile()').
 */
#ifdef WITH_DOUBLE
#define FMT_FLT_INPUT "%lf"
#else
#define FMT_FLT_INPUT "%f"
#endif //WITH_DOUBLE
