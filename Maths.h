#pragma once

#include <math.h>
#include <float.h>
#include <stdint.h>

#include "Const.h"
#include "Config.h"


int64_t GCD_Euclid (int64_t a, int64_t b);

#define round(x)	    ( (int64_t)(((x)>=EPS)?(((x)+(0.5))):(((x)-(0.5)))) ) //Half away from zero 'rounding' with casting to int
#define lerp(a,b,t)	    ( (((entry_t)(1.0))-(t))*(a)+(t)*(b) )
#define __swap__(a,b,T) do { T (t)=(a); (a)=(b); (b)=(t); } while (false)
#define __swap_i(a,b)   do { __swap__((a),(b),int); } while (false)
#define __swap_f(a,b)	do { __swap__((a),(b),float_t); } while (false)
#define __swap_d(a,b)	do { __swap__((a),(b),double_t); } while (false)
#define ispowerof2_i(x) ( ((x)&&(((x)&(~(x)+1))==(x))) )
#define square_i(x)     ( ((x)==0)?(0):((x)*(x)) )
#define square_fd(x)	( (x)*(x) )

//TODO: (more) proper floats comparison
#define equals_f(a,b)   ( (fabs(a-b)<=EPS*max_3(((float_t)(1.0)),fabs(a),fabs(b))) )
#define equals_d(a,b)   ( (fabs(a-b)<=EPS*max_3(((double_t)(1.0)),fabs(a),fabs(b))) )

#define iszero_d(a)     ( fabs(a)<=EPS )
#define iszero_f(a)     ( fabsf(a)<=EPS )
#define isnotzero_d(a)  ( fabs(a)>EPS )
#define isnotzero_f(a)  ( fabsf(a)>EPS )

#ifdef DOUBLE_PRECISION
#define equals(a,b)     equals_d((a),(b))
#define iszero(a)       iszero_d((a))
#define isnotzero(a)    isnotzero_d((a))
#define swap(a,b)       __swap_d((a),(b))
#define abs(a)          fabs((a))
#define hypot(a,b)      hypot((a),(b))
#else
#define equals(a,b)     equals_f((a),(b))
#define iszero(a)       iszero_f((a))
#define isnotzero(a)    isnotzero_f((a))
#define swap(a,b)       __swap_f((a),(b))
#define abs(a)          fabsf((a))
#define hypot(a,b)      hypotf((a),(b))
#endif // DOUBLE_PRECISION

//#ifdef __GNUC__
#ifndef max
#define max(a,b)	    ( ((a)<(b))?(b):(a) )
#endif // max
#ifndef min
#define min(a,b)	    ( ((a)<(b))?(a):(b) )
#endif // min
//#endif // __GNUC__

#define max_3(a,b,c)    ( max(max((a),(b)),(c)) )
#define min_3(a,b,c)    ( min(min((a),(b)),(c)) )
