#pragma once

#include <math.h>
#include <float.h>
#include <stdint.h>

#include "Const.h"
#include "Config.h"

int64_t GCD_Euclid (int64_t a, int64_t b);

#define round(x)	    ( (int64_t)(((x)>=EPS)?(((x)+(0.5))):(((x)-(0.5)))) ) //Half away from zero 'rounding' with casting to int
#define lerp(a,b,t)	    ( (((entry_t)(1.0))-(t))*(a)+(t)*(b) )
#define swap_i(a,b)     ( ((a)!=(b))?(((a)^=(b)^=(a)^=(b))):(0) )
#define swap_f(a,b)	    do { float_t t=(a); (a)=(b); (b)=(t); } while (false)
#define swap_d(a,b)	    do { double_t t=(a); (a)=(b); (b)=(t); } while (false)
#define ispowerof2_i(x) ( ((x)&&(((x)&(~(x)+1))==(x))) )
#define square_i(x)     ( (x==0)?(0):((x)*(x)) )
#define square_fd(x)	( (x)*(x) )

//TODO: (more) proper floats comparison
#define equals_f(a,b)   ( ((fabsf((b)-(a)))<(F_EPS*fabsf((b)+(a)))) )
#define equals_d(a,b)   ( ((fabs((b)-(a)))<(D_EPS*fabs((b)+(a)))) )
//#define equals_d(a,b)   ( (a)==(b) )
//#define equals_d(a,b)   ( fabs((a)-(b))<=D_EPS )
//#define equals_d(a,b)   ( (fabs(a-b)<=D_EPS*max_3(1.0,fabs(a),fabs(b))) )

#define iszero_d(a)     ( fabs(a)<=D_EPS )
#define iszero_f(a)     ( fabsf(a)<=F_EPS )
#define isnotzero_d(a)  ( fabs(a)>D_EPS )
#define isnotzero_f(a)  ( fabsf(a)>F_EPS )

#ifdef DOUBLE_PRECISION
#define equals(a,b)     equals_d((a),(b))
#define iszero(a)       iszero_d((a))
#define isnotzero(a)    isnotzero_d((a))
#define swap(a,b)       swap_d((a),(b))
#define abs(a)          fabs((a))
#define hypot(a,b)      hypot((a),(b))
#else
#define equals(a,b)     equals_f((a),(b))
#define iszero(a)       iszero_f((a))
#define isnotzero(a)    isnotzero_f((a))
#define swap(a,b)       swap_f((a),(b))
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
