#pragma once

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>

#include "Const.h"
#include "Config.h"


int64_t GCD_Euclid (int64_t a, int64_t b);

//Half away from zero rounding with casting to int64_t
#define round(x)		( (int64_t) (((x) >= (entry_type) 0.0) ? (((x) + (entry_type) 0.5)) : (((x) - (entry_type) 0.5))) )

#define lerp(a,b,t)		( ((entry_type) 1.0 - (t)) * (a) + (t) * (b) )

#define _swap_t(a,b,T)	do { T (t) = (a); (a) = (b); (b) = (t); } while (false)
#define _swap_i(a,b)	do { _swap_t((a), (b), int); } while (false)
#define _swap_f(a,b)	do { _swap_t((a), (b), float_t); } while (false)
#define _swap_d(a,b)	do { _swap_t((a), (b), double_t); } while (false)

#define ispowerof2_i(x)	( (x) && (((x) & (~(x) + 1)) == (x)) )

#define square_i(x)		( ((x) == 0) ? 0 : ((x) * (x)) )
#define square(x)		( (x) * (x) )

//TODO: (more) proper floats comparison
#define equals_d(a,b)	( fabs((a) - (b)) <= (EPS) * max_3((double_t) 1.0, fabs((a)), fabs((b))) )
#define equals_f(a,b)	( fabsf((a) - (b)) <= (EPS) * max_3((float_t) 1.0, fabsf((a)), fabsf((b))) )

#define iszero_d(a)		( fabs((a)) <= (EPS) )
#define iszero_f(a)		( fabsf((a)) <= (EPS) )

#define isnotzero_d(a)	( fabs((a)) > (EPS) )
#define isnotzero_f(a)	( fabsf((a)) > (EPS) )

#ifdef WITH_DOUBLE
#define equals(a,b)		( equals_d((a), (b)) )
#define iszero(a)		( iszero_d((a)) )
#define isnotzero(a)	( isnotzero_d((a)) )
#define swap(a,b)		do { _swap_d((a), (b)); } while (false)
#define abs(a)			( fabs((a)) )
#define hypot(a,b)		( hypot((a), (b)) )
#define sqrt(a)			( sqrt((a)) )
#define pow(a,b)		( pow((a), (b)) )
#else //WITH_DOUBLE
#define equals(a,b)		( equals_f((a), (b)) )
#define iszero(a)		( iszero_f((a)) )
#define isnotzero(a)	( isnotzero_f((a)) )
#define swap(a,b)		do { _swap_f((a), (b)); } while (false)
#define abs(a)			( fabsf((a)) )
#define hypot(a,b)		( hypotf((a), (b)) )
#define sqrt(a)			( sqrtf((a)) )
#define pow(a,b)		( powf((a), (b)) )
#endif //WITH_DOUBLE

#ifdef __GNUC__
#ifndef max
#define max(a,b)		( ((b) < (a)) ? (a) : (b) )
#endif //max
#ifndef min
#define min(a,b)		( ((b) < (a)) ? (b) : (a) )
#endif //min
#endif //__GNUC__

#define max_3(a,b,c)	( max(max((a), (b)), (c)) )
#define min_3(a,b,c)	( min(min((a), (b)), (c)) )
