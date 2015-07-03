#pragma once

#include <math.h>
#include <stdint.h>

#include "Const.h"
#include "Config.h"

int64_t GCD_euclid (int64_t a, int64_t b);

#define round(x)	    ( (int64_t)(((x)>=EPS)?(((x)+(0.5))):(((x)-(0.5)))) ) //Half away from zero 'rounding' with casting to int
#define lerp(a,b,t)	    ( ((1.0)-(t))*(a)+(t)*(b) )
#define swap_d(a,b)	    { double t=(a); (a)=(b); (b)=(t); }
#define swap_i(a,b)     ( ((a)!=(b))?(((a)^=(b)^=(a)^=(b))):(0) )
#define ispowerof2_i(x) ( ((x)&&(((x)&(~(x)+1))==(x))) )
#define square_i(x)     ( (!x)?(0):((x)*(x)) )
#define square_d(x)	    ( (x)*(x) )
#define equals_d(a,b)   ( ((fabs((b)-(a)))<(EPS*fabs((b)+(a)))) ) //TODO: (more) proper floats comparison
#define equals_ld(a,b)  ( ((fabsl((b)-(a)))<(EPS*fabsl((b)+(a)))) ) //TODO: (more) proper floats comparison
//#define equals_d(a,b) ( fabs(a-b)<=EPS )
//#define equals_d(a,b) ( (fabs(a-b)<=EPS*max(1.0f,fabs(a),fabs(b))) )
#define iszero(a)       ( fabs(a)<=EPS )
#define nonzero(a)      ( fabs(a)>EPS )

#ifdef __GNUC__
#define min(a,b)	( ((a)<(b))?(a):(b) )
#define max(a,b)	( ((a)>(b))?(a):(b) )
#endif
