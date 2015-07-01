#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "Config.h"


#ifdef CHECKS_ENABLED
#define Check$(c,m)   ( (c)?(1):( fprintf(stderr, "Warning: check `%s` in `%s` failed ==> `%s`\n", #c, __FUNCTION__, m), 0) )
#else
#define Check$(c,m)
#endif // CHECKS_ENABLED

#ifdef ASSERTS_ENABLED
#define Assert$(c,m)  ( (c)?(1):( fprintf(stderr, "Warning: assertion `%s` failed ==> `%s`\n\tin `%s` \n\tat <%s>:%d\n", #c, m, __FUNCTION__, __FILE__, __LINE__), abort(), 0) )
#else
#define Assert$(c,m)
#endif // ASSERTS_ENABLED

#define free$(p)      do { free(p); p = NULL; } while (0)

bool exists_d (long double *c, size_t start, long double value);
bool exists_u (size_t *c, size_t start, size_t end, size_t value);

bool nextCombination (size_t *index, ptrdiff_t k, ptrdiff_t n);

int64_t *AllocVec_i (size_t Size);
size_t *AllocVec_u (size_t Size);
long double *AllocVec_ld (size_t Size);

void fillSequential_u (size_t *index, size_t size, size_t start);

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

#ifdef __GNUC__
#define min(a,b)	( ((a)<(b))?(a):(b) )
#define max(a,b)	( ((a)>(b))?(a):(b) )
#endif

#define CRLF$ puts("")
#define HR$ puts("------------------------------------------------------------------------------")
