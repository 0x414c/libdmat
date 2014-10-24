#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>


/**
 \def	DOCHECKS

 \brief	Enable various checks?
 */
#define DOCHECKS

/**
 \def	DOASSERTS

 \brief	Enable asserts?
 */
#define DOASSERTS


#ifdef DOCHECKS
#define Check(c,m)   ( (c)?(1):( fprintf(stderr, "Warning: %s %s - %s\n", __FUNCTION__, #c, m), 0) )
#else
#define Check(c,m)      
#endif // DOCHECKS

#ifdef DOASSERTS
#define Assert(c,m)	 ( (c)?(1):( fprintf(stderr, "Assertion (%s) failed ==> (%s)\n\tin %s() \n\tat <%s>, line %d\n", #c, m, __FUNCTION__, __FILE__, __LINE__), abort()) )
#else
#define Assert(c,m)	 
#endif // DOASSERTS


#define Free(p)      do { free(p); p = NULL; } while (0)

bool exists_d (long double *c, size_t s, long double x);
bool exists_u (size_t *c, size_t s, size_t e, size_t x);

bool nextCombination (size_t *index, ptrdiff_t k, ptrdiff_t n);

int64_t *iAllocVec (size_t Size);
size_t *uAllocVec (size_t Size);
long double *ldAllocVec (size_t Size);

void fillIndex (size_t *index, size_t size, size_t start);

int64_t GCD_euclid (int64_t a, int64_t b);

#define round(x)	 ( (int64_t)(((x)>=EPS)?(((x)+0.5)):(((x)-0.5))) ) //Half away from zero 'rounding' with casting to int
#define lerp(a,b,t)	 ( (1.0-(t))*(a)+(t)*(b) )
#define swap_d(a,b)	 { double t=(a); (a)=(b); (b)=(t); }
#define swap_i(a,b)	 ( ((a)!=(b))?(((a)^=(b)^=(a)^=(b))):(0) )
#define powerof2(x)  ( ((x)&&(((x)&(~(x)+1))==(x))) )
#define square_i(x)  ( (!x)?(0):((x)*(x)) )
#define square_d(x)	 ( (x)*(x) )
#define equal_d(a,b) ( ((fabs((b)-(a)))<(EPS*fabs((b)+(a)))) ) //TODO: (more) proper floats comparison
#define equal_ld(a,b) ( ((fabsl((b)-(a)))<(EPS*fabsl((b)+(a)))) ) //TODO: (more) proper floats comparison
//#define equal_d(a,b) ( fabs(a-b)<=EPS )
//#define equal_d(a,b) ( (fabs(a-b)<=EPS*max(1.0f,fabs(a),fabs(b))) )

#ifdef __GNUC__
#define min(a,b)	( ((a)<(b))?(a):(b) )
#define max(a,b)	( ((a)>(b))?(a):(b) )
#endif

#define NL puts("")
#define HR puts("------------------------------------------------------------------------------")
