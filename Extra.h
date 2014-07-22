#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>


#define CHECK

#ifdef CHECK
#define Check(c,m)   ( (c)?(1):( fprintf(stderr, "Warning: %s - %s\n", #c, m), 0) )
#else
#define Check(c,m)   (c)   
#endif

#define Assert(c,m)	 ( (c)?(1):( fprintf(stderr, "Assertion (%s) failed ==> (%s)\n\tin %s() \n\tat <%s>, line %d\n", #c, m, __FUNCTION__, __FILE__, __LINE__), abort()) )

#define Free(p)      { free(p); p = NULL; }

bool exists_d (long double *c, size_t s, long double x);
bool exists_u (size_t *c, size_t s, size_t e, size_t x);

bool nextCombination (size_t *index, ptrdiff_t k, ptrdiff_t n);

int64_t *iAllocVec (size_t Size);
size_t *uAllocVec (size_t Size);
long double *ldAllocVec (size_t Size);

void fillIndex (size_t *index, size_t Size);

#define round(x)	 ( (int32_t)(((x)>=EPS)?(((x)+0.5)):(((x)-0.5))) ) //Half away from zero 'rounding'
#define swap_d(a,b)	 { double t=a; a=b; b=t; }
#define swap_i(a,b)	 ( ((a)!=(b))?(((a)^=(b)),((b)^=(a)),((a)^=(b))):(0) )
#define powerof2(x)  ( ((x)&&((x&(~x+1))==x)) )
#define square_i(x)  ( (!x)?(0):((x)*(x)) )
#define square_d(x)	 ( (x)*(x) )
#define equal_d(a,b) ( ((fabs(b-a))<(EPS*fabs(b+a)))) //TODO: (more) proper floats comparison
//#define eq_d(a,b)   ( a==b )

#ifdef __GNUC__
#define min(a,b)	( ((a)<(b))?(a):(b) )
#define max(a,b)	( ((a)>(b))?(a):(b) )
#endif

#define NL puts("")
#define HR puts("------------------------------------------------------------------------------")
