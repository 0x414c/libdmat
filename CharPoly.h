#pragma once

#include <stdint.h>
#include <stddef.h>

#include "Matrix.h"


#define MAXITERS (1000)
#define DELTA (10)
#define START (-100000)
#define END (100000)


int64_t *GetCharPolyCoeffs (Mat A);
long double *GetPolynomialRoots (int64_t *c, size_t s, size_t *rootsCount);
long double FindRoot (int64_t *c, size_t s, long double a, long double b, size_t maxIterCount);
long double EvalPolyAt (long double x, int64_t *c, size_t s);

void printCharacteristicEquation (int64_t *c, size_t s, FILE *file);
