#pragma once

#include <stdint.h>
#include <stddef.h>

#include "Matrix.h"


Mat* EigenDcmp_n (Mat A);
Mat GetEigenvectors (Mat A, long double eigenvalue);
void fillEigenvectorMatrix (Mat EV, Mat OUT, size_t s, size_t e);
void fillEigenvalueMatrix (Mat L, long double *eVal, size_t c, size_t k);

Mat EDMatPow (Mat* res, size_t n);
Mat Spectrum (Mat* evd);
double SpectralRadius (Mat Sp);

int64_t GetPrincipalMinorsSum (Mat A, size_t order);
void strikeOut (Mat A, size_t d);
