#pragma once

#include "Matrix.h"


Mat CholeskyDcmp (Mat A);
Mat Solve_cholesky (Mat L, Mat B);
double Det_cholesky (Mat L);
