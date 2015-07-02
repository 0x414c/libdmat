#pragma once

#include "Matrix.h"


Mat Dcmp_Cholesky (Mat A);
Mat Solve_Cholesky (Mat L, Mat B);
entry_t Det_cholesky (Mat L);
