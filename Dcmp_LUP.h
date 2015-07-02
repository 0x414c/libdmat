#pragma once

#include <stdbool.h>

#include "Matrix.h"


Mat *LUDcmp_gauss (Mat A);
Mat *LUDcmp_crout (Mat A);
Mat Pivotize_LU(Mat A);
double Det_LUP(Mat *LUP);
Mat Solve_LUP(Mat *lup, Mat B);

bool isSingular_LUP(Mat *LUP);
