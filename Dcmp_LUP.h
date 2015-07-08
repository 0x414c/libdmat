#pragma once

#include <stdbool.h>

#include "Matrix.h"


Mat *Dcmp_LU_Gauss (Mat A);
Mat *Dcmp_LU_Crout (Mat A);
Mat Pivotize_LU (Mat A);
double Det_LUP (Mat *LUP);
Mat Solve_LUP (Mat *LUP, Mat B);

bool isSingular_LUP (Mat *LUP);
