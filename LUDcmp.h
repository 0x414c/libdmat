#pragma once

#include <stdbool.h>

#include "Matrix.h"


Mat *LUDcmp_gauss (Mat A);
Mat *LUDcmp_crout (Mat A);
Mat LUPivotize (Mat A);
double Det_lup (Mat *lup);
Mat Solve_lup (Mat *lup, Mat B);

bool isSingular_lup (Mat *lup);
