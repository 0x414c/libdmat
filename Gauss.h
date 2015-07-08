#pragma once

#include "Matrix.h"


void toRowEchelonForm (Mat A);
void toRowEchelonForm_reference (Mat A);
void toReducedRowEchelonForm (Mat A);

entry_t Det_Gauss (Mat A);
entry_t Det_Bareiss (Mat A);

Mat Solve_GaussJordan (Mat A, Mat B);
Mat Solve_Gauss (Mat A, Mat B);
