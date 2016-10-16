#pragma once

#include "Matrix.h"


void toUpperTriangularForm (Mat A);
void toUpperTriangularForm_ref (Mat A);

void toRowEchelonForm (Mat A);

void toRowEchelonForm_Bareiss (Mat A);
void toReducedRowEchelonForm_Bareiss (Mat A);

void toReducedRowEchelonForm (Mat A);

entry_type Det_Gauss (Mat A);
entry_type Det_Bareiss (Mat A);

#define Det$(A) ( Det_Gauss((A)) )

Mat Solve_Gauss (Mat A, Mat B);
Mat Solve_GaussJordan (Mat A, Mat B);

Mat Solve_Bareiss (Mat A, Mat B);
Mat Solve_Montante (Mat A, Mat B);
