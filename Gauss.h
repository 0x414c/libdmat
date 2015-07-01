#pragma once

#include "Matrix.h"


void toRowEchelonForm (Mat A);
void toRowEchelonForm_r (Mat A);
void toReducedRowEchelonForm (Mat A);

double Det_gauss (Mat A);
double Det_bareiss (Mat A);

Mat Solve_gaussjordan (Mat A, Mat B);
Mat Solve_gauss (Mat A, Mat B);

Mat GaussianSolve_h (Mat A);
void simpleSolver_h (double **a, size_t size, double *x);
void undeterminedSolver_h (Mat RREF, Mat A, Mat R);
