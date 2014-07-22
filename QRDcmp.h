#pragma once

#include "Matrix.h"

Mat *QRDcmp_householder (Mat A);
Mat Solve_qr (Mat *qr, Mat B);
Mat QRSolve_t (Mat A, Mat B);
double Det_qr (Mat *qr);
