#pragma once

#include "Matrix.h"

Mat *Dcmp_QR_Householder (Mat A);
Mat Solve_QR (Mat *QR, Mat B);
entry_t Det_QR (Mat *QR);
