#pragma once

#include "Matrix.h"

Mat *Dcmp_QR_Householder (Mat A);
Mat Solve_QR (Mat *QR, Mat B);
entry_type Det_QR (Mat *QR);
