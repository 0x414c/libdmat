#pragma once

#include <math.h>
#include <stdbool.h>

#include "Matrix.h"


#define THRESHOLD 16


size_t Rank (Mat RREF);
double Trace (Mat A);

bool IsIdentity (Mat A);
//int IsSingular (Mat A);
bool IsSymmetric (Mat A);
bool IsEqual (Mat A, Mat B);

void toInverse (Mat A);
void toTransposed_square (Mat A);
Mat Transposed (Mat A);
void toTransposed (Mat *A);

Mat MatMul_naive (Mat A, Mat B);
Mat MatMul_naive_recursive (Mat A, Mat B);
Mat MatMul_strassen (Mat A, Mat B);
Mat MatMul_strassen_optimized (Mat A, Mat B);

size_t fixSize (size_t Size);
Mat MatPow (Mat A, size_t deg);
void matMul (Mat *A, Mat B);

#define MatMul(A,B) MatMul_naive(A,B)

Mat KroneckerProd (Mat A, Mat B);
Mat KroneckerSum (Mat A, Mat B);

double OneNorm (Mat A);
double TwoNorm (Mat A);
double InfinityNorm (Mat A);
double EuclideanNorm (Mat A);
double ConditionNumber (Mat A);

// 'templates' for element-wise funcNametions
#define for_i(A) for (size_t i = 0; i < A->rowsCount; i++)
#define for_j(A) for (size_t j = 0; j < A->colsCount; j++)

#define tElementWise2(funcName, operator)							        \
	void __##funcName (Mat A, Mat B) {							            \
	for_i(A) for_j(A) A->a[i][j] = A->a[i][j] operator B->a[i][j]; return; }\
void __add(); void __sub(); void __mul(); void __div();	

#define tElementWise3(funcName, operator)							        \
	void ___##funcName (Mat A, Mat B, Mat C) {							    \
	for_i(A) for_j(B) C->a[i][j] = A->a[i][j] operator B->a[i][j]; return; }\
void ___add(); void ___sub(); void ___mul(); void ___div();

#define tScalar(funcName, operator)                                  \
	void _s##funcName (Mat A, double x) {				       	     \
	for_i(A) for_j(A) A->a[i][j] = A->a[i][j] operator (x); return; }\
void _sadd(); void _ssub(); void _smul(); void _sdiv();	 
