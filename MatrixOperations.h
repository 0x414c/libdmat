#pragma once

#include <math.h>
#include <stdbool.h>

#include "Matrix.h"


#define MM_SIZE_THRESHOLD 16


size_t Rank (Mat RREF);
entry_t Trace (Mat A);

bool IsIdentity (Mat A);
bool IsSingular (Mat A);
bool IsSymmetric (Mat A);
bool IsSkewSymmetric (Mat A);
bool IsEqual (Mat A, Mat B);

Mat Inverse (Mat A);
void toInverse (Mat* A);
void toTransposed_square (Mat A);
Mat Transposed (Mat A);
void toTransposed (Mat *A);

Mat MatMul_naive (Mat A, Mat B);
Mat MatMul_naive_recursive (Mat A, Mat B);
Mat MatMul_strassen (Mat A, Mat B);
Mat MatMul_strassen_optimized (Mat A, Mat B);

size_t _fixSize (size_t Size);
Mat MatPow (Mat A, size_t deg);
void matMul (Mat *A, Mat B);

#define MatMul$(A,B) MatMul_naive(A,B)

Mat KroneckerProd (Mat A, Mat B);
Mat KroneckerSum (Mat A, Mat B);

entry_t OneNorm (Mat A);
entry_t TwoNorm (Mat A);
entry_t InfinityNorm (Mat A);
entry_t EuclideanNorm (Mat A);
entry_t ConditionNumber (Mat A);

// 'templates' for element-wise functions
#define for_i$(A) for (size_t i = 0; i < A->rowsCount; i++)
#define for_j$(A) for (size_t j = 0; j < A->colsCount; j++)

#define TElementWise_MatrixMatrix1$(funcName, operator)						  \
	void _mm1_##funcName (Mat A, Mat B) {							          \
	for_i$(A) for_j$(A) A->a[i][j] = A->a[i][j] operator B->a[i][j]; return; }\
void _mm1_add(); void _mm1_sub(); void _mm1_mul(); void _mm1_div();

#define TElementWise_MatrixMatrix2$(funcName, operator)						  \
	void _mm2_##funcName (Mat A, Mat B, Mat C) {							  \
	for_i$(A) for_j$(B) C->a[i][j] = A->a[i][j] operator B->a[i][j]; return; }\
void _mm2_add(); void _mm2_sub(); void _mm2_mul(); void _mm2_div();

#define TElementWise_MatrixScalar$(funcName, operator)                 \
	void _ms_##funcName (Mat A, double x) {				       		   \
	for_i$(A) for_j$(A) A->a[i][j] = A->a[i][j] operator (x); return; }\
void _ms_add(); void _ms_sub(); void _ms_mul(); void _ms_div();
