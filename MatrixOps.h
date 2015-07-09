#pragma once

#include <math.h>
#include <stdbool.h>

#include "Matrix.h"
#include "Config.h"


size_t Rank (Mat RREF);
entry_t Trace (Mat A);

bool IsIdentity (Mat A);
bool IsSingular (Mat A);
bool IsSymmetric (Mat A);
bool IsSkewSymmetric (Mat A);
bool IsEntriesEqual (Mat A, Mat B);
bool IsDimsEqual (Mat A, Mat B);

Mat Inverse (Mat A);
void toInverse (Mat* A);
void toTransposed_square (Mat A);
Mat Transposed (Mat A);
void toTransposed (Mat *A);

void matMul_inplace (Mat *A, Mat B);
Mat MatMul_naive (Mat A, Mat B);
Mat MatMul_naive_recursive (Mat A, Mat B);
Mat MatMul_Strassen (Mat A, Mat B);
Mat MatMul_Strassen_optimized (Mat A, Mat B);
size_t _fixSize (size_t Size);
Mat MatPow (Mat A, size_t deg);

#define MatMul$(A,B) MatMul_naive((A),(B))

Mat KroneckerProd (Mat A, Mat B);
Mat KroneckerSum (Mat A, Mat B);

entry_t OneNorm (Mat A);
entry_t TwoNorm (Mat A);
entry_t InfinityNorm (Mat A);
entry_t EuclideanNorm (Mat A);
entry_t ConditionNumber (Mat A);

entry_t DiagProd (Mat A);

Mat MatLerp_entrywise (Mat A, Mat B, entry_t t);

// 'templates' for element-wise functions
#define _foreach$(idx,mat,cmp)	for (size_t (idx) = 0; (idx) < (mat)->cmp; (idx)++)

#define TElementWise_MatrixMatrix1$(funcName, operator)				\
	void _mm1_##funcName (Mat A, Mat B) {							\
		_foreach$(i,A,rowsCount)									\
		_foreach$(j,A,colsCount) 									\
		A->a[i][j] = A->a[i][j] operator B->a[i][j]; }				\

extern void _mm1_add();
extern void _mm1_sub();
extern void _mm1_mul();
extern void _mm1_div();

#define TElementWise_MatrixMatrix2$(funcName, operator)				\
	void _mm2_##funcName (Mat A, Mat B, Mat C) {					\
		_foreach$(i,A,rowsCount)									\
		_foreach$(j,A,colsCount)									\
		(C->a[i][j]) = (A->a[i][j]) operator (B->a[i][j]); }		\

extern void _mm2_add();
extern void _mm2_sub();
extern void _mm2_mul();
extern void _mm2_div();

#define TElementWise_MatrixScalar$(funcName, operator)              \
	void _ms_##funcName (Mat A, double x) {					    	\
		_foreach$(i,A,rowsCount)									\
		_foreach$(j,A,colsCount)									\
		(A->a[i][j]) = (A->a[i][j]) operator (x); }					\

extern void _ms_add();
extern void _ms_sub();
extern void _ms_mul();
extern void _ms_div();
