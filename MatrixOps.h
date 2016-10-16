#pragma once

#include <math.h>
#include <stdbool.h>

#include "Matrix.h"
#include "Config.h"


size_t Rank (Mat RREF);
entry_type Trace (Mat A);

bool IsIdentity (Mat A);
bool IsSingular (Mat A);
bool IsSymmetric (Mat A);
bool IsSkewSymmetric (Mat A);
bool IsEntriesEqual (Mat A, Mat B);
bool IsDimsEqual (Mat A, Mat B);

Mat Inverse_GaussJordan (Mat A);
//Mat Inverse_Strassen (Mat A);
void toInverse (Mat* A);

#define Inverse$(A) ( Inverse_GaussJordan((A)) )

Mat Transpose (Mat A);
void toTransposed (Mat *A);
void toTransposed_square (Mat A);

#define Transpose$(A) ( Transpose((A)) )

Mat MatMul_naive (Mat A, Mat B);
Mat MatMul_naive_recursive (Mat A, Mat B);
void matMul_inplace (Mat *A, Mat B);

Mat MatMul_Strassen (Mat A, Mat B);
//Mat MatMul_Strassen_optimized (Mat A, Mat B);
size_t _adjustSize (size_t size);

#define MatMul$(A,B) ( MatMul_naive((A), (B)) )

Mat MatPow_naive (Mat A, size_t pow);

#define MatPow$(A,e) ( MatPow_naive((A), (e)) )

Mat KroneckerProd (Mat A, Mat B);
Mat KroneckerSum (Mat A, Mat B);

entry_type OneNorm (Mat A);
entry_type TwoNorm (Mat A);
entry_type InfinityNorm (Mat A);
entry_type EuclideanNorm (Mat A);
entry_type ConditionNumber (Mat A);

entry_type DiagProd (Mat A);

Mat MatLerp_entrywise (Mat A, Mat B, entry_type t);

#define _foreach$(idx,data,cmp) for (size_t (idx) = 0; (idx) < (data)->cmp; (idx)++)

//'templates' for entrywise functions
#define EntryWiseFunc_MatrixMatrix_Inplace$(funcName, operator)		\
	void _##funcName##_m_m_inplace (Mat A, Mat B) {					\
		_foreach$(i, A, rowsCount)									\
		_foreach$(j, A, colsCount) 									\
		A->data[i][j] = A->data[i][j] operator B->data[i][j];		\
	}																\

extern void _add_m_m_inplace();
extern void _sub_m_m_inplace();
extern void _mul_m_m_inplace();
extern void _div_m_m_inplace();

#define EntryWiseFunc_MatrixMatrix$(funcName, operator)				\
	void _##funcName##_m_m_m (Mat A, Mat B, Mat C) {				\
		_foreach$(i, A, rowsCount)									\
		_foreach$(j, A, colsCount)									\
		C->data[i][j] = A->data[i][j] operator B->data[i][j];		\
	}																\

extern void _add_m_m_m();
extern void _sub_m_m_m();
extern void _mul_m_m_m();
extern void _div_m_m_m();

#define EntryWiseFunc_MatrixScalar_Inplace$(funcName, operator)		\
	void _##funcName##_m_s_inplace (Mat A, long double x) {			\
		_foreach$(i, A, rowsCount)									\
		_foreach$(j, A, colsCount)									\
		A->data[i][j] = A->data[i][j] operator x;					\
	}																\

extern void _add_m_s_inplace();
extern void _sub_m_s_inplace();
extern void _mul_m_s_inplace();
extern void _div_m_s_inplace();
