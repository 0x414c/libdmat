#pragma once

#include <stddef.h>
#include <stdbool.h>
		  

typedef double entry_t; //TODO:

typedef struct _Mat_struct {
	entry_t **a;	
	size_t rowsCount;
	size_t colsCount;
	size_t rank;	
	bool isSingular;
	bool isSPD;
	bool isRankDeficient;
	int permutationSign;
	entry_t trace;
	entry_t det;
} *Mat;
