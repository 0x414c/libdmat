#pragma once

#include <stdint.h>
#include <stddef.h>

#include "Extras.h"
#include "Maths.h"
#include "Const.h"


bool nextPermutation (size_t *index, ptrdiff_t k, ptrdiff_t n);

bool exists_d (double *c, size_t start, double value);
bool exists_u (size_t *c, size_t start, size_t end, size_t value);

int64_t *AllocVec_i (size_t Size);
size_t *AllocVec_u (size_t Size);
double *AllocVec_d (size_t Size);

void fillSequential_u (size_t *index, size_t size, size_t start);
