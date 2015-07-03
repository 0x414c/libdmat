#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "Config.h"


#ifdef CHECKS_ENABLED
#define Check$(c,m)   ( (c)?(true):( fprintf(stderr, "Warning: check `%s` in `%s` failed ==> `%s`\n", #c, __FUNCTION__, m), false) )
#else
#define Check$(c,m)
#endif // CHECKS_ENABLED

#ifdef ASSERTS_ENABLED
#define Assert$(c,m)  ( (c)?(true):( fprintf(stderr, "Warning: assertion `%s` failed ==> `%s`\n\tin `%s` \n\tin file %s\n\tat line %d\n", #c, m, __FUNCTION__, __FILE__, __LINE__), abort(), false) )
#else
#define Assert$(c,m)
#endif // ASSERTS_ENABLED

#define free$(p)      do { free(p); p = NULL; } while (0)

bool exists_d (long double *c, size_t start, long double value);
bool exists_u (size_t *c, size_t start, size_t end, size_t value);

bool nextPermutation (size_t *index, ptrdiff_t k, ptrdiff_t n);

int64_t *AllocVec_i (size_t Size);
size_t *AllocVec_u (size_t Size);
long double *AllocVec_ld (size_t Size);

void fillSequential_u (size_t *index, size_t size, size_t start);

#define CRLF$ puts("")
#define HR$ puts("------------------------------------------------------------------------------")
