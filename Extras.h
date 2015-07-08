#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

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


#define free$(p)      do { free(p); p = NULL; } while (false)

void __swap (void *x, void *y, void* tmp, size_t size);

#define swap$(a, b)   do { __swap(&(a), &(b), (char[(sizeof(a)==sizeof(b))?(ptrdiff_t)sizeof(a):-1]) {0}, sizeof(a)); } while (false)


#define CRLF$ puts("")
#define HR$ puts("------------------------------------------------------------------------------")
