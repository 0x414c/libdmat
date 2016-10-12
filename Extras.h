#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "Config.h"


#ifdef WITH_CHECKS
#define Check$(c,m)		( (c) ? (true) : (fprintf(ERROUT, "Check `%s' have failed\n  with message `%s'\n  in function `%s'\n", #c, (m), __func__), false) )
#else //WITH_CHECKS
#define Check$(c, m)
#endif //WITH_CHECKS

#ifdef WITH_ASSERTS
#define Assert$(c,m)	( (c) ? (true) : (fprintf(ERROUT, "Assertion `%s' have failed\n  with message `%s'\n  in function `%s'\n  from file `%s'\n  at line %d\n", #c, (m), __func__, __FILE__, __LINE__), abort(), false) )
#else //WITH_ASSERTS
#define Assert$(c, m)
#endif //WITH_ASSERTS


#define free$(p)		do { free(p); (p) = (NULL); } while (false)

void _swap_v (void *x, void *y, void *tmp, size_t size);

#define swap$(a, b)		do { _swap_v(&(a), &(b), (char[(sizeof(a) == sizeof(b)) ? (ptrdiff_t) sizeof(a) : -1]) {0}, sizeof(a)); } while (false)


#define CRLF$			do { puts(""); } while (false)
#define HR$				do { puts("------------------------------------------------------------------------------"); } while (false)
