#pragma once

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "Config.h"


#ifdef CHECKS_ENABLED
#define Check$(c,m)		( (c) ? (true) : (fprintf(stderr, "Check `%s' have failed:\n  with message `%s'\n  in function `%s'\n", #c, (m), __FUNCTION__), false) )
#else
#define Check$(c, m)
#endif //CHECKS_ENABLED

#ifdef ASSERTS_ENABLED
#define Assert$(c,m)	( (c) ? (true) : (fprintf(stderr, "Assertion `%s' have failed:\n  with message `%s'\n  in function `%s'\n  from file `%s'\n  at line %d\n", #c, (m), __FUNCTION__, __FILE__, __LINE__), abort(), false) )
#else
#define Assert$(c, m)
#endif //ASSERTS_ENABLED


#define free$(p)		do { free(p); (p) = (NULL); } while (false)

void _swap_v (void *x, void *y, void *tmp, size_t size);

#define swap$(a, b)		do { _swap_v(&(a), &(b), (char[(sizeof(a) == sizeof(b)) ? (ptrdiff_t) sizeof(a) : -1]) {0}, sizeof(a)); } while (false)


#define CRLF$			do { puts(""); } while (false)
#define HR$				do { puts("------------------------------------------------------------------------------"); } while (false)
