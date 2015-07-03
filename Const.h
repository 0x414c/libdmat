#pragma once

#include <float.h>

#include "Config.h"


#ifdef DOUBLE_PRECISION
#define EPS (DBL_EPSILON)
#else
#define EPS (FLT_EPSILON)
#endif // DOUBLE_PRECISION
