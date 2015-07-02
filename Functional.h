#pragma once

#include <stddef.h>

#include "Matrix.h"


typedef entry_t (*Function2_u_u_e) (size_t, size_t);

entry_t _add (size_t x, size_t y);
entry_t _mul (size_t x, size_t y);
entry_t _sub (size_t x, size_t y);
entry_t _div (size_t x, size_t y);
