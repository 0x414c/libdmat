#pragma once

#include <stddef.h>

#include "Datatypes.h"


typedef entry_type (*function_s_s_e_t) (size_t, size_t);

entry_type _add_s_s_e (size_t x, size_t y);
entry_type _mul_s_s_e (size_t x, size_t y);
entry_type _sub_s_s_e (size_t x, size_t y);
entry_type _div_s_s_e (size_t x, size_t y);
