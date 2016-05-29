#pragma once

#include <stddef.h>

#include "Datatypes.h"


typedef entry_t (*function_s_s_e_t) (size_t, size_t);

entry_t _add_s_s_e (size_t x, size_t y);
entry_t _mul_s_s_e (size_t x, size_t y);
entry_t _sub_s_s_e (size_t x, size_t y);
entry_t _div_s_s_e (size_t x, size_t y);
