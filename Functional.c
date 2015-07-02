#include <stddef.h>

#include "Functional.h"


entry_t _add (size_t x, size_t y) {
    return x + y;
}

entry_t _mul (size_t x, size_t y) {
    return x * y;
}

entry_t _sub (size_t x, size_t y) {
    return x - y;
}

entry_t _div (size_t x, size_t y) {
    return x / y;
}
