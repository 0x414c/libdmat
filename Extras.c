#include <stddef.h>
#include <string.h>

#include "Extras.h"


__inline__ void _swap_v (void *x, void *y, void *tmp, size_t size) {
    memcpy(tmp, x, size);
    memcpy(x, y, size);
    memcpy(y, tmp, size);
}
