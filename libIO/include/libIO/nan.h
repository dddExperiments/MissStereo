#ifndef NAN_H
#define NAN_H

#include <math.h>

/// Not a number, only for setting a variable.
static const float NaN=sqrt(-1.0f);

/// Test if \a x is \em not NaN
inline bool is_number(float x) {
    return (x==x);
}

#endif
