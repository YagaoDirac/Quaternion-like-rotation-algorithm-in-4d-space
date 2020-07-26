#pragma once

#include <cmath>

#define __D4D_REAL_FLOAT__ 1

#ifdef __D4D_REAL_FLOAT__
typedef float real;
#endif // __D4D_REAL_FLOAT__

real D4D_default_epsilon();
real D4D_pi();



/*DO NOT replace this function with any const value. Or you need to prepare a cpp file for it.
In that way, files should look like:

//.h
extern const real D4D_default_epsilon;

//.cpp
const real D4D_default_epsilon(0.000000001f);//or any other value you prefer.
*/