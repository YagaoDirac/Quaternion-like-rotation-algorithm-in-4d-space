#pragma once

#include "../include/D4D config.h"

real D4D_default_epsilon() { return real(0.000000001f); };
real D4D_pi() 
{
#ifdef __D4D_REAL_FLOAT__
	return std::acosf(-1);
#else
	return std::acos(-1);
#endif
};




