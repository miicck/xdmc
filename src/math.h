#ifndef __MATH__
#define __MATH__

// Distance used for finite-difference calculations
const double FD_EPS = 0.001;

inline int sign(double val)
{
	return (0 < val) - (val < 0);
}

#endif
