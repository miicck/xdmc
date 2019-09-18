#ifndef __MATH__
#define __MATH__

#include "params.h"

inline int sign(double val)
{
    return (0 < val) - (val < 0);
}

inline double fexp(double x)
{
    return exp(x);
}

inline double coulomb(double q1, double q2, double r)
{
    return q1*q2/(r+params::coulomb_softening);
}

#endif
