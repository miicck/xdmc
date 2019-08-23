#ifndef __MATH__
#define __MATH__

inline int sign(double val)
{
    return (0 < val) - (val < 0);
}

inline double fexp(double x)
{
    return exp(x);
}

#endif
