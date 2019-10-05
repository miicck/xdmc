#include "dmc_math.h"
#include "params.h"

int sign(double val)
{
    // Returns the sign of non-zero values
    // or 0 for exactly 0 values
    return (0 < val) - (val < 0);
}

double coulomb(double q1, double q2, double r)
{
    // The coulomb interation used throughout
    // the code, including any softening
    return q1*q2/(r+params::coulomb_softening);
}

unsigned factorial(unsigned n)
{
    // Returns n!
    if (n == 0) return 1;
    return n * factorial(n-1);
}
