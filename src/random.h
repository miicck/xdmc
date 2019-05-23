#ifndef __RANDOM__
#define __RANDOM__

#include <cmath>

// Generate a uniform random number \in [0,1].
inline double rand_uniform()
{
        return (rand()%RAND_MAX)/double(RAND_MAX);
}

// Generate a zero-mean noramlly distributed number
// with the specified variance using a Box-Muller transform.
inline double rand_normal(double var)
{
        double u1 = rand_uniform();
        double u2 = rand_uniform();
        return sqrt(-2*var*log(u1)) * sin(2*PI*u2);
}

#endif
