/*
    DMCPP
    Copyright (C) 2019 Michael Hutcheon (email mjh261@cam.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License see <https://www.gnu.org/licenses/>.
*/

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


