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

#ifndef __PARTICLE__
#define __PARTICLE__

#include <string>
#include "constants.h"

// A quantum mechanical particle, as described by its quantum numbers
// and a position which diffuses according to DMC.
class particle
{
public:
    particle();
    ~particle();
    static int constructed_count;

    double sq_distance_to(particle* other); // Returns | this->coords - other->coords |^2
    double interaction(particle* other);    // The interaction energy with some other particle
    int exchange_symmetry(particle* other); // 1, 0 or -1 depending on spin statistics
    void exchange(particle* other);         // Swap my coordinates with those of other

    void sample_wavefunction();   // Called when a request is sent to sample a walker wvfn.
    void diffuse(double tau);     // Called when a config diffuses in DMC
    particle* copy();             // Should return a (deep) copy of this particle

    double charge  = 0;           // The charge of this particle (electron charge = -1)
    double mass    = 0;           // The mass of this particle (electron mass = 1)
    int half_spins = 0;           // This spin of this particle (+/- 1 => spin = +/- 1/2)

    std::string name = "Particle";
    std::string one_line_description();

    // The location of this particle
    double* coords;
};

#endif
