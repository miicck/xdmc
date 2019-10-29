/*
    XDMC
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

#include <string>
#include <sstream>

#include "params.h"
#include "particle.h"
#include "random.h"
#include "dmc_math.h"

int particle::constructed_count;

particle::particle()
{
    // Track the number of particles that have been constructed
    ++ constructed_count;
    coords = new double[params::dimensions];
}

particle::~particle()
{
    // TRack the number of particles that have been deconstructed
    -- constructed_count;
    delete[] coords;
}

particle* particle :: copy()
{
    // Create a copy of this particle by copying each of its
    // quantum numbers and it's location
    particle* p = new particle();

    p->half_spins = this->half_spins;
    p->charge     = this->charge;
    p->mass       = this->mass;

    for (unsigned i=0; i<params::dimensions; ++i)
        p->coords[i] = this->coords[i];

    return p;
}

std::string particle :: one_line_description()
{
    // Returns a one line description of the particle
    // for writing to output files
    std::stringstream des;
    des << name;
    des << " [charge: " << this->charge 
        << " mass: " << this->mass;
    if (this->half_spins % 2 == 0)
        des << " spin: " << this->half_spins/2 << "]";
    else
        des << " spin: " << this->half_spins << "/2]";
    return des.str();
}

int particle :: exchange_symmetry(particle* other)
{
    // Returns 
    //     0 for non-identical particles
    //     1 for identical bosons
    //    -1 for identical fermions

    // Check if identical (same quantum numbers)
    if (this->half_spins != other->half_spins)      return 0;
    if (fabs(this->mass   - other->mass)   > 10e-4) return 0;
    if (fabs(this->charge - other->charge) > 10e-4) return 0;

    // Return exchange symmetry
    if (this->half_spins % 2 == 0) return 1; // Boson
    return -1; // Fermion
}

void particle :: exchange(particle* other)
{
    // Swap the coordinates of this with other
    for (unsigned i=0; i<params::dimensions; ++i)
    {
        double tmp = this->coords[i];
        this->coords[i]  = other->coords[i];
        other->coords[i] = tmp;
    }
}

double particle :: sq_distance_to(particle* other)
{
    // Unpacking for speed
    // (about twice as fast by my tests)
    if (params::dimensions == 1)
    {
        double dx = this->coords[0] - other->coords[0];
        return dx*dx;
    }
    else if (params::dimensions == 2)
    {
        double dx0 = this->coords[0] - other->coords[0];
        double dx1 = this->coords[1] - other->coords[1];
        return dx0*dx0 + dx1*dx1;
    }
    else if (params::dimensions == 3)
    {
        double dx0 = this->coords[0] - other->coords[0];
        double dx1 = this->coords[1] - other->coords[1];
        double dx2 = this->coords[2] - other->coords[2];
        return dx0*dx0 + dx1*dx1 + dx2*dx2;
    }

    // Returns | this->coords - other->coords |^2
    double r2 = 0;
    for (unsigned i=0; i<params::dimensions; ++i)
    {
        double dxi = this->coords[i] - other->coords[i];
        r2 += dxi * dxi;
    }
    return r2;
}

double particle::interaction(particle* other)
{
    if (fabs(this->charge)  < 10e-10) return 0;
    if (fabs(other->charge) < 10e-10) return 0;

    double r = sqrt(this->sq_distance_to(other));
    return coulomb(this->charge, other->charge, r);
}

void particle :: diffuse(double tau)
{
    // Diffuse the particle by moving each
    // coordinate by an amount sampled from
    // a normal distribution with variance tau.
    for (unsigned i=0; i<params::dimensions; ++i)
        this->coords[i] += rand_normal(tau);
}

void particle :: sample_wavefunction()
{
    // Output my coordinates in the form x1,x2, ... xn
    // (with no trailing comma)
    for (unsigned i=0; i<params::dimensions; ++i)
    {
        params::wavefunction_file << this->coords[i];
        if (i < params::dimensions - 1)
            params::wavefunction_file << ",";
    }
    params::wavefunction_file << ";";
}

