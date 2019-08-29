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

#include <sstream>

#include "constants.h"
#include "random.h"
#include "math.h"
#include "walker.h"
#include "params.h"

// Used to track the number of constructred walkers
// to ensure that we delete them all again properly
int walker :: constructed_count = 0;

walker :: walker()
{
    // Default constructor:
    // Setup the walker with the particles
    // that describe the system.
    ++ constructed_count;
    this->weight = 1;
    for (unsigned i=0; i<params::template_system.size(); ++i)
        this->particles.push_back(params::template_system[i]->copy());
}

walker :: walker(std::vector<particle*> particles)
{
    // Constructor to create a walker explicitly from
    // a collection of particles
    ++ constructed_count;
    this->weight = 1;
    this->particles = particles;
}

walker :: ~walker()
{
    // Clear up memory (delete all the particles).
    -- constructed_count;
    for (unsigned i=0; i<particles.size(); ++i)
        delete(particles[i]);
}


walker* walker :: copy()
{
    // Return an exact copy of this walker
    // (copy each of the particles and the weight)
    std::vector<particle*> copied_particles;
    for (unsigned i=0; i<particles.size(); ++i)
        copied_particles.push_back(particles[i]->copy());
    walker* copy = new walker(copied_particles);
    copy->weight = this->weight;
    return copy;
}

walker* walker :: branch_copy()
{
    // Return a branched version of this walker
    // (weight is already accounted for by branching
    // step, the sign is all that continues)
    walker* bcopy = this->copy();
    bcopy->weight = sign(bcopy->weight);
    return bcopy;
}

std::string walker :: summary()
{
    // Return a summary of this walker
    std::stringstream ss;
    ss << "Weight: " << this->weight << "\n";
    for (unsigned i=0; i<particles.size(); ++i)
    {
        ss << "    Particle " << i << ":";
        for (unsigned j=0; j<params::dimensions; ++j)
            ss << particles[i]->coords[j] << " ";
        ss << "\n";
    }
    return ss.str();
}

void walker :: reflect_to_irreducible()
{
    // Reflect this walker using fermionic exchanges until
    // we're in the irreducible section of configuration space
    while(true)
    {
        bool swap_made = false;
        for (unsigned i=0; i<particles.size()-1; ++i)
        {
            particle* p1 = particles[i];
            particle* p2 = particles[i+1];
            if (p1->exchange_symmetry(p2) == 0)
                continue;

            // These particles should be swapped if they're
            // in the wrong order (sort by increasing coordinates)
            bool swap = true;
            for (unsigned j=0; j < params::dimensions; ++j)
                if (p1->coords[j] < p2->coords[j])
                {
                    // These are in the right order
                    swap = false;
                    break;
                }

            if (swap)
            {
                // Swap these particles
                particles[i]   = p2;
                particles[i+1] = p1;
                swap_made = true;
            }
        }

        // No swaps => particles in correct order
        if (!swap_made) break;
    }
}

double walker :: potential()
{
    // No need to reevaluate the potential
    if (!potential_dirty)
        return last_potential;
    
    // Evaluate the potential of the system
    // in the configuration described by this walker
    last_potential = 0;
    for (unsigned i = 0; i < particles.size(); ++i)
    {
        // Sum up external potential contributions
        for (unsigned j=0; j<params::potentials.size(); ++j)
            last_potential += params::potentials[j]->potential(particles[i]);

        // Particle-particle interactions
        // note j<i => no double counting
        for (unsigned j=0; j<i; ++j)
            last_potential += particles[i]->interaction(particles[j]);
    }

    potential_dirty = false;
    return last_potential;
}

void walker :: diffuse(double tau=params::tau)
{
    // Diffuse all of the particles
    for (unsigned i=0; i<particles.size(); ++i)
        particles[i]->diffuse(tau);
    
    // Particles have moved => potential has changed
    potential_dirty = true;
}

void walker :: exchange()
{
    // Apply random exchange moves to particles
    // Note: because only identical particles
    // are exchanged, the potential remains the
    // same => we do not need to set the potential_dirty
    // flag.

    // No exchanges possible
    if (params::exchange_values.size() == 0) return;

    // Make each type of exchange with equal probability
    //double ne = (double)params::exchange_values.size();
    //if (rand_uniform() < 1/(ne+1)) return;
    if (rand_uniform() > params::exchange_prob) return;

    // Pick a random exchangable pair
    int i = rand() % params::exchange_values.size();
    particle* p1 = particles[params::exchange_pairs[2*i]];
    particle* p2 = particles[params::exchange_pairs[2*i+1]];

    // Exchange them 
    this->weight *= double(params::exchange_values[i]);
    p1->exchange(p2);
}

double walker :: sq_distance_to(walker* other)
{
    // Return the squared distance in configuration space
    // between these two walkers: |x_this - x_other|^2 
    double r2 = 0;
    for (unsigned i=0; i<particles.size(); ++i)
        r2 += particles[i]->sq_distance_to(other->particles[i]);
    return r2;
}

double walker :: diffusive_greens_function(walker* other)
{
    // Evaluate the diffusive greens function of this walker
    // at the configuration of the other walker
    double r2 = this->sq_distance_to(other);
    return fexp(-r2/(2*params::tau));
}

void walker :: write_wavefunction()
{
    // Write the walker wavefunction in the form
    // [weight: x1, y1, z1 ...; x2, y2, z2 ...; ...]
    // where x1 is the x coord of the first particle etc
    params::wavefunction_file << this->weight << ":";
    for (unsigned i=0; i<particles.size(); ++i)
    {
        for (unsigned j=0; j<params::dimensions; ++j)
        {
            params::wavefunction_file
                << particles[i]->coords[j];
            if (j != params::dimensions - 1)
                params::wavefunction_file << ",";
        }
        if (i != particles.size() - 1)
            params::wavefunction_file << ";";
    }
}
