/*

    XDMC
    Copyright (C) Michael Hutcheon (email mjh261@cam.ac.uk)

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
#include <mpi.h>

#include "catch.h"
#include "constants.h"
#include "random.h"
#include "dmc_math.h"
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

walker* walker :: mpi_copy(walker* to_copy, int root_pid)
{
    // Create a copy of a given walker across mpi 
    // processes (i.e copy a walker from the given
    // root process).
    //
    // to_copy 
    //         = the walker to copy (on the root process)
    //         = nullptr (all other processes)

    walker* copy;
    if (params::pid == root_pid)
        // Simply copy the walker
        copy = to_copy->copy();
    else
        // Create a new walker
        copy = new walker();

    // Distribute the walker coordinates across processes
    for (unsigned i=0; i<copy->particles.size(); ++i)
        MPI_Bcast(copy->particles[i]->coords, params::dimensions, MPI_DOUBLE, root_pid, MPI_COMM_WORLD);

    // Distribute the walker weight across processes
    MPI_Bcast(&copy->weight, 1, MPI_DOUBLE, root_pid, MPI_COMM_WORLD);
    return copy;
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

int walker :: reflect_to_irreducible()
{
    // Reflect this walker using exchange symmetry until
    // we're in the irreducible section of configuration space
    // doesn't change the weight at all (even for fermionic exchanges)
    // Returns the number of exchanges made

    int exchanges_made = 0;

    while(true)
    {
        bool swap_made = false;

        for (unsigned i=0; i<params::exchange_groups.size(); ++i)
        {
            exchange_group* group = params::exchange_groups[i];
            for (unsigned j=1; j<group->particles.size(); ++j)
            {
                particle* p1 = particles[group->particles[j-1]];
                particle* p2 = particles[group->particles[j]];

                // If particles in the exchange group
                // are in the wrong order, swap them
                if (!p1->in_coord_order(p2))
                {
                    particles[group->particles[j-1]] = p2;
                    particles[group->particles[j]]   = p1;
                    swap_made = true;
                    ++exchanges_made;
                }
            }
        }

        // No swaps => we're in the irreducible wedge
        if (!swap_made)
            break;
    }

    return exchanges_made;
}

bool walker :: is_irreducible()
{
    // Returns true if the walker is in the 
    // irreducible wedge of configuration space.
    for (unsigned i=0; i<params::exchange_groups.size(); ++i)
    {
        exchange_group* group = params::exchange_groups[i];
        for (unsigned j=1; j<group->particles.size(); ++j)
        {
            particle* p1 = particles[group->particles[j-1]];
            particle* p2 = particles[group->particles[j]];

            // If particles in the exchange group
            // are in the wrong order => we're not
            // in the irreducible wedge.
            if (!p1->in_coord_order(p2))
                return false;
        }
    }
    return true;
}

bool walker :: crossed_nodal_surface(walker* other)
{
    // Returns true if, to get to the other walker
    // we must cross a nodal surface

    // Cant tell in > 1D
    if (params::dimensions > 1)
        throw "Can't tell if a walker crossed the nodal surface in > 1D!";

    for (unsigned n=0; n < params::exchange_groups.size(); ++n)
    {
        exchange_group* eg = params::exchange_groups[n];

        // Only fermionic exchanges define nodes
        if (eg->sign >= 0) continue;

        // Run over exchange pairs
        for (unsigned m=0; m<eg->pairs.size(); ++m)
        {
            unsigned i = eg->pairs[m].first;
            unsigned j = eg->pairs[m].second;

            // Check if this pair has crossed it's nodal surface
            particle* pi1 = particles[i];
            particle* pj1 = particles[j];
            particle* pi2 = other->particles[i];
            particle* pj2 = other->particles[j];

            double d1 = pi1->coords[0] - pj1->coords[0];
            double d2 = pi2->coords[0] - pj2->coords[0];

            if (sign(d1) != sign(d2))
                return true;
        }
    }

    return false;
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

    for(unsigned n=0; n<params::exchange_groups.size(); ++n)
    {
        // Make an exchange with probability exchange_prob
        if (rand_uniform() > params::exchange_prob) continue;
        exchange_group* eg = params::exchange_groups[n];

        if (params::full_exchange)
        {
            // Pick a random permutation
            unsigned i = rand() % eg->perms->size();

            // Record where the old particles were
            std::vector<particle*> old_particles;
            for (unsigned j=0; j<particles.size(); ++j)
                old_particles.push_back(particles[j]);

            // Put them into their permuted positions
            for (unsigned j=0; j<eg->perms->elements(); ++j)
            {
                unsigned k_unperm = (*eg->perms)[0][j];
                unsigned k_perm   = (*eg->perms)[i][j];
                particles[k_perm] = old_particles[k_unperm];
            }

            // Update the weight according to the sign of the permutation
            this->weight *= eg->weight_mult(i);
        }
        else
        {
            // Pick a random exchangable pair (i.e exchange operator)
            unsigned i = rand() % eg->pairs.size();
            particle* p1 = particles[eg->pairs[i].first];
            particle* p2 = particles[eg->pairs[i].second];

            // Exchange them 
            this->weight *= double(eg->sign);
            p1->exchange(p2);
         }
    }
}

void walker :: change_sign()
{
    // Pick a random exchange group with negative sign
    unsigned group = rand() % params::exchange_groups.size();
    exchange_group* eg = params::exchange_groups[group];
    if (eg->sign >= 0) throw "Not implemented!";

    // Pick a random pair of particles from that exchange group
    unsigned pair = rand() % eg->pairs.size();
    particle* p1  = particles[eg->pairs[pair].first];
    particle* p2  = particles[eg->pairs[pair].second];

    // Make the corresponding exchange move
    this->weight *= double(eg->sign);
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

double walker :: diffusive_greens_function(walker* other, double tau)
{
    // Evaluate the diffusive greens function of this walker
    // at the configuration of the other walker
    double r2 = this->sq_distance_to(other);
    return fexp(-r2/(2*tau))/sqrt(2*PI*tau);
}

double* walker :: exchange_diffusive_gf(walker* other, double tau)
{
    // Evaluate the exchange-diffusive greens function
    // sum_{P_i} G(other, P_i this, tau) \sign(P_i).
    // Note this will ignore the weight of other, and
    // use only the configuration.
    double* ret = new double[2];
    ret[0] = 0;
    ret[1] = 0;
    for (unsigned n=0; n<params::exchange_groups.size(); ++n)
    {
        exchange_group* eg = params::exchange_groups[n];

        double r2_unpermuted = 0;
        for (unsigned i=0; i<particles.size(); ++i)
        {
            // Check if the i^th particle is
            // permuted by this group 
            bool is_permuted = false;
            for (unsigned j=0; j<eg->perms->elements(); ++j)
                if ((*eg->perms)[0][j] == i)
                {
                    is_permuted = true;
                    break;
                }
            if (is_permuted)
                continue;

            // Sum r^2 for particles that are not permuted
            r2_unpermuted += particles[i]->sq_distance_to(other->particles[i]);
        }

        // Loop over permutations
        unsigned* unperm = (*eg->perms)[0];
        for (unsigned m=0; m<eg->perms->size(); ++m)
        {
            unsigned* perm = (*eg->perms)[m];
            double r2      = r2_unpermuted;

            // Add the contribution to r^2 from the permuted particles
            for (unsigned i=0; i<eg->perms->elements(); ++i)
            {
                unsigned j = perm[i];
                unsigned k = unperm[i];
                r2   += particles[j]->sq_distance_to(other->particles[k]);
            }

            // Sum the greens functions
            double gf = fexp(-r2/(2*tau))/sqrt(2*PI*tau);
            if (eg->weight_mult(m) > 0) ret[0] += gf;
            else ret[1] += gf;
        }
    }
    return ret;
}

void walker :: write_coords(output_file& file)
{
    // Write the walker wavefunction in the form
    // [weight: x1, y1, z1 ...; x2, y2, z2 ...; ...]
    // where x1 is the x coord of the first particle etc
    file << this->weight << ":";
    for (unsigned i=0; i<particles.size(); ++i)
    {
        for (unsigned j=0; j<params::dimensions; ++j)
        {
            file << particles[i]->coords[j];
            if (j != params::dimensions - 1)
                file << ",";
        }
        if (i != particles.size() - 1)
            file << ";";
    }
    file << "\n";
}

unsigned walker :: particle_count()
{
    // Return the number of particles
    return particles.size();
}

bool walker :: compare(walker* other)
{
    // Returns false if these walkers are not identical
    // in some way (for testing purposes)
    if (this->sq_distance_to(other) != 0) return false;
    if (this->weight != other->weight) return false;
    if (this->particle_count() != other->particle_count()) return false;
    return true;
}

TEST_CASE("Basic walker tests", "[walker]")
{
    // Create some walkers to test
    walker* w1 = new walker();
    walker* w2 = new walker();

    SECTION("Copy method")
    {
        walker* w = w1->copy();
        REQUIRE(w->compare(w1));
        delete w;
    }

    SECTION("MPI copy method")
    {
        walker* w = params::pid == 0 ? w1 : nullptr;
        walker* wc = walker::mpi_copy(w, 0);
        REQUIRE(wc->compare(w1));
    }

    delete w1;
    delete w2;
}

TEST_CASE("Irreducible wedge", "[walker]")
{
    SECTION("Reflect to irreducible")
    {
        // Test several random configs
        for (int n=0; n<10; ++n)
        {
            // Create a walker to test
            walker* w1 = new walker();
            
            // Put it in a pseudo-random config
            w1->diffuse(1.0);

            // Check reflect_to_irreducible works
            w1->reflect_to_irreducible();
            REQUIRE(w1->is_irreducible());

            delete w1;
        }
    }
}
