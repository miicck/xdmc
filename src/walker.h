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

#ifndef __WALKER__
#define __WALKER__

#include <vector>
#include <string>
#include "particle.h"
#include "params.h"

// The object used by the diffusion monte carlo algorithm
// to represent a snapshot of the system.
class walker
{
public:
    walker();
    ~walker();
    static int constructed_count;

    double weight = 1.0;
    unsigned particle_count();

    double potential();
    double sq_distance_to(walker* other);
    double diffusive_greens_function(walker* other, double tau=params::tau);
    double* exchange_diffusive_gf(walker* other, double tau=params::tau);

    bool crossed_nodal_surface(walker* other);
    bool compare(walker* other);

    void diffuse(double tau);
    void exchange();
    void change_sign();
    int reflect_to_irreducible();
    bool is_irreducible();

    walker* copy();
    walker* branch_copy();
    static walker* mpi_copy(walker* to_copy, int root_pid);

    void write_coords(output_file& file);
    std::string summary();
private:

    // The particles in this system snapshot
    std::vector<particle*> particles;

    // Create a walker from a given particle set 
    walker(std::vector<particle*> particles);

    // True if the potential needs re-evaluating
    bool potential_dirty = true;
    double last_potential = 0;
};

#endif

