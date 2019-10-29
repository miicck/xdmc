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

#ifndef __SIMULATION__
#define __SIMULATION__

#include <vector>
#include <fstream>

#include "particle.h"
#include "potential.h"
#include "output_file.h"
#include "dmc_math.h"

// This represents a group of particle indicies
// that can be exchanged with one another
class exchange_group
{
public:
    exchange_group();
    ~exchange_group();

    int sign;
    int weight_mult(unsigned permutation);
    std::vector<unsigned> particles;
    std::vector<std::pair<unsigned,unsigned>> pairs;
    permutations<unsigned>* perms;

    void add(unsigned index);
    void finalize();
    std::string one_line_summary();
};

// This namespace represents a complete specification for
// a simulation to run and progress, including ready-to-write
// output files and control parameters. 
namespace params
{
    //%%%%%%%%%%%%%%%%%%%%//
    // PROGRAM PARAMETERS //
    //%%%%%%%%%%%%%%%%%%%%//

    PYTHON_GEN_PARAMS_HERE

    // The external potentials applied to the system (additive)
    extern std::vector<external_potential*> potentials;

    // The system which will be copied to generate walkers
    // and exchange information about it
    extern std::vector<particle*> template_system;
    extern std::vector<exchange_group*> exchange_groups;

    // Output files
    extern output_file nodal_surface_file;
    extern output_file wavefunction_file;
    extern output_file evolution_file;
    extern output_file progress_file;
    extern output_file error_file;

    //%%%%%%%%%%%//
    // FUNCTIONS //
    //%%%%%%%%%%%//

    // Loads system from input, initializes MPI, opens output files etc.
    void load(int argc, char** argv);

    // Closes output files and frees template_system and potentials
    void free_memory();

    // Flush output files so we have information if a run terminates
    void flush();

    // Get the time since startup
    double time();

    // Get the time DMC iterations have been running
    double dmc_time();
};

#endif


