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

#include <mpi.h>
#include <vector>
#include <sstream>
#include <iterator>
#include <iostream>

#include "params.h"
#include "particle.h"
#include "walker.h"

using namespace params;

// The value of clock() at startup
int start_clock;

PYTHON_GEN_PARAMS_HERE

// Global param:: variables
std::vector<external_potential*> params::potentials;
std::vector<particle*>           params::template_system;
std::vector<exchange_group*>     params::exchange_groups;
output_file params::nodal_surface_file;
output_file params::wavefunction_file;
output_file params::evolution_file;
output_file params::progress_file;
output_file params::error_file;

// Split a string on whitespace
std::vector<std::string> split_whitespace(std::string to_split)
{
        // Iterator magic
        std::istringstream buffer(to_split);
        std::vector<std::string> split;
        std::copy(std::istream_iterator<std::string>(buffer),
                  std::istream_iterator<std::string>(),
                  std::back_inserter(split));
        return split;
}

void parse_particle(std::vector<std::string> split)
{
    // Parse a particle from the format
    // particle name mass charge half_spins x1 x2 x3 ...
    particle* p   = new particle();
    p->name       = split[1];
    p->mass       = std::stod(split[2]);
    p->charge     = std::stod(split[3]);
    p->half_spins = std::stoi(split[4]);
    for (unsigned i=0; i<dimensions; ++i)
        p->coords[i] = std::stod(split[5+i]);
    template_system.push_back(p);
}

void parse_atomic_potential(std::vector<std::string> split)
{
    // Parse an atomic potential from the format
    // atomic_potential charge x1 x2 x3 ...
    
    double charge  = std::stod(split[1]);
    double* coords = new double[dimensions];
    for (unsigned i=0; i<dimensions; ++i)
        coords[i] = std::stod(split[2+i]);

    // Create an atomic potential
    potentials.push_back(new atomic_potential(charge, coords));
}

// Check that the parameter set is sensible
bool check_params()
{
    if (params::template_system.size() == 0)
    {
        params::error_file << "Error: no particles found in system!\n";
        return false;
    }

    return true;
}

// Parse the input file.
bool read_input()
{
    // Open the input file
    std::ifstream input("input");
    if (!input.is_open())
    {
        params::error_file << "Error: could not read input file!\n";
        return false;
    }

    for (std::string line; getline(input, line); )
    {
        // Ignore comments
        if (line.rfind("!" , 0) == 0) continue;
        if (line.rfind("#" , 0) == 0) continue;
        if (line.rfind("//", 0) == 0) continue;

        // Split line into words
        std::vector<std::string> split = split_whitespace(line);
        if (split.size() == 0) continue;
        std::string tag = split[0];

        PYTHON_PARSE_PARAMS_HERE

        // Parse a particle input line
        if (tag == "particle")
            parse_particle(split);
    
        // Add a potential from a grid to the system
        else if (tag == "grid_potential")
            potentials.push_back(new grid_potential(split[1]));

        // Add a harmonic well to the system
        else if (tag == "harmonic_well")
            potentials.push_back(new harmonic_well(std::stod(split[1])));

        // Add an atomic potential to the system
        else if (tag == "atomic_potential")
            parse_atomic_potential(split);

        // Record misunderstandings
        else
            error_file << "Did not know what to do with input line '" 
                       << line << "'\n";
    }

    // Divide target population across processes
    params::target_population /= params::np;

    input.close();

    // Record exchange groups within the system
    bool in_group[template_system.size()];
    for (unsigned i=0; i<template_system.size(); ++i)
        in_group[i] = false;

    for (unsigned i=0; i<template_system.size(); ++i)
    {
        // Already in a group
        if (in_group[i]) continue; 

        particle* p1 = template_system[i];

        // Add p1 to its own group
        exchange_group* eg = new exchange_group();
        eg->add(i);
        in_group[i] = true;

        for (unsigned j=i+1; j<template_system.size(); ++j)
        {
            particle* p2 = template_system[j];
            if (p1->exchange_symmetry(p2) != 0)
            {
                // Add p2 to p1's group
                eg->add(j);
                in_group[j] = true;
            }
        }

        if (eg->particles.size() < 2)
        {
            // It's not an exchange group if there is only one!
            delete eg;
        }
        else
        {
            // Record this exchange group
            eg->finalize();
            params::exchange_groups.push_back(eg);
        }
    }

    // Check the parameter set
    return check_params();
}

// Construct/destruct an exchange_group
exchange_group ::  exchange_group() { perms = nullptr; }
exchange_group :: ~exchange_group() { if (perms != nullptr) delete perms; }

// Add a particle index to an exchange group
void exchange_group :: add(unsigned index) { particles.push_back(index); }

int exchange_group :: weight_mult(unsigned perm_index)
{
    // Bosonic exhcnage => weight stays the same
    if (this->sign == 1)
        return 1;
    
    // Odd fermionic permutations => weight -> weight * -1
    if (this->perms->sign(perm_index) == -1)
        return -1;
    
    // Even fermionic permutations => weight stays the same
    return 1;
}

void exchange_group :: finalize()
{
    // Work out exchange group permutations
    perms = new permutations<unsigned>(particles);

    // Work out exchange group sign and double check that it is consistent
    // across the group. Construct the exchange pairs
    sign = -2;
    for (unsigned i=0; i<particles.size(); ++i)
        for (unsigned j=0; j<i; ++j)
        {
            unsigned n = particles[i];
            unsigned m = particles[j];

            // Record the pair
            std::pair<int,int> pair(m, n);
            pairs.push_back(pair);

            // Check/record the sign
            particle* p1 = params::template_system[n];
            particle* p2 = params::template_system[m];
            if (sign < -1) sign = p1->exchange_symmetry(p2); 
            else if (sign != p1->exchange_symmetry(p2))
                throw "Error, inconsistent signs in exchange group!";
        }
}

std::string exchange_group :: one_line_summary()
{
    // Return a simple description of the exchange group
    std::stringstream s;
    s << "Sign = " << sign
      << ", permutations = " << perms->size() 
      << ", particles = {";
    for (unsigned i=0; i<particles.size(); ++i)
        s << particles[i] << " ";
    s << "\b}";
    return s.str();
}

void output_sim_details()
{
    // Output a summary of the parameter values used
    progress_file << "Parameter values\n";
    PYTHON_OUTPUT_PARAMS_HERE

    // Output a summary of potentials to the progress file
    progress_file << "Potentials\n";
    for (unsigned i=0; i<potentials.size(); ++i)
        progress_file << "    " << potentials[i]->one_line_description() << "\n";

    // Output a summary of particles to the progress file
    progress_file << "Particles\n";
    for (unsigned i=0; i<template_system.size(); ++i)
           progress_file << "    " << i << ": "
                 << template_system[i]->one_line_description() << "\n";

    // Output a summary of the exchange groups
    progress_file << "Exchange groups:\n";
    for (unsigned i=0; i<exchange_groups.size(); ++i)
    {
        exchange_group* eg = exchange_groups[i];
        progress_file << "    " << eg->one_line_summary() << "\n";
        for (unsigned j=0; j<eg->pairs.size(); ++j)
            progress_file << "        " 
                          << eg->pairs[j].first << " <--> " 
                          << eg->pairs[j].second << "\n";

        progress_file << "    Permutations:\n";
        for (unsigned j=0; j<eg->perms->size(); ++j)
        {
            progress_file << "        p = {";
            for (unsigned k=0; k<eg->perms->elements(); ++k)
                progress_file << (*eg->perms)[j][k] << " ";
            progress_file << "\b} " << " sign = " << eg->perms->sign(j) << "\n";
        }
    }

    progress_file << "\n";
}

bool params::load(int argc, char** argv)
{
    // Initialize mpi
    if (MPI_Init(&argc, &argv) != 0) exit(MPI_ERROR);

    // Get the number of processes and my id within them
    if (MPI_Comm_size(MPI_COMM_WORLD, &np)  != 0) exit(MPI_ERROR);
    if (MPI_Comm_rank(MPI_COMM_WORLD, &pid) != 0) exit(MPI_ERROR);

    // Seed random number generator
    srand(pid*clock());

    // Open various output files
    if (pid == 0)
    {
        // Files on the root process
        progress_file.open("progress");
        evolution_file.open("evolution");
    }

    // Files on all processes have their pid appended
    error_file.open("error_"+std::to_string(pid));
    error_file.auto_flush = true;
    wavefunction_file.open("wavefunction_"+std::to_string(pid));
    nodal_surface_file.open("nodal_surface_"+std::to_string(pid));

    // Read our input and setup parameters accordingly 
    bool input_success = read_input();

    // Check all processes succeeded
    bool all_success = false;
    MPI_Allreduce(&input_success, &all_success, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    if (!all_success)
    {
        progress_file << "Errors occured whilst reading input, stopping.\n";
        return false;
    }

    // Output parameters to the progress file
    output_sim_details();

    // Load successful
    return true;
}

void params::print_usage_info()
{
    // Print parameter usage info
    std::cout << "####################\n";
    std::cout << "# Input parameters #\n";
    std::cout << "####################\n";

    PYTHON_GEN_USAGE_INFO_HERE
}

void params::free_memory()
{
    // Close various output files
    progress_file.close();
    evolution_file.close();
    wavefunction_file.close();

    // Free memory used in exchange groups 
    for (unsigned i=0; i<exchange_groups.size(); ++i)
        delete exchange_groups[i];

    // Free memory in template_system
    for (unsigned i=0; i<template_system.size(); ++i)
        delete template_system[i];

    // Free memory in potentials
    for (unsigned i=0; i<potentials.size(); ++i)
        delete potentials[i];

    // Output info on objects that werent deconstructed properly
    if (walker::constructed_count != 0 || particle::constructed_count != 0)
    error_file << "PID: " << pid << " un-deleted objects:\n"
               << "  Walkers   : " << walker::constructed_count   << "\n"
               << "  Particles : " << particle::constructed_count << "\n";

    error_file.close();
    MPI_Finalize();
}

void params :: flush()
{
    // Flush all of the output files to disk
    error_file.flush();
    progress_file.flush();
    evolution_file.flush();
    wavefunction_file.flush();
}

double params :: time()
{
    // Return the time in seconds since startup
    return double(clock()-start_clock)/double(CLOCKS_PER_SEC);
}

double params :: dmc_time()
{
    // Return the time that the DMC algorithm has
    // been running (excluding setup time)
    return double(clock()-dmc_start_clock)/double(CLOCKS_PER_SEC);
}
