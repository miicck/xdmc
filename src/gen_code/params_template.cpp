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
std::vector<int> params::exchange_pairs;
std::vector<int> params::exchange_values;
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
    for (int i=0; i<dimensions; ++i)
        p->coords[i] = std::stod(split[5+i]);
    template_system.push_back(p);
}

void parse_atomic_potential(std::vector<std::string> split)
{
    // Parse an atomic potential from the format
    // atomic_potential charge x1 x2 x3 ...
    
    double charge  = std::stod(split[1]);
    double* coords = new double[dimensions];
    for (int i=0; i<dimensions; ++i)
        coords[i] = std::stod(split[2+i]);
    potentials.push_back(new atomic_potential(charge, coords));
}

// Parse the input file.
void read_input()
{
    std::ifstream input("input");
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

    // Work out exchange properties of the system
    for (unsigned i=0; i<template_system.size(); ++i)
    {
        particle* p1 = template_system[i];
        for (unsigned j=0; j<i; ++j)
        {
            particle* p2 = template_system[j];
            int exchange_sym = p1->exchange_symmetry(p2);

            // These particles cannot be exchanged 
            if (exchange_sym == 0) continue;

            // Record this exchangable pair
            exchange_pairs.push_back(j);
            exchange_pairs.push_back(i);
            exchange_values.push_back(exchange_sym);
        }
    }
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

    // Output a summary of exchange information
    progress_file << "Exchange pairs (sign, particle 1, particle 2)\n";
    for (unsigned i=0; i<exchange_values.size(); ++i)
           progress_file << "    " << exchange_values[i] << " "
                 << exchange_pairs[2*i]   << " "
                 << exchange_pairs[2*i+1] << "\n";

    progress_file << "\n";
}

void params::load(int argc, char** argv)
{
    // Get the start time so we can time stuff
    start_clock = clock();
    
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

    // Read our input and setup parameters accordingly 
    // do for each process sequentially to avoid access issues
    for (int pid_read = 0; pid_read < np; ++ pid_read)
    {
            if (pid == pid_read) read_input();
            MPI_Barrier(MPI_COMM_WORLD);
    }
    
    // Output parameters to the progress file
    output_sim_details();
}

void params::free_memory()
{
    // Close various output files
    progress_file.close();
    evolution_file.close();
    wavefunction_file.close();

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
