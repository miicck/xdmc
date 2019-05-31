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

#include "simulation.h"
#include "particle.h"
#include "walker.h"

simulation_spec simulation;

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

// Parse an atom from the input file. Creates a
// nucleus and the appropriate number of electrons
// in the template system.
void simulation_spec :: parse_atom(std::vector<std::string> split)
{
        double charge = std::stod(split[1]);
        int electrons = std::stoi(split[2]);
        double mass   = std::stod(split[3]);

	double* centre = new double[simulation.dimensions];
	for (int i=4; i<4+dimensions; ++i)
		centre[i] = std::stod(split[i]);

	// Create the nucleus
	potentials.push_back(new atomic_potential(charge, centre));

        // Add the specified number of electrons
        for (int i=0; i<electrons; ++i)
        {
                auto e = new electron();
                for (int i=0; i<dimensions; ++i)
                        e->coords[i] = centre[i];
                template_system.push_back(e);
        }
}

// Parse the input file.
void simulation_spec :: read_input()
{
        std::ifstream input("input");
        for (std::string line; getline(input, line); )
        {
                auto split = split_whitespace(line);
                if (split.size() == 0) continue;
                std::string tag = split[0];

                // Ignore comments
                if (tag.rfind("!" , 0) == 0) continue;
                if (tag.rfind("#" , 0) == 0) continue;
                if (tag.rfind("//", 0) == 0) continue;

                // Read in the dimensionality
                if (tag == "dimensions")
                        dimensions = std::stoi(split[1]);

                // Read in the number of DMC walkers and convert
                // to walkers-per-process
                else if (tag == "walkers")
                        target_population = std::stoi(split[1])/np;

                // Read in the number of DMC iterations
                else if (tag == "iterations")
                        dmc_iterations = std::stoi(split[1]);

                // Read in the DMC timestep
                else if (tag == "tau")
                        tau = std::stod(split[1]);

		// Turn off exhange moves
		else if (tag == "no_exchange")
			exchange_moves = false;

		// Turn off cancellation of walker weights
		else if (tag == "no_cancellation")
			cancellation = false;

                // Read in an atom
                else if (tag == "atom")
                        parse_atom(split);

                // Read in an electron
                else if (tag == "electron")
                        template_system.push_back(new electron());

                // Adds a noninteracting fermion into the system
                else if (tag == "nif")
                        template_system.push_back(new non_interacting_fermion());

                // Add a harmonic well to the system
                else if (tag == "harmonic_well")
                        potentials.push_back(new harmonic_well(std::stod(split[1])));
        }
        input.close();
}

void simulation_spec::load(int argc, char** argv)
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

        // Read our input and setup parameters accordingly 
        // do for each process sequentially to avoid access issues
        for (int pid_read = 0; pid_read < np; ++ pid_read)
        {
                if (pid == pid_read) read_input();
                MPI_Barrier(MPI_COMM_WORLD);
        }

	// Open various output files
	error_file.open(       "error_"        + std::to_string(pid));
	progress_file.open(    "progress_"     + std::to_string(pid));
	evolution_file.open(   "evolution_"    + std::to_string(pid));
	wavefunction_file.open("wavefunction_" + std::to_string(pid));
}

void simulation_spec::free_memory()
{
	// Close various output files
	progress_file.close();
	evolution_file.close();
	wavefunction_file.close();

	// Free memory in template_system
	for (int i=0; i<template_system.size(); ++i)
		delete template_system[i];

	// Free memory in potentials
	for (int i=0; i<potentials.size(); ++i)
		delete potentials[i];

	// Output info on objects that werent deconstructed properly
        if (walker::constructed_count != 0 || particle::constructed_count != 0)
        error_file << "PID: " << pid << " un-deleted objects:\n"
                   << "  Walkers   : " << walker::constructed_count   << "\n"
                   << "  Particles : " << particle::constructed_count << "\n";

	error_file.close();
	MPI_Finalize();
}

void simulation_spec :: flush()
{
	// Flush all of the output files to disk
	error_file.flush();
	progress_file.flush();
	evolution_file.flush();
	wavefunction_file.flush();
}

double simulation_spec :: time()
{
	// Return the time in seconds since startup
	return double(clock()-start_clock)/double(CLOCKS_PER_SEC);
}
