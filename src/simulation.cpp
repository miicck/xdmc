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

void simulation_spec :: parse_particle(std::vector<std::string> split)
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

void simulation_spec :: parse_atomic_potential(std::vector<std::string> split)
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

		// Read in a particle
		else if (tag == "particle")
			parse_particle(split);

                // Add a harmonic well to the system
                else if (tag == "harmonic_well")
                        potentials.push_back(new harmonic_well(std::stod(split[1])));

		// Add an atomic potential to the system
		else if (tag == "atomic_potential")
			parse_atomic_potential(split);

		// Should we write the wavefunctions this run
		else if (tag == "write_wavefunction")
			write_wavefunction = split[1] == "true";

		// Turn off exchange moves
		else if (tag == "no_exchange")
			exchange_moves = false;

		// Turn off walker cancellations
		else if (tag == "no_cancellation")
			make_cancellations = false;
        }
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

			// Record the number of exchange pairs
			if (exchange_sym < 0) ++fermionic_exchange_pairs;
			else ++bosonic_exchange_pairs;

			// Record this exchangable pair
			exchange_pairs.push_back(j);
			exchange_pairs.push_back(i);
			exchange_values.push_back(exchange_sym);
		}
	}
}

void simulation_spec :: output_sim_details()
{
	// Output system information
	progress_file << "System loaded\n";
	progress_file << "    Dimensions         : " << dimensions               << "\n";
	progress_file << "    Particles          : " << template_system.size()   << "\n";
	progress_file << "    Exchange pairs     : " << exchange_values.size()   << "\n"; 
	progress_file << "       Bosonic         : " << bosonic_exchange_pairs   << "\n";
	progress_file << "       Fermionic       : " << fermionic_exchange_pairs << "\n";
	progress_file << "    Total charge       : " << total_charge()           << "\n";
	progress_file << "    DMC timestep       : " << tau                      << "\n";

	progress_file << "    DMC walkers        : " 
		      << target_population*np << " (total) "
		      << target_population    << " (per process)\n";

	progress_file << "    DMC iterations     : " 
		      << dmc_iterations << " => Imaginary time in [0, " 
		      << dmc_iterations*tau << "]\n";

	progress_file << "    MPI processes      : " << np  << "\n";
	progress_file << "    Write wavefunction : " << write_wavefunction << "\n";
	progress_file << "    Exchange moves     : " << exchange_moves << "\n";
	progress_file << "    Cancellations      : " << make_cancellations << "\n";

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
	if (pid == 0)
	{
		// Files on the root process
		progress_file.open("progress");
		evolution_file.open("evolution");
	}

	// Files on all processes have their pid appended
	error_file.open("error_"+std::to_string(pid));
	wavefunction_file.open("wavefunction_"+std::to_string(pid));
	
	// Output parameters to the progress file
	output_sim_details();
}

double simulation_spec :: total_charge()
{
	// Work out the total charge on the system
	double sum = 0;
	for (unsigned i=0; i<template_system.size(); ++i)
		sum += template_system[i]->charge;
	return sum;
}

void simulation_spec::free_memory()
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
