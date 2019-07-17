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

// Class to help with I/O. Leaves the file
// closed until it's needed.
class output_file
{
public:
	output_file() { filename = "/dev/null"; }
	output_file(std::string fn) { filename = fn; }

	template<class T>
	output_file& operator<<(T t)
	{
		if (!file.is_open()) file.open(filename);
		file << t;
		if (auto_flush) flush();
		return (*this);
	}

	void open(std::string fn) { filename = fn; }
	void close() { if(file.is_open()) file.close(); }
	void flush() { if(file.is_open()) file.flush(); }

	// Set to true to automatically flush the file after each write
	bool auto_flush = false; 

private:
	std::string filename;
	std::ofstream file;
};

// This class represents a complete specification for
// a simulation to run and progress, including ready-to-write
// output files and control parameters. 
class simulation_spec
{
public:
        // Program parameters
        int pid                    = 0;    // The MPI PID of this process
        int np                     = 1;    // The number of MPI processes
        int dimensions             = 3;    // The dimensionality of our system
        int target_population      = 500;  // The number of DMC walkers per process
        int dmc_iterations         = 1000; // The number of DMC iterations to carry out
        double tau                 = 0.01; // The DMC timestep
	double tau_c_ratio         = 1.0;  // The ratio for tau_c to tau
        double trial_energy        = 0;    // Energy used to control the DMC population
	bool write_wavefunction    = true; // True if we should write wavefunction files
	bool exchange_moves        = true; // True if we should make exchange moves
	bool make_cancellations    = true; // True if we should make walker cancellations
	bool particle_interactions = true; // True if particle-particle interactions are on

        // The external potentials applied to the system (additive)
        std::vector<external_potential*> potentials;

        // The system which will be copied to generate walkers
        std::vector<particle*> template_system;

	// Exchange information about the system
	std::vector<int> exchange_pairs;
	std::vector<int> exchange_values;
	int fermionic_exchange_pairs = 0;
	int bosonic_exchange_pairs   = 0;

	// Derived information about the system
	double total_charge();

        // Output files
        output_file wavefunction_file;
        output_file evolution_file;
	output_file progress_file;
	output_file error_file;

	// Flush output files so we have information if a run terminates
	void flush();

	// Output details of the simulation to the progress file
	void output_sim_details();

	// Loads system from input, initializes MPI, opens output files etc.
        void load(int argc, char** argv);

	// Closes output files and frees template_system and potentials
	void free_memory();

	// Get the time since startup
	double time();

private:
	
	// Reads/parses the input file
	void read_input();

	// Parse a particle from an input line
	void parse_particle(std::vector<std::string> split);

	// Parse an atomic potential from an input file
	void parse_atomic_potential(std::vector<std::string> split);

	// The result of clock() called at load
	int start_clock;
};

extern simulation_spec simulation;

#endif


