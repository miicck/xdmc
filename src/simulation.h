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

class simulation_spec
{
public:
        // Program parameters
        int pid               = 0;      // The MPI PID of this process
        int np                = 1;      // The number of MPI processes
        int dimensions        = 3;      // The dimensionality of our system
        int target_population = 500;    // The number of DMC walkers per process
        int dmc_iterations    = 1000;   // The number of DMC iterations to carry out
        double tau            = 0.01;   // The DMC timestep
        double trial_energy   = 0;      // Energy used to control the DMC population

        // The system which will be copied to generate walkers
        std::vector<particle*> template_system;

        // The potentials applied to the system (additive)
        std::vector<external_potential*> potentials;

        // Output files
        std::ofstream wavefunction_file;
        std::ofstream evolution_file;
	std::ofstream progress_file;

	// Loads system from input, initializes MPI, opens output files etc.
        void load(int argc, char** argv);

	// Closes output files and frees template_system and potentials
	void free_memory();

private:
	
	// Reads/parses the input file
	void read_input();
	
	// Parse an atom from an input line split by whitespace
	void parse_atom(std::vector<std::string> split);
};

extern simulation_spec simulation;

#endif


