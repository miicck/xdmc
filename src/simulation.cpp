#include <mpi.h>
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

        // Create the nucleus
        auto n = new nucleus(charge, mass);
        for (int i=4; i<4+dimensions; ++i)
                n->coords[i] = std::stod(split[i]);

        template_system.push_back(n);

        // Add the specified number of electrons
        for (int i=0; i<electrons; ++i)
        {
                auto e = new electron();
                for (int i=0; i<dimensions; ++i)
                        e->coords[i] = n->coords[i];
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
                        target_population = std::stoi(split[1])/simulation.np;

                // Read in the number of DMC iterations
                else if (tag == "iterations")
                        dmc_iterations = std::stoi(split[1]);

                // Read in the DMC timestep
                else if (tag == "tau")
                        tau = std::stod(split[1]);

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
	// Initialize mpi
        if (MPI_Init(&argc, &argv) != 0) exit(MPI_ERROR);

        // Get the number of processes and my id within them
        if (MPI_Comm_size(MPI_COMM_WORLD, &np)  != 0) exit(MPI_ERROR);
        if (MPI_Comm_rank(MPI_COMM_WORLD, &pid) != 0) exit(MPI_ERROR);

        // Seed random number generator
        srand(simulation.pid*clock());

        // Read our input and setup parameters accordingly 
        // do for each process sequentially to avoid access issues
        for (int pid_read = 0; pid_read < simulation.np; ++ pid_read)
        {
                if (pid == pid_read) read_input();
                MPI_Barrier(MPI_COMM_WORLD);
        }

	// Open various output files
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
        //if (walker::count != 0 || particle::count != 0)
        //std::cout << "PID: " << pid << " un-deleted objects:\n"
        //         << "  Walkers   : " << walker::count   << "\n"
        //         << "  Particles : " << particle::count << "\n";

	MPI_Finalize();
}
