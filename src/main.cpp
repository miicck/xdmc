#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iterator>
#include <mpi.h>
#include <unistd.h>

#include "particle.h"
#include "potential.h"
#include "simulation.h"
#include "constants.h"
#include "random.h"
#include "walker.h"

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
void parse_atom(std::vector<std::string> split)
{
	double charge = std::stod(split[1]);
	int electrons = std::stoi(split[2]);
	double mass   = std::stod(split[3]);

	// Create the nucleus
	auto n = new nucleus(charge, mass);
	for (int i=4; i<4+simulation.dimensions; ++i)
		n->coords[i] = std::stod(split[i]);

	simulation.template_system.push_back(n);

	// Add the specified number of electrons
	for (int i=0; i<electrons; ++i)
	{
		auto e = new electron();
		for (int i=0; i<simulation.dimensions; ++i)
			e->coords[i] = n->coords[i];
		simulation.template_system.push_back(e);
	}
}

// Parse the input file.
void read_input()
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
			simulation.dimensions = std::stoi(split[1]);

		// Read in the number of DMC walkers and convert
		// to walkers-per-process
		else if (tag == "walkers") 
			simulation.target_population = std::stoi(split[1])/simulation.np;

		// Read in the number of DMC iterations
		else if (tag == "iterations")
			simulation.dmc_iterations = std::stoi(split[1]);

		// Read in the DMC timestep
		else if (tag == "tau")
			simulation.tau = std::stod(split[1]);

		// Read in an atom
		else if (tag == "atom")
			parse_atom(split);

		// Read in an electron
		else if (tag == "electron")
			simulation.template_system.push_back(new electron());

		// Adds a noninteracting fermion into the system
		else if (tag == "nif")
			simulation.template_system.push_back(new non_interacting_fermion());

		// Add a harmonic well to the system
		else if (tag == "harmonic_well")
			simulation.potentials.push_back(new harmonic_well(std::stod(split[1])));
	}

	std::cout <<   "  PID: " << simulation.pid 
		  << "\n    Dimensionallity   : " << simulation.dimensions
		  << "\n    Iterations        : " << simulation.dmc_iterations
	          << "\n    Target population : " << simulation.target_population
		  << "\n    Particles         : " << simulation.template_system.size()
		  << "\n    Potentials        : " << simulation.potentials.size() << "\n";

	input.close();
}

// Returns the number of walkers that should survive after a 
// walker moves from potential pot_before to pot_after. i.e returns
// 0 if the walker should die and 1 if it should live, 2 if another should
// spawn etc...
int walkers_surviving(double pot_before, double pot_after)
{
	double p = exp(-simulation.tau*(pot_before + pot_after - 2*simulation.trial_energy)/2);
	return std::min(int(floor(p+rand_uniform())), 3);
}

// The main DMC algorithm
int main(int argc, char** argv)
{
	// Initialize mpi
	if (MPI_Init(&argc, &argv) != 0) exit(MPI_ERROR);

	// For basic timing
	int start_clock = clock();

	// Get the number of processes and my id within them
	if (MPI_Comm_size(MPI_COMM_WORLD, &simulation.np)  != 0) exit(MPI_ERROR);
	if (MPI_Comm_rank(MPI_COMM_WORLD, &simulation.pid) != 0) exit(MPI_ERROR);

	if (simulation.pid == 0) 
		std::cout << simulation.np << " processes reporting for duty:\n";

	// Seed random number generator
	srand(simulation.pid*clock());

	// Read our input and setup parameters accordingly 
	// do for each process sequentially to avoid access issues
	for (int pid_read = 0; pid_read < simulation.np; ++ pid_read)
	{
		if (simulation.pid == pid_read) read_input();
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	// Ready output files
	simulation.open_output_files();

	// Our DMC walkers
	std::vector<walker*> walkers;

	// Initialize the walkers
	for (int i=0; i<simulation.target_population; ++i)
	{
		walker* w = new walker();
		walkers.push_back(w);

		// Carry out an initial large diffusion
		// to avoid unneccasary equilibriation
		// and exact particle overlap at the origin
		w->diffuse(1.0);
	}
	
	// Run our DMC iterations
	for (int iter = 1; iter <= simulation.dmc_iterations; ++iter)
	{
		// Output progress on root node
		if (simulation.pid == 0) std::cout << "\rIteration: " << iter << "/"
					<< simulation.dmc_iterations << "     ";

		// The new walkers that appear after the iteration
		std::vector<walker*> new_walkers;

		// Variables to accumulate
		double average_potential = 0;

		for (int n=0; n<walkers.size(); ++n)
		{
			walker* w = walkers[n];

			// Diffuse the walker and evaluate potential change
			double pot_before = w->potential();
			w->diffuse(simulation.tau);
			double pot_after  = w->potential();

			// Work out how many survivded the move
			int surviving = walkers_surviving(pot_before, pot_after);

			if (surviving == 0) delete(w); // Delete this walker if it died
			else new_walkers.push_back(w); // Otherwise it survives

			// Spawn new walkers
			for (int s=0; s < surviving-1; ++s) 
				new_walkers.push_back(w->copy());

			// Accumulate the average potential
			average_potential += surviving*pot_after;
		}

		// Set our walkers to the newly generated ones
		walkers = new_walkers;

		// Set trial energy to control population, but allow some fluctuation
		average_potential /= double(walkers.size());
		double log_pop_ratio = log(walkers.size()/double(simulation.target_population));
		simulation.trial_energy = average_potential - log_pop_ratio; 

		// Output the results of this iteration
		simulation.iterations_file 
			<< walkers.size() << "," 
			<< simulation.trial_energy << "\n";
	}
	if (simulation.pid == 0) std::cout << "\n";

	// Sample the final wavefunction to file
	for (int n=0; n<walkers.size(); ++n)
		walkers[n]->sample_wavefunction();

	// Clean up memory
	for (int i=0; i<walkers.size(); ++i) delete walkers[i];
	simulation.free_memory();

	// Output info on objects that werent deconstructed properly
	if (walker::count != 0 || particle::count != 0)
	std::cout << "PID: " << simulation.pid << " un-deleted objects:\n"
		  << "  Walkers   : " << walker::count   << "\n"
		  << "  Particles : " << particle::count << "\n";

	// Output success message
	if (simulation.pid == 0)
	{
		double time = (clock() - start_clock)/double(CLOCKS_PER_SEC);
		std::cout << "Done, total time: " << time << "s.\n";
	}

	// End mpi safely
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
