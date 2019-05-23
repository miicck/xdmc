#include <iostream>

#include "particle.h"
#include "potential.h"
#include "simulation.h"
#include "constants.h"
#include "random.h"
#include "walker.h"

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
	// Read input files, ready output files, initialize MPI etc.
	simulation.load(argc, argv);
	
	// For basic timing
	int start_clock = clock();

	// Our DMC walkers
	std::vector<walker*> walkers;

	// Initialize the walkers
	simulation.progress_file << "Initializing walkers...\n";
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
	simulation.progress_file << "Starting DMC simulation...\n";
	for (int iter = 1; iter <= simulation.dmc_iterations; ++iter)
	{
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

		// Output the main results of this iteration to track evolution
		simulation.evolution_file 
			<< walkers.size() << "," 
			<< simulation.trial_energy << "\n";
	}

	// Sample the final wavefunction to file
	for (int n=0; n<walkers.size(); ++n)
		walkers[n]->sample_wavefunction();

	// Clean up memory
	for (int i=0; i<walkers.size(); ++i) delete walkers[i];

	// Output success message
	double time = (clock() - start_clock)/double(CLOCKS_PER_SEC);
	simulation.progress_file << "Done, total time: " << time << "s.\n";

	// Free memory used in the simulation specification
	simulation.free_memory();
}
