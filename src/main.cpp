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

#include <vector>
#include <iostream>

#include "simulation.h"
#include "particle.h"
#include "walker.h"
#include "random.h"
#include "constants.h"

// Returns the number of walkers that should survive after a 
// walker moves from potential pot_before to pot_after. i.e returns
// 0 if the walker should die and 1 if it should live, 2 if another should
// spawn etc...
int walkers_surviving(double pot_before, double pot_after)
{
	if (std::isinf(pot_after)) return 0;
	if (std::isinf(pot_before)) return 0;

	double p = exp(-simulation.tau*(pot_before + pot_after - 2*simulation.trial_energy)/2);
	return std::min(int(floor(p+rand_uniform())), 3);
}

// Renormalize the set of walkers such that sum |w_i|^2 = 1
void renormalize(std::vector<walker*> &walkers)
{
	double av_weight_sq = 0;
	for (int n=0; n<walkers.size(); ++n)
		av_weight_sq += walkers[n]->weight * walkers[n]->weight;

	av_weight_sq /= double(walkers.size());
	for (int n=0; n<walkers.size(); ++n)
		walkers[n]->weight /= sqrt(av_weight_sq);
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

		// Carry out an initial diffusion
		// to avoid exact particle overlap
		w->diffuse();
	}
	
	// Run our DMC iterations
	simulation.progress_file << "Starting DMC simulation...\n";
	for (int iter = 1; iter <= simulation.dmc_iterations; ++iter)
	{
		// The new walkers that appear after the iteration
		std::vector<walker*> new_walkers;

		// Variables to accumulate
		double average_potential = 0;
		double average_weight = 0;
		double average_weight_sq = 0;

		// Carry out exchange moves on the walkers
		for (int n=0; n<walkers.size(); ++n)
			walkers[n]->exchange();

		// Apply cancellation of walkers
		for (int n=0; n<walkers.size(); ++n)
			for (int m=0; m<n; ++m)
				walkers[n]->cancel(walkers[m]);

		// Renormalize weights
		renormalize(walkers);

		for (int n=0; n<walkers.size(); ++n)
		{
			walker* w = walkers[n];

			// Diffuse the walker and evaluate potential change
			double pot_before = w->potential();
			w->diffuse();
			double pot_after  = w->potential();

			// Work out how many will survive the move
			int surviving = walkers_surviving(pot_before, pot_after);
			if (w->weight*w->weight < 0.1)
				surviving = 0;

			// Accumulate averages
			average_potential  += surviving*pot_after;
			average_weight     += surviving*w->weight;
			average_weight_sq  += surviving*w->weight*w->weight; 

			if (surviving == 0) delete(w); // Delete this walker if it died
			else new_walkers.push_back(w); // Otherwise it survives

			// Spawn new walkers
			for (int s=0; s < surviving-1; ++s) 
				new_walkers.push_back(w->copy());
		}

		// Set our walkers to the newly generated ones
		walkers = new_walkers;

		average_potential /= double(walkers.size());
		average_weight    /= double(walkers.size());
		average_weight_sq /= double(walkers.size());

		// Set trial energy to control population, but allow some fluctuation
		double log_pop_ratio = log(walkers.size()/double(simulation.target_population));
		simulation.trial_energy = average_potential - log_pop_ratio; 

		// Output the main results of this iteration to track evolution
		simulation.evolution_file 
			<< walkers.size()          << "," 
			<< simulation.trial_energy << ","
			<< average_weight          << ","
			<< average_weight_sq       << "\n";
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




















