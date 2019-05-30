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

#include "simulation.h"
#include "particle.h"
#include "walker.h"
#include "random.h"
#include "constants.h"

// Run the DMC calculation
void run_dmc()
{
	// Our DMC walkers
	simulation.progress_file << "Initializing walkers...\n";
	walker_collection walkers;
	
	// Run our DMC iterations
	simulation.progress_file << "Starting DMC simulation...\n";
	for (int iter = 1; iter <= simulation.dmc_iterations; ++iter)
	{
		simulation.progress_file 
			<< "Iteration: " << iter << "/" 
			<< simulation.dmc_iterations
			<< " Population: " << walkers.size()
			<< " Trial energy: " << simulation.trial_energy << "\n";

		// Carry out exchange moves on the walkers
		if (simulation.exchange_moves)
			walkers.make_exchange_moves();

		// Apply cancellation of walkers
		if (simulation.cancellation)
			walkers.apply_cancellations();

		// Apply diffusion-branching step
		walkers.diffuse_and_branch();

		// Set trial energy to control population, but allow some fluctuation
		double log_pop_ratio = log(walkers.size()/double(simulation.target_population));
		simulation.trial_energy = walkers.energy_estimate() - log_pop_ratio; 

		// Output the main results of this iteration to track evolution
		simulation.evolution_file 
			<< walkers.size()              << "," 
			<< simulation.trial_energy     << ","
			<< walkers.average_weight()    << ","
			<< walkers.average_weight_sq() << "\n";

		// Sample the wavefunction for the last 10 steps
		//if (simulation.dmc_iterations - iter < 10)
		walkers.sample_wavefunction();

		// Flush output files after every iteration
		simulation.flush();
	}

	// Output success message
	simulation.progress_file << "Done, total time: " << simulation.time() << "s.\n";
}

// Program entrypoint
int main(int argc, char** argv)
{
	// Read input files, ready output files, initialize MPI etc.
	simulation.load(argc, argv);

	// Run the DMC simulation
	run_dmc();

	// Free memory used in the simulation specification
	simulation.free_memory();
}
