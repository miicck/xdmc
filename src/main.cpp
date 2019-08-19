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

#include "params.h"
#include "particle.h"
#include "walker_collection.h"
#include "random.h"
#include "constants.h"

// Run the DMC calculation
void run_dmc()
{
	// Our DMC walkers
	params::progress_file << "Initializing walkers...\n";
	walker_collection walkers;
	
	// Run our DMC iterations
	params::progress_file << "Starting DMC simulation...\n";
	for (int iter = 1; iter <= params::dmc_iterations; ++iter)
	{
        // Reset expectation values
        walkers.expect_vals.reset();

		// Apply diffusion-branching step.
		walkers.diffuse_and_branch();

		// Carry out exchange moves on the walkers
		if (params::exchange_moves)
			walkers.make_exchange_moves();

		// Apply cancellation of walkers
        walkers.apply_cancellations();

        // Apply walker seperation-correction
        if (params::correct_seperations)
                walkers.correct_seperations();

        // Normalize expectation values
        walkers.expect_vals.normalize(walkers.size());

        // Set trial energy to control population
        double log_pop_ratio = log(double(walkers.size())/double(params::target_population));
        params::trial_energy = walkers.expect_vals.average_potential - log_pop_ratio;

		// Output information about the walkers
		// at this iteration
		walkers.write_output(iter);
	}

	// Output success message
	params::progress_file << "\nDone, total time: " << params::time() << "s.\n";
}

// Program entrypoint
int main(int argc, char** argv)
{
	// Read input files, ready output files, initialize MPI etc.
	params::load(argc, argv);

	// Run the DMC simulation
	run_dmc();

	// Free memory used in the simulation specification
	params::free_memory();
}
