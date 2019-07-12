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
#include "walker_collection.h"
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
		// Carry out exchange moves on the walkers
		if (simulation.exchange_moves)
			walkers.make_exchange_moves();

		// Apply cancellation of walkers
		if (simulation.make_cancellations)
			walkers.apply_cancellations();

		// Apply diffusion-branching step.
		walkers.diffuse_and_branch();

		// Output information about the walkers
		// at this iteration
		walkers.write_output(iter);
	}

	// Output success message
	simulation.progress_file << "\nDone, total time: " << simulation.time() << "s.\n";
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
