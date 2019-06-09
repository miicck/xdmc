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
#include "simulation.h"
#include "particle.h"
#include "walker.h"
#include "random.h"
#include "constants.h"

// Forward decleration
void mpi_reduce_iteration(walker_collection& walkers, int iter);

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
		// Apply diffusion-branching step.
		walkers.diffuse_and_branch();

		// Carry out exchange moves on the walkers
		walkers.make_exchange_moves();

		// Apply cancellation of walkers
		walkers.apply_cancellations();

		// Reduce this iteration and output
		mpi_reduce_iteration(walkers, iter);
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

void mpi_reduce_iteration(walker_collection& walkers, int iter)
{
	// Sum up walkers across processes
	int population = walkers.size();
	int population_red;
	MPI_Reduce(&population, &population_red, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	// Average trial energy across processes
	double triale = simulation.trial_energy;
	double triale_red;
	MPI_Reduce(&triale, &triale_red, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	triale_red /= double(simulation.np);

	// Avarage <weight> across processes
	double av_weight = walkers.average_weight();
	double av_weight_red;
	MPI_Reduce(&av_weight, &av_weight_red, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	av_weight_red /= double(simulation.np);

	// Average <|weight|> across processes
	double av_mod_weight = walkers.average_mod_weight();
	double av_mod_weight_red;
	MPI_Reduce(&av_mod_weight, &av_mod_weight_red, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	av_mod_weight_red /= double(simulation.np);

	// Output iteration information
	simulation.progress_file << "\nIteration " << iter << "/" << simulation.dmc_iterations << "\n";
	simulation.progress_file << "    Trial energy   : " << triale_red     << " Hartree\n";
	simulation.progress_file << "    Population     : " << population_red << "\n";

	if (iter == 1)
	{
		// Before the first iteration, output names of the
		// evolution file columns
		simulation.evolution_file
			<< "Population,"
			<< "Trial energy (Hartree),"
			<< "Average weight,"
			<< "Average |weight|\n";
	}

	// Output evolution information to file
	simulation.evolution_file
			<< population_red    << ","
			<< triale_red        << ","
			<< av_weight_red     << ","
			<< av_mod_weight_red << "\n";

	// Write the wavefunction to file
	if (simulation.write_wavefunction)
	{
		simulation.wavefunction_file << "# Iteration " << iter << "\n";
		walkers.write_wavefunction();
	}

	// Flush output files after every call
	simulation.flush();
}
