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

#include <sstream>
#include <mpi.h>

#include "random.h"
#include "math.h"
#include "walker_collection.h"
#include "simulation.h"

walker_collection :: walker_collection()
{
	// Reserve a reasonable amount of space to deal efficiently
	// with the fact that the population can fluctuate.
	walkers.reserve(simulation.target_population*4);

	// Initialize the set of walkers to the target
	// population size.
	for (int i=0; i<simulation.target_population; ++i)
	{
		walker* w = new walker();
		walkers.push_back(w);

		// Carry out initial diffusion to avoid
		// exact particle overlap on first iteration
		w->diffuse(simulation.pre_diffusion);
	}
}

walker_collection :: ~walker_collection()
{
	// Clean up memory
	for (int n=0; n<size(); ++n)
		delete (*this)[n];
}

double potential_greens_function(double pot_before, double pot_after)
{
	// Evaluate the potential-dependent part of the greens
	// function G_v(x,x',tau) = exp(-tau*(v(x)+v(x')-2E_T)/2)
        if (std::isinf(pot_after))  return 0;
        if (std::isinf(pot_before)) return 0;

        double av_potential = (pot_before + pot_after)/2;
        double exponent     = -simulation.tau*(av_potential - simulation.trial_energy);
	return exp(exponent);
}

int branch_from_weight(double weight)
{
	// Returns how many walkers should be produced
	// from one walker of the given weight
	int num = (int)floor(fabs(weight) + rand_uniform());
	return std::min(num, 3);
}

void walker_collection :: diffuse_and_branch()
{
	// Carry out diffusion and branching of the walkers
	// (this is done as one step because we need to know
	// the potential before and after diffusion to apply
	// the branching step)

	// Set the trial energy to control population
	double log_pop_ratio = log(sum_mod_weight()/double(simulation.target_population));
	simulation.trial_energy = average_potential() - log_pop_ratio;

	int nmax = size();
	for (int n=0; n < nmax; ++n)
	{
		walker* w = (*this)[n];

		// Diffuse according to the diffusive
		// part of the greens function
		// (recording the potential before/after)
		double pot_before = w->potential();
		w->diffuse(simulation.tau);
		double pot_after  = w->potential();

		// Multiply the weight by the potential-
		// dependent part of the greens function
		w->weight *= potential_greens_function(pot_before, pot_after);

		// Apply branching step, adding branched
		// survivors to the end of the collection
		int surviving = branch_from_weight(w->weight);
		for (int s=0; s<surviving; ++s)
			walkers.push_back(w->branch_copy());
	}

	// Delete the previous iterations walkers
	for (int n=0; n < nmax; ++n)
	{
		// Replace the nth walker with the
		// last walker, allowing us to shorten
		// the array from the end (delete the
		// nth walker in the process)
		delete walkers[n];
		walkers[n] = walkers.back();
		walkers.pop_back();
	}
}

double walker_collection :: sum_mod_weight()
{
	// Returns sum_i |w_i|
	// This is the effective population
	double sum = 0;
	for (int n=0; n<size(); ++n)
		sum += fabs((*this)[n]->weight);
	return sum;
}

double walker_collection :: average_mod_weight()
{
	// Returns (1/N) * sum_i |w_i|
	return sum_mod_weight()/double(size());
}

double walker_collection :: average_weight()
{
	// Returns (1/N) * sum_i w_i
	double av_weight = 0;
	for (int n=0; n<size(); ++n)
		av_weight += (*this)[n]->weight;
	av_weight /= double(size());
	return av_weight;
}

double walker_collection :: average_potential()
{
	// Returns (1/N) * sum_i |w_i|*v_i
	// Where v_i is the potential energy 
	// of the i^th walker.
	double energy = 0;
	for (int n=0; n<size(); ++n)
		energy += (*this)[n]->potential() * fabs((*this)[n]->weight);
	energy /= double(size());
	return energy;
}

void walker_collection :: make_exchange_moves()
{
	// Apply exchange moves to each of the walkers
	for (int n=0; n<size(); ++n)
		(*this)[n]->exchange();
}

void walker_collection :: apply_cancellations()
{
        // Record weights and cancellation probabilities
	double*  weights_before = new double [size()];
        for (int n=0; n<size(); ++n)
                weights_before[n] = (*this)[n]->weight;

        // Apply all pair-wise cancellations
	for (int n=0; n<size(); ++n)
        {
                walker* wn = (*this)[n];

                for (int m=0; m<n; ++m)
                {
                        walker* wm = (*this)[m];
                        double cp  = wn->cancel_prob(wm);
                        
                        // Cancel these walkers with probability cp
                        wn->weight *= 1.0 - cp;
                        wm->weight *= 1.0 - cp;
                }
        }

	// Evaluate a measure of the amount of cancellation
	// that occured
	cancellation_amount_last = 0;
	for (int n=0; n<size(); ++n)
	{
		double delta = (*this)[n]->weight - weights_before[n];
		cancellation_amount_last += fabs(delta);
	}

        // Clean up memory
        delete[] weights_before;
}

void walker_collection :: write_output(int iter)
{
        // Sum up walkers across processes
        int population = this->size();
        int population_red;
        MPI_Reduce(&population, &population_red, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        // Average trial energy across processes
        double triale = simulation.trial_energy;
        double triale_red;
        MPI_Reduce(&triale, &triale_red, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        triale_red /= double(simulation.np);

        // Avarage <weight> across processes
        double av_weight = this->average_weight();
        double av_weight_red;
        MPI_Reduce(&av_weight, &av_weight_red, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        av_weight_red /= double(simulation.np);

        // Average <|weight|> across processes
        double av_mod_weight = this->average_mod_weight();
        double av_mod_weight_red;
        MPI_Reduce(&av_mod_weight, &av_mod_weight_red, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        av_mod_weight_red /= double(simulation.np);

	// Sum cancellation amount across processes
	double cancel = this->cancellation_amount_last;
	double cancel_red;
	MPI_Reduce(&cancel, &cancel_red, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Calculate timing information
	double time_per_iter     = simulation.time()/iter;
	double seconds_remaining = time_per_iter * (simulation.dmc_iterations - iter);
	double percent_complete  = double(100*iter)/simulation.dmc_iterations;

        // Output iteration information
        simulation.progress_file << "\nIteration " << iter << "/" << simulation.dmc_iterations
				 << " (" << percent_complete << "%)\n";
	simulation.progress_file << "    Time running    : " << simulation.time()
				 << "s (" << time_per_iter << "s/iter)\n";
	simulation.progress_file << "    ETA             : "
				 << seconds_remaining << "s \n";
        simulation.progress_file << "    Trial energy    : " << triale_red     << " Hartree\n";
        simulation.progress_file << "    Population      : " << population_red << "\n";
	simulation.progress_file << "    Canceled weight : " << cancel_red     << "\n";

        if (iter == 1)
        {
                // Before the first iteration, output names of the
                // evolution file columns
                simulation.evolution_file
                        << "Population,"
                        << "Trial energy (Hartree),"
                        << "Average weight,"
                        << "Average |weight|,"
			<< "Cancelled weight\n";
        }

        // Output evolution information to file
        simulation.evolution_file
                        << population_red    << ","
                        << triale_red        << ","
                        << av_weight_red     << ","
                        << av_mod_weight_red << ","
			<< cancel_red        << "\n";

        // Write the wavefunction to file
        if (simulation.write_wavefunction)
        {
                simulation.wavefunction_file << "# Iteration " << iter << "\n";
		for (int n=0; n<size(); ++n)
		{
			(*this)[n]->write_wavefunction();
			simulation.wavefunction_file << "\n";
		}
        }

        // Flush output files after every call
        simulation.flush();
}
