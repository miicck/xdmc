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
#include "params.h"

walker_collection :: walker_collection()
{
    // Reserve a reasonable amount of space to deal efficiently
    // with the fact that the population can fluctuate.
    walkers.reserve(2*int(params::target_population));

    // Initialize the set of walkers to the target
    // population size.
    for (unsigned i=0; i<params::target_population; ++i)
    {
        walker* w = new walker();
        w->diffuse(params::pre_diffusion);
        w->reflect_to_irreducible();
        walkers.push_back(w);
    }
}

walker_collection :: ~walker_collection()
{
    // Clean up memory
    for (unsigned n=0; n<walkers.size(); ++n)
        delete walkers[n];
}

walker_collection* walker_collection :: copy()
{
    // Create an exact copy of this collection
    // of walkers
    std::vector<walker*> copied_walkers;
    for (unsigned n=0; n<walkers.size(); ++n)
        copied_walkers.push_back(walkers[n]->copy());
    return new walker_collection(copied_walkers);
}

void walker_collection :: propagate()
{
    // Apply the stages of walker propagation
    make_exchange_moves();
    make_diffusive_moves();
    apply_renormalization();
    branch();
}

void walker_collection :: make_exchange_moves()
{
    // Apply exchange moves to each of the walkers
    for (unsigned n=0; n<walkers.size(); ++n)
        walkers[n]->exchange();
}

double walker_collection :: diffused_wavefunction(walker* c)
{
    // Evaluate the diffused wavefunction 
    // at the configuration of c: 
    // \psi_D(c) = \sum_i w_i G_D(c, x_i, dt)

    double psi_d = 0;
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];
        psi_d += w->weight * w->diffusive_greens_function(c);
    }
    return psi_d;
}

double potential_greens_function(double pot_before, double pot_after)
{
    // Evaluate the potential-dependent part of the greens
    // function G_v(x,x',tau) = exp(-tau*(v(x)+v(x'))/2)
    if (std::isinf(pot_after))  return 0;
    if (std::isinf(pot_before)) return 0;
    return fexp( -params::tau * (pot_before + pot_after)/2.0 );
}

void walker_collection :: make_diffusive_moves()
{
    // Record locations before diffusion
    walker_collection* walkers_last = this->copy();

    // Carry out diffusion of the walkers
    int accepted_moves = 0;
    int rejected_moves = 0;
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w_after  = walkers[n];
        walker* w_before = w_after->copy();

        // Make some proposed move
        // (note this is a metropolis proposal step, so
        //  does not have to be a greens function; any
        //  ergodic step will do - this one just ensures
        //  a reasonable acceptance probability).
        w_after->diffuse(params::tau*10);

        // Record the diffusive wavefunction and the potential
        // before and after the proposed move
        double psi_before = walkers_last->diffused_wavefunction(w_before);
        double psi_after  = walkers_last->diffused_wavefunction(w_after);

        double pot_before = w_before->potential();
        double pot_after  = w_after->potential();

        // Accept the move with probability given by
        // the diffused wavefunction
        double accept_prob = fabs(psi_after/psi_before);

        if (rand_uniform() < accept_prob)
        {
            // Accept new walker
            // inherit sign of diffused wavefunction
            // apply potential part of greens function
            walkers[n]       = w_after;
            w_after->weight  = sign(psi_after);
            w_after->weight *= potential_greens_function(pot_after, pot_before);
            accepted_moves  += 1;
            delete w_before;
        }
        else
        {
            // Revert to previous walker
            // inherit sign of diffused wavefunction
            // apply potential part of greens function
            walkers[n]        = w_before;
            w_before->weight  = sign(psi_before);
            w_before->weight *= potential_greens_function(pot_before, pot_before);
            rejected_moves   += 1;
            delete w_after;
        }

    }

    // Calculate accpetance ratio
    params::acceptance_ratio = double(accepted_moves)/
                               double(accepted_moves + rejected_moves);

    // Free memory
    delete walkers_last;
}

void walker_collection :: apply_renormalization()
{
    // Set the trial energy and apply renormalization
    // so that the population fluctuates around the 
    // target population

    // The population at the start of the iteration
    double pop_before_propagation = double(walkers.size());

    // The effective population now, after the cumulative
    // effect of this iterations greens functions
    // (i.e cancellation, diffusion, potential etc...)
    double pop_after_propagation  = sum_mod_weight();

    // Set trial energy to minimize fluctuations
    params::trial_energy  = log(pop_before_propagation / pop_after_propagation)/params::tau;

    // Bias towards target population
    params::trial_energy -= log(pop_before_propagation / params::target_population);

    // Apply normalization greens function
    double gn = fexp(params::trial_energy * params::tau);
    for (unsigned n=0; n<walkers.size(); ++n)
        walkers[n]->weight *= gn;
}

int branch_from_weight(double weight)
{
    // Returns how many walkers should be produced
    // from one walker of the given weight
    return (int)floor(fabs(weight) + rand_uniform());
}

void walker_collection :: branch()
{
    // Carry out branching of the walkers
    unsigned nmax = walkers.size();
    for (unsigned n=0; n < nmax; ++n)
    {
        walker* w = walkers[n];

        // Apply branching step, adding branched
        // survivors to the end of the collection
        int surviving = branch_from_weight(w->weight);
        for (int s=0; s<surviving; ++s)
            walkers.push_back(w->branch_copy());
    }

    // Delete the previous iterations walkers
    for (unsigned n=0; n < nmax; ++n)
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
    for (unsigned n=0; n<walkers.size(); ++n)
        sum += fabs(walkers[n]->weight);
    return sum;
}

double walker_collection :: average_mod_weight()
{
    // Returns (1/N) * sum_i |w_i|
    return sum_mod_weight()/double(walkers.size());
}

double walker_collection :: average_weight()
{
    // Returns (1/N) * sum_i w_i
    double av_weight = 0;
    for (unsigned n=0; n<walkers.size(); ++n)
        av_weight += walkers[n]->weight;
    av_weight /= double(walkers.size());
    return av_weight;
}

double walker_collection :: average_potential()
{
    // Returns (1/W) * sum_i |w_i|*v_i.
    // Where v_i is the potential energy 
    // of the i^th walker (which has
    // weight w_i) and W = sum_i |w_i|
    double pot = 0;
    double weight = 0;
    for (unsigned n=0; n<walkers.size(); ++n)
    {
        pot    += walkers[n]->potential() * fabs(walkers[n]->weight);
        weight += fabs(walkers[n]->weight);
    }
    pot /= weight;
    return pot;
}

double mpi_average(double val)
{
    // Get the average of val across proccesses on pid 0
    double res;
    MPI_Reduce(&val, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    res /= double(params::np);
    return res;
}

double mpi_sum(double val)
{
    // Get the sum of val across proccesses on pid 0
    double res;
    MPI_Reduce(&val, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return res;
}

void walker_collection :: write_output()
{
    // Sum various things across processes
    double population_red    = mpi_sum(double(walkers.size()));

    // Average various things across processes
    double triale_red        = mpi_average(params::trial_energy);
    double av_weight_red     = mpi_average(average_weight());
    double potential_red     = mpi_average(average_potential());
    double acceptance_red    = mpi_average(params::acceptance_ratio);

    // Calculate timing information
    double time_per_iter     = params::dmc_time()/params::dmc_iteration;
    double secs_remain       = time_per_iter * (params::dmc_iterations - params::dmc_iteration);
    double percent_complete  = double(100*params::dmc_iteration)/params::dmc_iterations;

    // Output iteration information
    params::progress_file << "\nIteration " << params::dmc_iteration 
                          << "/"  << params::dmc_iterations
                          << " (" << percent_complete << "%)\n";
    params::progress_file << "    Time running       : " << params::time()
                          << "s (" << time_per_iter << "s/iter)\n";
    params::progress_file << "    ETA                : " << secs_remain    << "s \n";
    params::progress_file << "    Trial energy       : " << triale_red     << " Hartree\n";
    params::progress_file << "    <V>                : " << potential_red  << " Hartree\n";
    params::progress_file << "    Population         : " << population_red
                          << " (" << population_red/params::np << " per process) "<< "\n";
    params::progress_file << "    Acceptance ratio   : " << acceptance_red*100 << "%\n";

    if (params::dmc_iteration == 1)
    {
        // Before the first iteration, output names of the
        // evolution file columns
        params::evolution_file
                << "Population,"
                << "Trial energy,"
                << "<V>,"
                << "Acceptance ratio,"
                << "Average weight\n";
    }

    // Output evolution information to file
    params::evolution_file
        << population_red    << ","
        << triale_red        << ","
        << potential_red     << ","
        << acceptance_red    << ","
        << av_weight_red     << "\n";

    // Write the wavefunction to file
    if (params::write_wavefunction)
    {
        params::wavefunction_file << "# Iteration " << params::dmc_iteration << "\n";
        for (unsigned n=0; n<walkers.size(); ++n)
        {
            walkers[n]->write_wavefunction();
            params::wavefunction_file << "\n";
        }
    }

    // Flush output files after every call
    params::flush();
}
