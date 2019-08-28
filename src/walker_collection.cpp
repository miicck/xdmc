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
    walkers.reserve(int(params::target_population*
                        params::max_pop_ratio)+10);

    // Initialize the set of walkers to the target
    // population size.
    walker* w = new walker();
    double diff_amt = params::pre_diffusion/2.0;
    w->diffuse(diff_amt);

    for (unsigned i=0; i<params::target_population; ++i)
    {
        // Initialize walkers with a large
        // antisymmetric componenet if we have
        // fermionic exchanges.
        w->exchange();
        walker* wcopy = w->branch_copy();
        wcopy->diffuse(diff_amt);
        walkers.push_back(wcopy);
    }

    delete w;
}

walker_collection :: ~walker_collection()
{
    // Clean up memory
    for (unsigned n=0; n<size(); ++n)
        delete (*this)[n];
}

walker_collection* walker_collection :: copy()
{
    // Create an exact copy of this collection
    // of walkers
    std::vector<walker*> copied_walkers;
    for (unsigned n=0; n<size(); ++n)
        copied_walkers.push_back((*this)[n]->copy());
    return new walker_collection(copied_walkers);
}

double potential_greens_function(double pot_before, double pot_after)
{
    // Evaluate the potential-dependent part of the greens
    // function G_v(x,x',tau) = exp(-tau*(v(x)+v(x')-2E_T)/2)
    if (std::isinf(pot_after))  return 0;
    if (std::isinf(pot_before)) return 0;
    return fexp( -params::tau * (pot_before + pot_after)/2.0 );
}

void walker_collection :: diffuse()
{
    // Carry out diffusion and branching of the walkers
    for (unsigned n=0; n < size(); ++n)
    {
        walker* w = (*this)[n];

        // Diffuse according to the diffusive
        // part of the greens function
        // (recording the potential before/after)
        double pot_before = w->potential();
        w->diffuse(params::tau);
        double pot_after  = w->potential();

        // Multiply the weight by the potential-
        // dependent part of the greens function
        w->weight *= potential_greens_function(pot_before, pot_after);
    }
}

int branch_from_weight(double weight)
{
    // Returns how many walkers should be produced
    // from one walker of the given weight
    int num = (int)floor(fabs(weight) + rand_uniform());
    if (num > 3)
    {
        // Warn user about large weights, often a sign of population explosion
        params::error_file << "Warning: Large weight detected; " << weight 
                           << " which will branch into "      << num 
                           << " walkers.\n";
    }
    return num;
}

void walker_collection :: clip_weight()
{
    // Adjust weights if expected population after
    // branching is unacceptable
    double amw        = this->average_mod_weight();
    double ex_weight  = amw*size();
    double max_weight = params::target_population * params::max_pop_ratio;
    double min_weight = params::target_population * params::min_pop_ratio;
    if (ex_weight < min_weight || ex_weight > max_weight)
    {
        // Warn the user this has happened
        params::error_file << "Warning: clipping applied at iteration "
                           << params::dmc_iteration 
                           << " this will introduce bias!\n";

        for (unsigned n=0; n<size(); ++n)
            (*this)[n]->weight /= amw;
    }
}

void walker_collection :: apply_renormalization()
{
    // Set the trial energy and apply renormalization
    // so that the population fluctuates around the 
    // target population

    // The population at the start of the iteration
    double pop_before_propagation = double(size());

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
    for (unsigned n=0; n<size(); ++n)
        (*this)[n]->weight *= gn;

    // Clip weight so that the population
    // remains between the min and max allowed values
    clip_weight();
}

void walker_collection :: branch()
{
    // Update trial energy/renormalize weights
    this->apply_renormalization();

    // Carry out branching of the walkers
    unsigned nmax = size();
    for (unsigned n=0; n < nmax; ++n)
    {
        walker* w = (*this)[n];

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
    for (unsigned n=0; n<size(); ++n)
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
    for (unsigned n=0; n<size(); ++n)
        av_weight += (*this)[n]->weight;
    av_weight /= double(size());
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
    for (unsigned n=0; n<size(); ++n)
    {
        pot    += (*this)[n]->potential() * fabs((*this)[n]->weight);
        weight += fabs((*this)[n]->weight);
    }
    pot /= weight;
    return pot;
}

double walker_collection :: average_kinetic()
{
    // Returns (1/W) * sum_i |w_i|*k_i.
    // Where k_i is the kinetic energy 
    // of the i^th walker (which has
    // weight w_i) and W = sum_i |w_i|
    double kin = 0;
    double weight = 0;
    for (unsigned n=0; n<size(); ++n)
    {
        kin    += (*this)[n]->kinetic() * fabs((*this)[n]->weight);
        weight += fabs((*this)[n]->weight);
    }
    kin /= weight;
    return kin;
}

void walker_collection :: make_exchange_moves()
{
    // Don't make exchange moves if they're turned off
    if (!params::exchange_moves)
        return;

    // Apply exchange moves to each of the walkers
    for (unsigned n=0; n<size(); ++n)
        (*this)[n]->exchange();
}

void walker_collection :: apply_cancellations(walker_collection* walkers_last)
{
    double weight_before = sum_mod_weight();

    // Apply the selected cancellation scheme.
    // Using walkers_last to construct trial wavefunctions
    // where a trial wavefunction is required.
    if      (params::cancel_scheme == "voronoi")
        this->apply_voronoi_cancellations();
    else if (params::cancel_scheme == "pairwise")
        this->apply_pairwise_cancellations();
    else if (params::cancel_scheme == "diffusive")
        this->apply_diffusive_cancellations(walkers_last);
    else if (params::cancel_scheme == "none")
        return;
    else
        params::error_file << "Unknown cancellation scheme: "
                              << params::cancel_scheme << "\n";

    // Record the amount of cancelled weight
    params::cancelled_weight = weight_before - sum_mod_weight();
}

void walker_collection :: apply_diffusive_cancellations(walker_collection* walkers_last)
{
    double* new_weights = new double[size()];

    for (unsigned n=0; n<size(); ++n)
    {
        walker* wn = (*this)[n];

        // Evaluate total diffusive Greens function
        // at the new configuration of wn:
        //   gf = \sum_{m!=n} w_m G_D(x_n', x_m, dt)
        //      + sss w_n G_D(x_n', x_n, dt)
        // where sss is the self sign strength.
        double gf = 0;
        for (unsigned m=0; m<walkers_last->size(); ++m)
        {
            double factor = 1.0;
            if (m == n) factor = params::self_sign_strength;
            walker* wm = (*walkers_last)[m];
            gf += factor * wm->diffusive_greens_function(wn) * wm->weight;
        }

        // Kill the walker if it ended up
        // in an opposite-sign region
        if (sign(wn->weight) == sign(gf))
            new_weights[n] = wn->weight;
        else
            new_weights[n] = 0.0;
    }

    // Apply new weights
    for (unsigned n=0; n<size(); ++n)
        (*this)[n]->weight = new_weights[n];

    // Free memory
    delete[] new_weights;
}

void walker_collection :: apply_voronoi_cancellations()
{
    // Kill walkers who have more nearest-neighbours of
    // the opposite sign than their own sign
    // Uses enough nearest neighbours to construct an
    // n-simplex in configuration space.

    // Allocate space for new weights
    double* new_weights = new double[size()];

    // Work out number of nearest neighbours
    unsigned nn = params::dimensions*params::template_system.size() + 1;
    if (size() < 1 + nn)
    {
        params::error_file << "Error in voronoi: population very small!\n";
        return;
    }

    for (unsigned n=0; n<size(); ++n)
    {
        walker* wn = (*this)[n];

        // Find the nn nearest neighbours
        double*   sq_distances = new double[nn];
        unsigned* nn_indicies  = new unsigned[nn];

        for (unsigned i=0; i<nn; ++i)
        {
            nn_indicies[i]  = -1;
            sq_distances[i] = INFINITY;
        }

        for (unsigned m=0; m<size(); ++m)
        {
            if (m == n) continue;

            walker* wm = (*this)[m];
            double sq_dis = wn->sq_distance_to(wm);
            
            for (unsigned i=0; i<nn; ++i)
                if (sq_dis < sq_distances[i])
                {
                    // Shift records along by one to make
                    // space for new entry
                    for (unsigned j=nn-1; j>i; --j)
                    {
                        nn_indicies[j]  = nn_indicies[j-1];
                        sq_distances[j] = sq_distances[j-1];
                    }

                    // Fill space with new entry
                    nn_indicies[i]  = m;
                    sq_distances[i] = sq_dis;
                    break;
                }
        }

        // Record amount of +ve and -ve
        // weight nearby
        double positive_weight = 0;
        double negative_weight = 0;

        for (unsigned i=0; i<nn; ++i)
        {
            if (nn_indicies[i] < 0) 
            {
                params::error_file << "Neg nn ind @ " << n << ", " << i << "\n";
                continue;
            }

            double nnw = (*this)[nn_indicies[i]]->weight;
            if (nnw < 0) negative_weight -= nnw;
            else positive_weight += nnw;
        }

        // Kill the walker if his team is outnumbered
        new_weights[n] = wn->weight;
        if (wn->weight < 0)
        {
            // Walker is -ve, kill it if +ve is winning
            if (positive_weight > negative_weight)
                new_weights[n] = 0;
        }
        else
        {
            // Walker is +ve, kill it if -ve is winning
            if (negative_weight > positive_weight)
                new_weights[n] = 0;
        }

        // Free memory
        delete[] sq_distances;
        delete[] nn_indicies;
    }

    // Don't apply cancellations if a large fraction
    // of the walkers will die
    double total_weight = 0;
    for (unsigned n=0; n<size(); ++n)
        total_weight += fabs(new_weights[n]);
    double frac_lost = 1 - total_weight/double(size());
    if (frac_lost > 0.5) return;

    // Apply new weights
    for (unsigned n=0; n<size(); ++n)
        (*this)[n]->weight = new_weights[n];

    // Free memory
    delete[] new_weights;
}

void walker_collection :: apply_pairwise_cancellations()
{
    // Apply all pair-wise cancellations
    for (unsigned n=0; n<size(); ++n)
    {
        walker* wn = (*this)[n];

        for (unsigned m=0; m<n; ++m)
        {
            walker* wm = (*this)[m];
            double cp  = wn->cancel_prob(wm);
            
            // Cancel these walkers with probability cp
            wn->weight *= 1.0 - cp;
            wm->weight *= 1.0 - cp;
        }
    }
}

void walker_collection :: correct_seperations()
{
    // Don't do anything if seperation correction is off
    if (!params::correct_seperations)
        return;
    
    // Correct walker seperations in a pairwise manner
    for (unsigned n=0; n<size(); ++n)
        for (unsigned m=0; m<n; ++m)
            (*this)[n]->drift_away_from((*this)[m]);
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
    double population_red    = mpi_sum(double(this->size()));
    double cancelled_weight  = mpi_sum(params::cancelled_weight);

    // Average various things across processes
    double triale_red        = mpi_average(params::trial_energy);
    double av_weight_red     = mpi_average(this->average_weight());
    double potential_red     = mpi_average(this->average_potential());

    // Calculate timing information
    double time_per_iter     = params::time()/params::dmc_iteration;
    double secs_remain       = time_per_iter * (params::dmc_iterations - params::dmc_iteration);
    double percent_complete  = double(100*params::dmc_iteration)/params::dmc_iterations;

    // Output iteration information
    params::progress_file << "\nIteration " << params::dmc_iteration 
                          << "/" << params::dmc_iterations
                          << " (" << percent_complete << "%)\n";
    params::progress_file << "    Time running       : " << params::time()
                          << "s (" << time_per_iter << "s/iter)\n";
    params::progress_file << "    ETA                : " << secs_remain    << "s \n";
    params::progress_file << "    Trial energy       : " << triale_red     << " Hartree\n";
    params::progress_file << "    <V>                : " << potential_red  << " Hartree\n";
    params::progress_file << "    Population         : " << population_red
                          << " (" << population_red/params::np << " per process) "<< "\n";
    params::progress_file << "    Cancelled weight   : " << cancelled_weight << "\n";

    if (params::dmc_iteration == 1)
    {
        // Before the first iteration, output names of the
        // evolution file columns
        params::evolution_file
                << "Population,"
                << "Trial energy,"
                << "<V>,"
                << "Cancelled weight,"
                << "Average weight\n";
    }

    // Output evolution information to file
    params::evolution_file
        << population_red    << ","
        << triale_red        << ","
        << potential_red     << ","
        << cancelled_weight  << ","
        << av_weight_red     << "\n";

    // Write the wavefunction to file
    if (params::write_wavefunction)
    {
        params::wavefunction_file << "# Iteration " << params::dmc_iteration << "\n";
        for (unsigned n=0; n<size(); ++n)
        {
            (*this)[n]->write_wavefunction();
            params::wavefunction_file << "\n";
        }
    }

    // Flush output files after every call
    params::flush();
}
