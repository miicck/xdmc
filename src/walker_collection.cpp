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

void expectation_values :: reset()
{
    // Reset the expectation values ready
    // for accumulation
    cancellation_amount = 0;
    clipped_weight    = 0;
    average_potential = 0;
}

void expectation_values :: normalize(unsigned walker_count)
{
    // After accumulating expectation values
    // renormalize them
    double wc = double(walker_count);
    average_potential /= wc;
}

walker_collection :: walker_collection()
{
    // Reserve a reasonable amount of space to deal efficiently
    // with the fact that the population can fluctuate.
    walkers.reserve(int(params::target_population*
                        params::max_pop_ratio)+10);

    // Initialize the set of walkers to the target
    // population size.
    for (int i=0; i<params::target_population; ++i)
    {
        walker* w = new walker();
        walkers.push_back(w);

        // Carry out initial diffusion to avoid
        // exact particle overlap on first iteration
        w->diffuse(params::pre_diffusion);
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
    double exponent     = -params::tau*(av_potential - params::trial_energy);
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
    // Adjust weights if expected population after
    // branching is unacceptable
    double amw        = this->average_mod_weight();
    double ex_weight  = amw*size();
    double max_weight = params::target_population * params::max_pop_ratio;
    double min_weight = params::target_population * params::min_pop_ratio;
    if (ex_weight < min_weight || ex_weight > max_weight)
        for (int n=0; n<size(); ++n)
        {
            // Record the amount of weight "clipped" to keep the
            // population between target_population * min_pop_ratio
            // and target_population * max_pop_ratio
            walker* wn        = (*this)[n];
            double  w_before  = wn->weight;
            wn->weight       /= amw;
            this->expect_vals.clipped_weight += w_before - wn->weight;
        }

    // Carry out diffusion and branching of the walkers
    int nmax = size();
    for (int n=0; n < nmax; ++n)
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

        // Apply branching step, adding branched
        // survivors to the end of the collection
        int surviving = branch_from_weight(w->weight);
        for (int s=0; s<surviving; ++s)
            walkers.push_back(w->branch_copy());

        // Accumulate expectation values
        expect_vals.average_potential += pot_after * surviving;
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
    // Apply the selected cancellation scheme
    if      (params::cancel_scheme == "voronoi")
        this->apply_voronoi_cancellations();
    else if (params::cancel_scheme == "pairwise")
        this->apply_pairwise_cancellations();
    else if (params::cancel_scheme == "diffusive")
        this->apply_diffusive_cancellations();
    else if (params::cancel_scheme == "none")
        return;
    else
        params::error_file << "Unknown cancellation scheme: "
                              << params::cancel_scheme << "\n";

}

void walker_collection :: apply_diffusive_cancellations()
{
    double* new_weights = new double[size()];

    for (int n=0; n<size(); ++n)
    {
        walker* wn = (*this)[n];

        // Evaluate total diffusive Greens function
        // at the configuration of wn
        // gf = \sum_m w_m G_D(x_n, x_m, dt)
        double gf = 0;
        for (int m=0; m<size(); ++m)
        {
            if (m == n) continue;
            walker* wm = (*this)[n];
            gf += wm->diffusive_greens_function(wn) * wm->weight;
        }

        // Kill the walker if it ended up
        // in an opposite-sign region
        if (sign(wn->weight) == sign(gf))
            new_weights[n] = 1.0;
        else
            new_weights[n] = 0.0;
    }

    // Apply new weights, tracking the amount cancelled
    expect_vals.cancellation_amount = 0;
    for (int n=0; n<size(); ++n)
    {
        walker* wn = (*this)[n];
        expect_vals.cancellation_amount += fabs(wn->weight - new_weights[n]);
        wn->weight = new_weights[n];
    }

    // Free memory
    delete[] new_weights;
}

void walker_collection :: apply_voronoi_cancellations()
{
    // Allocate space for new weights
    double* new_weights = new double[size()];

    // Calculate new weights
    int nn = params::dimensions + 1;
    if (size() < 1 + nn)
    {
        params::error_file << "Error in voronoi: population very small!\n";
        return;
    }

    for (int n=0; n<size(); ++n)
    {
        walker* wn = (*this)[n];

        // Find the d+1 nearest neighbours
        double* sq_distances = new double[nn];
        int*    nn_indicies  = new int[nn];

        for (int i=0; i<nn; ++i)
        {
            nn_indicies[i]  = -1;
            sq_distances[i] = INFINITY;
        }

        for (int m=0; m<size(); ++m)
        {
            if (m == n) continue;

            walker* wm = (*this)[m];
            double sq_dis = wn->sq_distance_to(wm);
            
            for (int i=0; i<nn; ++i)
                if (sq_dis < sq_distances[i])
                {
                    // Shift records along by one to make
                    // space for new entry
                    for (int j=nn-1; j>i; --j)
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

        for (int i=0; i<nn; ++i)
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
            if (positive_weight > negative_weight + 1)
                new_weights[n] = 0;
        }
        else
        {
            // Walker is +ve, kill it if -ve is winning
            if (negative_weight > positive_weight + 1)
                new_weights[n] = 0;
        }

        // Free memory
        delete[] sq_distances;
        delete[] nn_indicies;
    }

    // Don't apply cancellations if a large fraction
    // of the walkers will die
    double total_weight = 0;
    for (int n=0; n<size(); ++n)
        total_weight += fabs(new_weights[n]);
    double frac_lost = 1 - total_weight/double(size());
    if (frac_lost > 0.5) return;

    // Apply new weights, tracking the amount cancelled
    expect_vals.cancellation_amount = 0;
    for (int n=0; n<size(); ++n)
    {
        walker* wn = (*this)[n];
        expect_vals.cancellation_amount += fabs(wn->weight - new_weights[n]);
        wn->weight = new_weights[n];
    }

    // Free memory
    delete[] new_weights;
}

void walker_collection :: apply_pairwise_cancellations()
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
    expect_vals.cancellation_amount = 0;
    for (int n=0; n<size(); ++n)
    {
        double delta = (*this)[n]->weight - weights_before[n];
        expect_vals.cancellation_amount += fabs(delta);
    }

    // Clean up memory
    delete[] weights_before;
}

void walker_collection :: correct_seperations()
{
    // Correct walker seperations in a pairwise manner
    for (int n=0; n<size(); ++n)
        for (int m=0; m<n; ++m)
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

void walker_collection :: write_output(int iter)
{
    // Sum various things across processes
    double population_red    = mpi_sum(double(this->size()));
    double cancel_red        = mpi_sum(this->expect_vals.cancellation_amount);

    // Average various things across processes
    double triale_red        = mpi_average(params::trial_energy);
    double av_pot_red        = mpi_average(this->expect_vals.average_potential);
    double av_weight_red     = mpi_average(this->average_weight());

    // Calculate timing information
    double time_per_iter     = params::time()/iter;
    double secs_remain       = time_per_iter * (params::dmc_iterations - iter);
    double percent_complete  = double(100*iter)/params::dmc_iterations;

    // Output iteration information
    params::progress_file << "\nIteration " << iter << "/" << params::dmc_iterations
                             << " (" << percent_complete << "%)\n";
    params::progress_file << "    Time running    : " << params::time()
                             << "s (" << time_per_iter << "s/iter)\n";
    params::progress_file << "    ETA             : " << secs_remain    << "s \n";
    params::progress_file << "    Trial energy    : " << triale_red     << " Hartree\n";
    params::progress_file << "    <V>             : " << av_pot_red     << " Hartree\n";
    params::progress_file << "    Population      : " << population_red
                             << " (" << population_red/params::np << " per process) "<< "\n";
    params::progress_file << "    Canceled weight : " << cancel_red     << "\n";

    if (iter == 1)
    {
        // Before the first iteration, output names of the
        // evolution file columns
        params::evolution_file
                << "Population,"
                << "Trial energy (Hartree),"
                << "<V> (Hartree),"
                << "Average weight,"
                << "Cancelled weight\n";
    }

    // Output evolution information to file
    params::evolution_file
        << population_red    << ","
        << triale_red        << ","
        << av_pot_red        << ","
        << av_weight_red     << ","
        << cancel_red        << "\n";

    // Write the wavefunction to file
    if (params::write_wavefunction)
    {
        params::wavefunction_file << "# Iteration " << iter << "\n";
        for (int n=0; n<size(); ++n)
        {
            (*this)[n]->write_wavefunction();
            params::wavefunction_file << "\n";
        }
    }

    // Flush output files after every call
    params::flush();
}
