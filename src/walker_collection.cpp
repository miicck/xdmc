/*

    XDMC
    Copyright (C) Michael Hutcheon (email mjh261@cam.ac.uk)

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
#include <cmath>
#include <iostream>
#include <mpi.h>

#include "catch.h"
#include "random.h"
#include "dmc_math.h"
#include "walker_collection.h"
#include "params.h"
#include "mpi_utils.h"
#include "utils.h"

bool walker_collection :: propagate(walker_collection* walkers_last)
{
    // Apply the stages of walker propagation
    // returns false if this iteration should be
    // reverted, because of population explosion etc...

    // Diffusive moves involving G_D
    make_diffusive_moves(walkers_last);

    // Exchange moves
    make_exchange_moves();

    // Renormalize by applying exp(E_T \delta\tau)
    // (work out E_T as well)
    apply_renormalization();

    // Check for population explosion
    for (unsigned n=0; n<walkers.size(); ++n)
        if (fabs(walkers[n]->weight) > params::max_weight)
            return false;

    branch();

    // Check for population collapse
    if (walkers.size() == 0)
        return false;

    return true;
}

void walker_collection :: make_diffusive_moves(walker_collection* walkers_last)
{
    // Reset things
    if (params::write_nodal_surface)
        params::nodal_surface_file << "# Iteration " << params::dmc_iteration << "\n";
    params::cancelled_weight = 0.0;
    
    // Carry out the specified diffusion scheme
    if (params::diffusion_scheme == "exact_1d")
        diffuse_exact_1d();
    else if (params::diffusion_scheme == "max_seperation")
        diffuse_max_seperation(walkers_last);
    else if (params::diffusion_scheme == "max_seperation_mpi")
        diffuse_max_seperation_mpi(walkers_last);
    else if (params::diffusion_scheme == "stochastic_nodes")
        diffuse_stochastic_nodes(walkers_last);
    else if (params::diffusion_scheme == "stochastic_nodes_permutations")
        diffuse_stochastic_nodes_permutations(walkers_last);
    else if (params::diffusion_scheme == "stochastic_nodes_mpi")
        diffuse_stochastic_nodes_mpi(walkers_last);
    else if (params::diffusion_scheme == "exchange_diffuse")
        exchange_diffuse(walkers_last);
    else if (params::diffusion_scheme == "bosonic")
        diffuse_bosonic(walkers_last);
    else if (params::diffusion_scheme == "stay_irreducible")
        diffuse_in_irreducible(walkers_last);
    else
        throw "Unkown diffusion scheme";
}

void walker_collection :: estimate_tau_nodes()
{
    // Estimate the new value for tau_nodes
    double new_tau = params::tau_nodes;
    if (params::tau_nodes_estimator == "min_sep") new_tau = tau_nodes_min_sep();
    else if (params::tau_nodes_estimator == "min_sep_mpi") new_tau = tau_nodes_min_sep_mpi();
    else if (params::tau_nodes_estimator != "none") throw "Unkown tau_nodes estimator!";

    // Don't propagate NaN, or Inf values
    if (std::isfinite(new_tau))
        params::tau_nodes = new_tau;

    params::tau_nodes = mpi_average(params::tau_nodes);
}

void walker_collection :: make_exchange_moves()
{
    if (params::correct_average_weight)
    {
        // Work out the probability for +ve and -ve walkers
        // to change sign respectively, in order to keep the
        // average sign close to 0
        double positive_weight = 0;
        double negative_weight = 0;
        for (unsigned n=0; n<walkers.size(); ++n)
        {
            walker* w = walkers[n];
            if (w->weight > 0) positive_weight += w->weight;
            else negative_weight -= w->weight;
        }

        double positive_exchange_prob = (positive_weight - negative_weight)/positive_weight;
        double negative_exchange_prob = (negative_weight - positive_weight)/negative_weight;
        if (positive_exchange_prob < 0) positive_exchange_prob = 0;
        if (negative_exchange_prob < 0) negative_exchange_prob = 0;

        // Damp the rate of correction, so we don't end up
        // in an oscillatory regime
        positive_exchange_prob *= 0.9;
        negative_exchange_prob *= 0.9;

        // Apply exchange moves with the above probabilities
        for (unsigned n=0; n<walkers.size(); ++n)
        {
            walker* w = walkers[n];
            double prob = w->weight < 0 ? negative_exchange_prob : positive_exchange_prob;
            if (rand_uniform() < prob)
                w->change_sign();
        }
    }

    // Apply exchange moves to each of the walkers
    for (unsigned n=0; n<walkers.size(); ++n)
        walkers[n]->exchange();
}

double walker_collection :: diffused_wavefunction(
    walker* c, double tau=params::tau, int self_index=-1)
{
    // Evaluate the diffused wavefunction 
    // at the configuration of c: 
    // \psi_D(c) = \sum_i w_i G_D(c, x_i, dt)
    double psi_d = 0;
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];

        // Treat my own contribution to the
        // greens function differently
        double  amp = 1.0;
        if (int(n) == self_index) 
            amp = params::self_gf_strength;

        psi_d += amp * w->weight * w->diffusive_greens_function(c, tau);
    }
    return psi_d;
}

double* walker_collection :: diffused_wavefunction_signed(
    walker* c, double tau=params::tau, int self_index=-1)
{
    // Evaluate the diffused wavefunction as components whos sign
    // matches c and those that do not
    double* ret = new double[2];
    ret[0] = 0; // Same sign
    ret[1] = 0; // Opposite sign
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];

        // Treat my own contribution to the
        // greens function differently
        double  amp = 1.0;
        if (int(n) == self_index) 
            amp = params::self_gf_strength;

        double gf = amp * fabs(w->weight) * w->diffusive_greens_function(c, tau);
        if (sign(w->weight/c->weight) == 1) ret[0] += gf;
        else ret[1] += gf;
    }
    return ret;
}

double* walker_collection :: exchange_diffused_wfn_signed(
    walker* c, double tau=params::tau, int self_index=-1)
{
    // Evaluate the exchange-diffused wavefunction as components whos sign
    // matches c and those that do not
    double* ret = new double[2];
    ret[0] = 0; // Same sign
    ret[1] = 0; // Opposite sign
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w  = walkers[n];

        // Treat my own contribution to the
        // greens function differently
        double  amp = 1.0;
        if (int(n) == self_index) 
            amp = params::self_gf_strength;

        double* gf = w->exchange_diffusive_gf(c, tau);
        gf[0]     *= amp*fabs(w->weight);
        gf[1]     *= amp*fabs(w->weight);

        if (sign(w->weight/c->weight) == 1)
        {
            // w and c are same sign =>
            // gf[0] is same-sign contribution
            // gf[1] is opposite-sign contribution
            ret[0] += gf[0];
            ret[1] += gf[1]; 
        }
        else
        {
            // w and c are opposite signs =>
            // gf[0] is opposite-sign contribution
            // gf[1] is same-sign contribution
            ret[0] += gf[1]; 
            ret[1] += gf[0];
        }

        // Free memory
        delete[] gf;
    }
    return ret;
}

double potential_greens_function(double pot_before, double pot_after)
{
    // Evaluate the potential-dependent part of the greens
    // function G_v(x,x',tau) = exp(-tau*(v(x)+v(x'))/2)
    if (std::isinf(pot_after))  return 0;
    if (std::isinf(pot_before)) return 0;
    return fexp( -params::tau * (pot_before + pot_after)/2.0 );
}

void walker_collection :: diffuse_exact_1d()
{
    // Error if dimensions of system != 1
    if (params::dimensions != 1)
        throw "Dimension != 1 in exact 1d diffusion!";

    // Carry out diffusion of the walkers
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        // Diffuse the walker
        walker* w = walkers[n];
        walker* w_before = w->copy();
        w->diffuse(params::tau);

        // Kill walkers crossing the nodal surface
        if (w_before->crossed_nodal_surface(w))
        {
            // Record the nodal surface
            if (params::write_nodal_surface)
                w->write_coords(params::nodal_surface_file);

            // Kill the walker
            w->weight = 0;
            params::cancelled_weight += 1;
        }

        // Apply potential part of greens function
        double pot_before = w_before->potential();
        double pot_after  = w->potential();
        w->weight        *= potential_greens_function(pot_before, pot_after);

        delete w_before;
    }
}

void walker_collection :: diffuse_bosonic(walker_collection* walkers_last)
{
    // Carry out normal bosonic DMC diffusion
    // where each walker diffuses independently
    // according to the diffusive greens function
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];
        w->diffuse(params::tau);

        // Apply potential part of greens function
        double pot_before = walkers_last->walkers[n]->potential();
        double pot_after  = w->potential();
        w->weight        *= potential_greens_function(pot_before, pot_after);
    }
}

void walker_collection :: diffuse_max_seperation(walker_collection* walkers_last)
{
    // Carry out diffusion of the walkers in a manner
    // that will result in the maximum seperation of 
    // +ve wlakers to -ve walkers.
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];
        w->diffuse(params::tau);
        double* psi = walkers_last->diffused_wavefunction_signed(w, params::tau, int(n));
        double* psi_nodes = walkers_last->diffused_wavefunction_signed(w, params::tau_nodes, int(n));

        if (psi[0] < psi[1] || psi_nodes[0] < psi_nodes[1])
        {
            // Record the nodal surface
            if (params::write_nodal_surface)
                w->write_coords(params::nodal_surface_file);

            // w has strayed into the wrong neighbourhood, kill them
            w->weight = 0;
            params::cancelled_weight += 1;
        }
        else
        {
            // Account for cancellations
            params::cancelled_weight += fabs(w->weight * psi[1]/psi[0]);
            w->weight *= 1 - psi[1]/psi[0];
        }

        // Free memory
        delete[] psi;
        delete[] psi_nodes;

        // Apply potential part of greens function
        double pot_before = walkers_last->walkers[n]->potential();
        double pot_after  = w->potential();
        w->weight        *= potential_greens_function(pot_before, pot_after);
    }
}

void walker_collection :: diffuse_max_seperation_mpi(walker_collection* walkers_last)
{
    // Carry out diffusion of walkers, evaluating a stochastic nodal
    // surface using all the walkers across processes

    // Loop over processes
    for (int pid=0; pid<params::np; ++pid)
    {
        // Get the number of walkers on this process
        int walker_count = params::pid == pid ? walkers.size() : 0;
        MPI_Bcast(&walker_count, 1, MPI_INT, pid, MPI_COMM_WORLD);
        
        // Propagate the walkers on this process
        for (int n=0; n<walker_count; ++n)
        {
            // w is the n^th walker on the pid^th process
            walker* w = params::pid == pid ? walkers[n] : nullptr;

            // On the pid^th process, diffuse the walker
            if (params::pid == pid)
                w->diffuse(params::tau);

            // Get a copy of the n^th walker on the 
            // pid^th process, after diffusion.
            // (the copy will be valid on all processes)
            walker* w_after = walker::mpi_copy(w, pid);

            // Compute the across-process wavefunction after diffusion
            double* psi_after_pid = walkers_last->
                diffused_wavefunction_signed(w_after, params::tau_nodes, -1);
            double* psi_after = new double[2];
            MPI_Reduce(psi_after_pid, psi_after, 2, MPI_DOUBLE, MPI_SUM, pid, MPI_COMM_WORLD);

            // On the pid^th process, apply the cancellation function
            // w -> w * f_+/-
            if (params::pid == pid)
            {
                if (psi_after[0] < psi_after[1])
                {
                    // Record the nodal surface
                    if (params::write_nodal_surface)
                        w->write_coords(params::nodal_surface_file);

                    // Kill the walker
                    w->weight = 0;
                    params::cancelled_weight += 1;
                }
                else
                {
                    params::cancelled_weight += fabs(w->weight * psi_after[1]/psi_after[0]);
                    w->weight *= 1 - psi_after[1]/psi_after[0];
                }
            }

            // Free memory
            delete w_after;
            delete[] psi_after_pid;
            delete[] psi_after;
        }
    }

    // Apply the potential part of the greens function
    // (which is independent of the other processes)
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        double pot_before   = walkers_last->walkers[n]->potential();
        double pot_after    = walkers[n]->potential();
        walkers[n]->weight *= potential_greens_function(pot_before, pot_after);
    }
}


void walker_collection :: exchange_diffuse(walker_collection* walkers_last)
{
    // Carry out diffusion of the walkers in a manner
    // that will result in the maximum seperation of 
    // +ve walkers to -ve walkers, taking into account
    // the exchanged images of the walkers.
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];
        w->diffuse(params::tau);
        double* psi = walkers_last->exchange_diffused_wfn_signed(w, params::tau, int(n));

        if (psi[0] < psi[1])
        {
            // w has strayed into the wrong neighbourhood, kill them

            // Record the nodal surface
            if (params::write_nodal_surface)
                w->write_coords(params::nodal_surface_file);

            // Kill the walker
            w->weight = 0;
            params::cancelled_weight += 1;
        }
        else
        {
            // Account for cancellations
            params::cancelled_weight += fabs(w->weight * psi[1]/psi[0]);
            w->weight *= 1 - psi[1]/psi[0];
        }

        // Free memory
        delete psi;

        // Apply potential part of greens function
        double pot_before = walkers_last->walkers[n]->potential();
        double pot_after  = w->potential();
        w->weight        *= potential_greens_function(pot_before, pot_after);
    }
}

void walker_collection :: diffuse_stochastic_nodes(walker_collection* walkers_last)
{
    // Carry out diffusion of walkers, killing any that cross the
    // stochastic nodal surface set up last iteration
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];

        double psi_before = walkers_last->
            diffused_wavefunction(w, params::tau_nodes, int(n));

        w->diffuse(params::tau);

        double psi_after  = walkers_last->
            diffused_wavefunction(w, params::tau_nodes, int(n));

        if (sign(psi_before) != sign(psi_after))
        {
            // Record the nodal surface
            if (params::write_nodal_surface)
                w->write_coords(params::nodal_surface_file);

            // w has strayed into the wrong neighbourhood, kill them
            w->weight = 0;
            params::cancelled_weight += 1;
        }

        // Apply potential part of greens function
        double pot_before = walkers_last->walkers[n]->potential();
        double pot_after  = w->potential();
        w->weight        *= potential_greens_function(pot_before, pot_after);
    }
}

void walker_collection :: diffuse_stochastic_nodes_mpi(walker_collection* walkers_last)
{
    // Carry out diffusion of walkers, evaluating a stochastic nodal
    // surface using all the walkers across processes

    // Loop over processes
    for (int pid=0; pid<params::np; ++pid)
    {
        // Get the number of walkers on this process
        int walker_count = params::pid == pid ? walkers.size() : 0;
        MPI_Bcast(&walker_count, 1, MPI_INT, pid, MPI_COMM_WORLD);
        
        // Propagate the walkers on this process
        for (int n=0; n<walker_count; ++n)
        {
            // Get a copy of the i^th walker on the 
            // pid^th process, before diffusion
            walker* w = params::pid == pid ? walkers[n] : nullptr;
            walker* w_before = walker::mpi_copy(w, pid);

            // Compute the across-process wavefunction before diffusion
            double psi_before_pid = walkers_last->
                diffused_wavefunction(w_before, params::tau_nodes, -1); 
            double psi_before;
            MPI_Reduce(&psi_before_pid, &psi_before, 1, MPI_DOUBLE, MPI_SUM, pid, MPI_COMM_WORLD);

            // On the pid^th process, diffuse the walker
            if (params::pid == pid)
                w->diffuse(params::tau);

            // Get a copy of the i^th walker on the 
            // pid^th process, after diffusion
            walker* w_after = walker::mpi_copy(w, pid);

            // Compute the across-process wavefunction after diffusion
            double psi_after_pid = walkers_last->
                diffused_wavefunction(w_after, params::tau_nodes, -1);
            double psi_after;
            MPI_Reduce(&psi_after_pid, &psi_after, 1, MPI_DOUBLE, MPI_SUM, pid, MPI_COMM_WORLD);

            // On the pid^th process, kill the walker if it crossed the nodal surface
            if (params::pid == pid)
                if (sign(psi_before) != sign(psi_after))
                {
                    // Record the nodal surface
                    if (params::write_nodal_surface)
                        w->write_coords(params::nodal_surface_file);

                    // Kill the walker
                    w->weight = 0;
                    params::cancelled_weight += 1;
                }

            // Free memory
            delete w_before;
            delete w_after;
        }
    }

    // Apply the potential part of the greens function
    // (which is independent of the other processes)
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        double pot_before   = walkers_last->walkers[n]->potential();
        double pot_after    = walkers[n]->potential();
        walkers[n]->weight *= potential_greens_function(pot_before, pot_after);
    }
}

void walker_collection :: diffuse_stochastic_nodes_permutations(walker_collection* walkers_last)
{
    // Carry out diffusion of walkers, killing any that cross the
    // stochastic nodal surface set up last iteration
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];

        w->diffuse(params::tau);

        double* psi_after  = walkers_last->
            exchange_diffused_wfn_signed(w, params::tau_nodes, int(n));

        if (psi_after[0] < psi_after[1])
        {
            // Record the nodal surface
            if (params::write_nodal_surface)
                w->write_coords(params::nodal_surface_file);

            // w has strayed into the wrong neighbourhood, kill them
            w->weight = 0;
            params::cancelled_weight += 1;
        }

        delete[] psi_after;

        // Apply potential part of greens function
        double pot_before = walkers_last->walkers[n]->potential();
        double pot_after  = w->potential();
        w->weight        *= potential_greens_function(pot_before, pot_after);
    }
}

void walker_collection :: diffuse_in_irreducible(walker_collection* walkers_last)
{
    // Carry out diffusion of walkers, killing any that leave 
    // the irreducible wedge of configuration space
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        walker* w = walkers[n];
        w->diffuse(params::tau);

        if (w->is_irreducible())
        {
            // Apply potential part of greens function
            double pot_before = walkers_last->walkers[n]->potential();
            double pot_after  = w->potential();
            w->weight        *= potential_greens_function(pot_before, pot_after);
        }
        else
            w->weight = 0;
    }
}

double walker_collection :: distance_to_nearest_opposite(walker* w)
{
    double min_dis = -1;

    // Returns the distance to the nearest walker in
    // this collection to w, which has opposite sign to w
    for (unsigned n=0; n < walkers.size(); ++n)
    {
        // These walkers are the same sign
        if (sign(walkers[n]->weight) == sign(w->weight))
            continue;

        // Record the distnace if this is the closest so far
        double dis = w->sq_distance_to(walkers[n]);
        if (min_dis < 0 || dis < min_dis)
            min_dis = dis;
    }

    return sqrt(min_dis);
}

double walker_collection :: tau_nodes_min_sep_mpi()
{
    // Estimate tau_nodes from the minimum seperation
    // between any +ve and any -ve walker
    double average_min_dis = 0;
    int population = 0;

    // Loop over processes
    for (int pid=0; pid<params::np; ++pid)
    {
        // Get the number of walkers on this process
        int walker_count = params::pid == pid ? walkers.size() : 0;
        MPI_Bcast(&walker_count, 1, MPI_INT, pid, MPI_COMM_WORLD);
        population += walker_count;
        
        // For each walker on process pid
        for (int n=0; n<walker_count; ++n)
        {
            // Get a copy of the i^th walker on the pid^th process
            walker* w = params::pid == pid ? walkers[n] : nullptr;
            walker* w_copy = walker::mpi_copy(w, pid);

            // Get the minimum distance between w and any walker of the
            // opposite sign, across processes
            double min_dis_this = distance_to_nearest_opposite(w_copy);
            double min_dis;
            MPI_Allreduce(&min_dis_this, &min_dis, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

            // Accumulate the average minimum-distance
            average_min_dis += min_dis;

            delete w_copy;
        }
    }

    return average_min_dis / (2.0 * double(population));
}

double walker_collection :: tau_nodes_min_sep()
{
    double average_min_dis = 0;

    for (unsigned n=0; n<walkers.size(); ++n)
    {
        walker* w = walkers[n];
        average_min_dis += distance_to_nearest_opposite(w);
    }

    average_min_dis /= double(walkers.size());
    return average_min_dis / 2.0;
}

void walker_collection :: apply_renormalization()
{
    // Apply the selected renormalization scheme
    // <=> energy estimator
    if (params::energy_estimator == "growth")
        renormalize_growth();
    else if (params::energy_estimator == "potential")
        renormalize_potential();
    else
        throw "Unkown energy estimator!";

    // Revert trial energies found to be nan/inf
    static double last_non_nan = 0;
    if (std::isnan(params::trial_energy) || std::isinf(params::trial_energy))
        params::trial_energy = last_non_nan;
    else
        last_non_nan = params::trial_energy;

    params::trial_energy = mpi_average(params::trial_energy);
}

unsigned target_population()
{
    // Get the target population, which potentially varies with iteration
    unsigned target_population = params::target_population;
    if (params::end_target_population > 0)
    {
        double frac = params::dmc_iteration / double(params::dmc_iterations);
        double pop  = params::target_population * (1-frac) + params::end_target_population * frac;
        target_population = int(pop);
    }
    return target_population;
}

void walker_collection :: renormalize_potential()
{
    // Set the trial energy with reference to the
    // potential energy
    params::trial_energy = mpi_average(average_potential());

    // Bias towards target population
    params::trial_energy -= log(mpi_sum(sum_mod_weight()) / target_population());

    // Apply normalization greens function
    double gn = fexp(params::trial_energy * params::tau);
    for (unsigned n=0; n<walkers.size(); ++n)
        walkers[n]->weight *= gn;
}

void walker_collection :: renormalize_growth()
{
    // Set the trial energy and apply renormalization
    // so that the population fluctuates around the 
    // target population. Do this by employing the 
    // growth estimator of the energy

    // The population at the start of the iteration
    double pop_before_propagation = mpi_sum(double(walkers.size()));

    // The effective population now, after the cumulative
    // effect of this iterations greens functions
    // (i.e cancellation, diffusion, potential etc...)
    double pop_after_propagation  = mpi_sum(sum_mod_weight());

    // Set trial energy to minimize fluctuations
    double new_trial_energy = log(pop_before_propagation / pop_after_propagation)/params::tau;

    // Bias towards target population
    new_trial_energy -= log(pop_before_propagation / target_population());

    params::trial_energy = params::trial_energy * params::growth_mixing_factor
                         + new_trial_energy * (1.0 - params::growth_mixing_factor);

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

walker_collection :: walker_collection()
{
    // Per-process target population
    unsigned per_process_pop = params::target_population / params::np;

    // Reserve a reasonable amount of space to deal efficiently
    // with the fact that the population can fluctuate.
    walkers.reserve(2*int(per_process_pop));

    // Initialize the set of walkers to the target population size.
    for (unsigned i=0; i<per_process_pop; ++i)
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

double walker_collection :: sum_mod_weight()
{
    // Returns sum_i |w_i|
    // This is the effective population
    double sum = 0;
    for (unsigned n=0; n<walkers.size(); ++n)
        sum += fabs(walkers[n]->weight);
    return sum;
}

double walker_collection :: positive_weight()
{
    // Returns the sum of all positive weights
    double positive_weight = 0;
    for (unsigned n=0; n<walkers.size(); ++n)
        if (walkers[n]->weight > 0)
            positive_weight += walkers[n]->weight;
    return positive_weight;
}

double walker_collection :: negative_weight()
{
    // Returns the modulus of the sum of all negative weights
    double negative_weight = 0;
    for (unsigned n=0; n<walkers.size(); ++n)
        if (walkers[n]->weight < 0)
            negative_weight -= walkers[n]->weight;
    return negative_weight;
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

void walker_collection :: write_output(bool reverted)
{
    // Sum various things across processes
    double population_red    = mpi_sum(double(walkers.size()));
    double canc_weight_red   = mpi_sum(params::cancelled_weight);
    int    reverted_red      = mpi_sum(int(reverted));
    double canc_weight_perc  = 100.0*canc_weight_red/double(population_red);
    double pos_weight_red    = mpi_sum(positive_weight());
    double neg_weight_red    = mpi_sum(negative_weight());
    double total_weight_red  = pos_weight_red - neg_weight_red;
    double av_weight_red     = total_weight_red / population_red;

    // Average various things across processes
    double triale_red        = mpi_average(params::trial_energy);
    double tau_nodes_red     = mpi_average(params::tau_nodes);

    // Calculate timing information
    double time_per_iter     = params::dmc_time()/params::dmc_iteration;
    double percent_complete  = double(100*params::dmc_iteration)/params::dmc_iterations;
    int secs_remain = int(time_per_iter * (params::dmc_iterations - params::dmc_iteration));

    // Output iteration information
    params::progress_file << "\nIteration " << params::dmc_iteration 
        << "/"                              << params::dmc_iterations
        << " ("                             << percent_complete
        << "% imaginary time = "       << params::tau*params::dmc_iteration << ")\n"
        << "    Time running       : " << params::time()
        << "s ("                       << time_per_iter                 << "s/iter)\n"
        << "    ETA                : " << seconds_to_human(secs_remain) << "\n"
        << "    Trial energy       : " << triale_red                    << " Hartree\n"
        << "    Population         : " << population_red                
        << " ("                        << population_red/params::np     << " per process) \n"
        << "    Cancelled weight   : " << canc_weight_red
        << " ("                        << canc_weight_perc              << "% of the total weight)\n"
        << "    Reverted on        : " << reverted_red 
        << "/"                         << params::np                    << " processes\n"
        << "    Nodal timestep     : " << tau_nodes_red                 << " a.u\n";

    if (params::dmc_iteration == 1)
    {
        // Before the first iteration, output names of the
        // evolution file columns
        params::evolution_file
                << "Population,"
                << "Trial energy,"
                << "Positive weight,"
                << "Negative weight,"
                << "Average weight,"
                << "Cancelled weight,"
                << "Tau_nodes\n";
    }

    // Output evolution information to file
    params::evolution_file
        << population_red                  << ","
        << triale_red                      << ","
        << pos_weight_red                  << ","
        << neg_weight_red                  << ","
        << av_weight_red                   << ","
        << canc_weight_red                 << ","
        << tau_nodes_red                   << "\n";

    // Write the wavefunction to file
    if (params::write_wavefunction)
    {
        params::wavefunction_file << "# Iteration " << params::dmc_iteration << "\n";
        for (unsigned n=0; n<walkers.size(); ++n)
            walkers[n]->write_coords(params::wavefunction_file);
    }

    // Flush output files after every call
    params::flush();
}

bool walker_collection :: compare(walker_collection* other_walkers)
{
    // Compare two collections of walkers, returns false if 
    // they differ in any way (for testing purposes)
    for (unsigned i=0; i<walkers.size(); ++i)
    {
        walker* w1 = walkers[i];
        walker* w2 = other_walkers->walkers[i];
        if (!w1->compare(w2))
            return false;
    }

    return true;
}

TEST_CASE("Walker collection tests", "[walker_collection]")
{
    // Create collection of walkers
    walker_collection* c = new walker_collection();

    SECTION("Copy method")
    {
        // Test the copy method
        walker_collection* c_copy = c->copy();
        REQUIRE(c->compare(c_copy));
        delete c_copy;
    }

    // Free memory
    delete c;
}
