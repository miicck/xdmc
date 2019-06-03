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

#include "random.h"
#include "walker.h"
#include "simulation.h"

//%%%%%%%%//
// walker //
//%%%%%%%%//

int walker :: constructed_count = 0;

walker :: walker()
{
	// Default constructor:
	// Setup the walker with the particles
	// that describe the system.
	++ constructed_count;
	this->weight = 1;
	this->particles.clear();
	for (int i=0; i<simulation.template_system.size(); ++i)
		this->particles.push_back(simulation.template_system[i]->copy());
}

walker :: walker(std::vector<particle*> particles)
{
	// Constructor to create a walker explicitly from
	// a collection of particles
	++ constructed_count;
	this->weight = 1;
	this->particles = particles;
}

walker :: ~walker()
{
	// Clear up memory (delete all the particles).
	-- constructed_count;
	for (int i=0; i<particles.size(); ++i)
		delete(particles[i]);
}


walker* walker :: copy()
{
	// Return a copy of this walker (copy each of the particles)
	std::vector<particle*> copied_particles;
	for (int i=0; i<particles.size(); ++i)
		copied_particles.push_back(particles[i]->copy());
	walker* copy = new walker(copied_particles);
	copy->weight = this->weight;
	return copy;
}

double walker :: potential()
{
	// No need to reevaluate the potential
	if (!potential_dirty)
		return last_potential;
	
	// Evaluate the potential of the system
	// in the configuration described by this walker
	last_potential = 0;
	for (int i = 0; i < particles.size(); ++i)
	{
		// Sum up external potential contributions
		for (int j=0; j<simulation.potentials.size(); ++j)
			last_potential += simulation.potentials[j]->potential(particles[i]);

		// Particle-particle interactions
		// note j<i => no double counting
		for (int j=0; j<i; ++j)
			last_potential += particles[i]->interaction(particles[j]);
	}

	potential_dirty = false;
	return last_potential;
}

void walker :: diffuse(double tau=simulation.tau)
{
	// Diffuse all of the particles (classical
	// particles will automatically not diffuse)
	for (int i=0; i<particles.size(); ++i)
		particles[i]->diffuse(tau);
	
	// Particles have moved => potential has changed
	potential_dirty = true;
}

void walker :: exchange()
{
	// Apply random exchange moves to particles
	// Note: because only identical particles
	// are exchanged, the potential remains the
	// same => we do not need to set the potential_dirty
	// flag.

	// No exchanges possible
	if (simulation.exchange_values.size() == 0) return;

	// Apply exchanges stochastically
	if (rand_uniform() < 0.5) return;

	// Pick a random exchangable pair
	int i = rand() % simulation.exchange_values.size();
	particle* p1 = particles[simulation.exchange_pairs[2*i]];
	particle* p2 = particles[simulation.exchange_pairs[2*i+1]];

	// Exchange them 
	this->weight *= double(simulation.exchange_values[i]);
	for (int c=0; c<simulation.dimensions; ++c)
	{
		double tmp = p1->coords[c];
		p1->coords[c] = p2->coords[c];
		p2->coords[c] = tmp;
	}
}

void walker :: cancel(walker* other)
{
	// Apply cancellation of two walkers

	// Cancellation is only necassary if there is
	// a fermionic exchange in the system (as this is
	// the only way that walker signs can change)
	if (!simulation.has_fermionic_exchange) return;

	// Caclulate the probability that these 
	// two walkers would overlap in the next iteration
	double p = 1.0;
	for (int i=0; i<particles.size(); ++i)
		p *= particles[i]->overlap_prob(other->particles[i]);

	// With this probability, cancel the walkers
	if (rand_uniform() < p)
	{
		double av_weight = (this->weight + other->weight)/2.0;
		this->weight  = av_weight;
		other->weight = av_weight;
	}
}

void walker :: sample_wavefunction()
{
	simulation.wavefunction_file << this->weight << ":";
	for (int i=0; i<particles.size(); ++i)
		particles[i]->sample_wavefunction();
	simulation.wavefunction_file << "\n";
}

//%%%%%%%%%%%%%%%%%%%//
// walker_collection //
//%%%%%%%%%%%%%%%%%%%//

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
		w->exchange();
		w->diffuse(1.0);
	}
}

walker_collection :: ~walker_collection()
{
	// Clean up memory
	for (int n=0; n<size(); ++n)
		delete (*this)[n];
}

int walkers_surviving(double pot_before, double pot_after, double weight)
{
	// Returns the number of walkers that should survive after a
	// walker with the given weight moves from potential pot_before to pot_after.
	// Returns 0 if the walker should die and 1 if it should live, 2 if another should
	// spawn etc...

        if (std::isinf(pot_after))  return 0;
        if (std::isinf(pot_before)) return 0;

        double av_potential = (pot_before + pot_after)/2;
        double exponent     = -simulation.tau*(av_potential - simulation.trial_energy);
        double p            = fabs(weight) * exp(exponent);
        int surviving       = std::min(int(floor(p+rand_uniform())), 3);

        return surviving;
}

void walker_collection :: diffuse_and_branch()
{
	// Carry out diffusion and branching of the walkers
	// (this is done as one step because we need to know
	// the potential before and after diffusion to apply
	// the branching step)

	// Normalize the walkers so that the population stays
	// roughly constant
	this->renormalize();

	int nmax = size();
	for (int n=0; n < nmax; ++n)
	{
		walker* w = (*this)[n];

		double pot_before = w->potential();
		w->diffuse();
		double pot_after = w->potential();

		int surviving = walkers_surviving(pot_before, pot_after, w->weight);

		if (surviving == 0)
		{
			// Delete this walker
			// (swap it to the end first for efficient
			// removal)
			std::swap(walkers[n], walkers.back());
			walkers.pop_back();
			delete(w);

			// Deal with the fact this messed up the
			// thing we were iterating over
			--n;
			--nmax;
		}

		// Add neccasary copies of this walker at the end of the collection
		for (int s=0; s<surviving-1; ++s)
			walkers.push_back(w->copy());
	}
}

double walker_collection :: average_weight_sq()
{
	// Returns (1/N) * sum_i |w_i|^2
	double av_weight_sq = 0;
	for (int n=0; n<size(); ++n)
		av_weight_sq += (*this)[n]->weight*(*this)[n]->weight;
	av_weight_sq /= double(size());
	return av_weight_sq;
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

double walker_collection :: energy_estimate()
{
	// Returns (1/N) * sum_i |w_i|*v_i
	// Where v_i is the potential energy 
	// of the i^th walker.
	double energy = 0;
	for (int n=0; n<size(); ++n)
		energy += (*this)[n]->potential();
	energy /= double(size());
	return energy;
}

void walker_collection :: renormalize()
{
	// Renormalize the walker weights such that
	// sum_i |w_i|^2 = 1
	double norm = sqrt(this->average_weight_sq());
	for (int n=0; n<size(); ++n)
		(*this)[n]->weight /= norm;
}

void walker_collection :: sample_wavefunction()
{
	// Sample the wavefunction to a file
	for (int n=0; n<size(); ++n)
		(*this)[n]->sample_wavefunction();
	simulation.wavefunction_file << "#" << "\n";
}

void walker_collection :: make_exchange_moves()
{
	// Apply exchange moves to each of the walkers
	for (int n=0; n<size(); ++n)
		(*this)[n]->exchange();
}

void walker_collection :: apply_cancellations()
{
	// Apply cancellations to walker weights
	for (int n=0; n<size(); ++n)
		for (int m=0; m<n; ++m)
			(*this)[n]->cancel((*this)[m]);
}
