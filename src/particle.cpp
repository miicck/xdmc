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
#include "random.h"

int particle::constructed_count;

particle::particle()
{
	++ constructed_count;
	coords = new double[simulation.dimensions];
}

particle::~particle()
{
	-- constructed_count;
	delete[] coords;
}

particle* particle :: copy()
{
	particle* p = new particle();

	p->half_spins = this->half_spins;
	p->charge     = this->charge;
	p->mass       = this->mass;

	for (int i=0; i<simulation.dimensions; ++i)
		p->coords[i] = this->coords[i];

	return p;
}

int particle :: exchange_symmetry(particle* other)
{
	// Returns 
	//     0 for non-identical particles
	//     1 for identical bosons
	//    -1 for identical fermions

	// Check if identical (same quantum numbers)
	if (this->half_spins != other->half_spins)      return 0;
	if (fabs(this->mass   - other->mass)   > 10e-4) return 0;
	if (fabs(this->charge - other->charge) > 10e-4) return 0;

	// Return exchange symmetry
	if (this->half_spins % 2 == 0) return 1; // Boson
	return -1; // Fermion
}

double coulomb(particle* a, particle* b)
{
	if (fabs(a->charge) < 10e-10) return 0;
	if (fabs(b->charge) < 10e-10) return 0;

	double r = 0;
	for (int i=0; i<simulation.dimensions; ++i)
	{
		double dxi = a->coords[i] - b->coords[i];
		r += dxi * dxi;
	}
	r = sqrt(r);
	return (a->charge*b->charge)/r;
}

double particle::interaction(particle* other)
{
	return coulomb(this, other);
}

double particle :: overlap_prob(particle* clone)
{
	double r = 0;
	for (int i=0; i<simulation.dimensions; ++i)
	{
		double dxi = this->coords[i] - clone->coords[i];
		r += dxi * dxi;
	}
	return exp(-r/(4*simulation.tau));
}

void particle :: diffuse(double tau)
{
	// Diffuse the particle by moving each
	// coordinate by an amount sampled from
	// a normal distribution with variance tau.
	for (int i=0; i<simulation.dimensions; ++i)
		this->coords[i] += rand_normal(tau);
}

void particle :: sample_wavefunction()
{
	// Output my coordinates in the form x1,x2, ... xn
	// (with no trailing comma)
	for (int  i=0; i<simulation.dimensions; ++i)
	{
		simulation.wavefunction_file << this->coords[i];
		if (i < simulation.dimensions - 1)
			simulation.wavefunction_file << ",";
	}
	simulation.wavefunction_file << ";";
}
