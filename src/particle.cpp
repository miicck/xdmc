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

#include <iostream>

#include "simulation.h"
#include "particle.h"
#include "random.h"

int particle::count;

double particle::interaction(particle* other)
{
	// Coulomb
	double r = 0;
	for (int i=0; i<simulation.dimensions; ++i)
	{
		double dxi = this->coords[i] - other->coords[i];
		r += dxi * dxi;
	}
	r = sqrt(r);
	return (this->charge()*other->charge())/r;
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

particle::particle()
{
	++ count;
	coords = new double[simulation.dimensions];
}

particle::~particle()
{
	-- count;
	delete[] coords;
}

void quantum_particle :: diffuse(double tau)
{
	// Diffuse the particle by moving each
	// coordinate by an amount sampled from
	// a normal distribution with variance tau.
	for (int i=0; i<simulation.dimensions; ++i)
		this->coords[i] += rand_normal(tau);
}

void quantum_particle :: sample_wavefunction()
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

particle* electron :: copy()
{
	// Copy this electron.
	electron* ret = new electron();
	for (int i=0; i<simulation.dimensions; ++i)
		ret->coords[i] = this->coords[i];
	return ret;
}

particle* non_interacting_fermion :: copy()
{
	particle* ret = new non_interacting_fermion();
	for (int i=0; i<simulation.dimensions; ++i)
		ret->coords[i] = this->coords[i];
	return ret;
}

nucleus :: nucleus(double atomic_number, double mass_number)
{
	this->atomic_number = atomic_number;
	this->mass_number   = mass_number;
}

particle* nucleus :: copy()
{
	// Copy this nucleus
	nucleus* ret = new nucleus(atomic_number, mass_number);
	for (int i=0; i<simulation.dimensions; ++i)
		ret->coords[i] = this->coords[i];
	return ret;
}


