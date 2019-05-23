#include <cmath>

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
