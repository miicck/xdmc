#include <cmath>

#include "simulation.h"
#include "particle.h"

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
