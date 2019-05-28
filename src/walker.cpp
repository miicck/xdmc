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


#include "walker.h"
#include "simulation.h"

int walker :: count = 0;

double walker :: potential()
{
	// Evaluate the potential of the system
	// in the configuration described by this walker
	double ret = 0;
	for (int i = 0; i < particles.size(); ++i)
	{
		// Sum up external potential contributions
		for (int j=0; j<simulation.potentials.size(); ++j)
			ret += simulation.potentials[j]->potential(particles[i]);

		// Particle-particle interactions
		// note j<i => no double counting
		for (int j=0; j<i; ++j)
			ret += particles[i]->interaction(particles[j]);
	}
	return ret;
}

void walker :: diffuse(double tau)
{
	// Diffuse all of the particles (classical
	// particles will automatically not diffuse)
	for (int i=0; i<particles.size(); ++i)
		particles[i]->diffuse(tau);
}

walker :: walker()
{
	// Setup the walker with the particles
	// that describe the system.
	++ count;
	this->particles.clear();
	for (int i=0; i<simulation.template_system.size(); ++i)
		this->particles.push_back(simulation.template_system[i]->copy());
}

walker :: walker(std::vector<particle*> particles)
{
	++ count;
	this->particles = particles;
}

walker :: ~walker()
{
	// Clear up memory (delete all the particles).
	-- count;
	for (int i=0; i<particles.size(); ++i)
		delete(particles[i]);
}


walker* walker :: copy()
{
	// Return a copy of this walker (copy each of the particles)
	std::vector<particle*> copied_particles;
	for (int i=0; i<particles.size(); ++i)
		copied_particles.push_back(particles[i]->copy());
	return new walker(copied_particles);
}

void walker :: sample_wavefunction()
{
	for (int i=0; i<particles.size(); ++i)
		particles[i]->sample_wavefunction();
	simulation.wavefunction_file << "\n";
}



