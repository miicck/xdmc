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
	{
		particles[i]->sample_wavefunction();
		if (i < particles.size() - 1)
			simulation.wavefunction_file << ";";
	}
	simulation.wavefunction_file << "\n";
}

