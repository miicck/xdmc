#include "simulation.h"
#include "potential.h"

double harmonic_well::potential(particle* p)
{
	double r2 = 0;
	for (int i=0; i<simulation.dimensions; ++i)
		r2 += p->coords[i]*p->coords[i];
	return 0.5*r2*omega*omega;
}
