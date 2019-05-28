/*
    DMCPP
    test
    Copyright (C) 2019 Michael Hutcheon

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
#include "potential.h"

double harmonic_well::potential(particle* p)
{
	double r2 = 0;
	for (int i=0; i<simulation.dimensions; ++i)
		r2 += p->coords[i]*p->coords[i];
	return 0.5*r2*omega*omega;
}
