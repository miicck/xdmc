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

#ifndef __POTENTIAL__
#define __POTENTIAL__

#include "particle.h"

class external_potential
{
public:
        virtual double potential(particle* p)=0;
};

class harmonic_well : public external_potential
{
public:
        harmonic_well(double omega) { this->omega = omega; }
        virtual double potential(particle* p);
private:
        double omega = 1;
};

#endif
