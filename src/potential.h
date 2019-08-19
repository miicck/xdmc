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

#ifndef __POTENTIAL__
#define __POTENTIAL__

#include <string>
#include "particle.h"

class external_potential
{
public:
    virtual double potential(particle* p)=0;
    virtual std::string one_line_description()=0;
    virtual ~external_potential() { }
};

class grid_potential : public external_potential
{
public:
    grid_potential(std::string filename);
    virtual double potential(particle* p);
    virtual std::string one_line_description();
    virtual ~grid_potential() { delete[] this->data; }
private:
    int grid_size;
    double extent;
    double* data;
};

class harmonic_well : public external_potential
{
public:
    harmonic_well(double omega) { this->omega = omega; }
    virtual double potential(particle* p);
    virtual std::string one_line_description();
private:
    double omega = 1;
};

class atomic_potential : public external_potential
{
public:
    atomic_potential(double charge, double* coords)
    {
        this->coords = coords;
        this->charge = charge;
    }
    
    virtual ~atomic_potential()
    {
        delete[] this->coords;
    }

    virtual double potential(particle* p);
    virtual std::string one_line_description();
private:
    double  charge;
    double* coords;
};

#endif


