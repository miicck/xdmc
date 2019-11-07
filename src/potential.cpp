/*

    XDMC
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

#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include "dmc_math.h"
#include "params.h"
#include "potential.h"
#include "constants.h"

std::string grid_potential :: one_line_description()
{
    std::stringstream des;
    des << "Grid potential (extent = " << extent << " grid size = " << grid_size << ")";
    return des.str();
}

grid_potential :: grid_potential(std::string filename)
{
    std::ifstream file(filename);

    // Read the number of dimensions, the grid size and the extent from file
    unsigned read_dimensions;
    file >> read_dimensions >> this->grid_size >> this->extent;

    // Check dimensionality is correct
    if (read_dimensions != params::dimensions)
    {
        params::error_file << "Incorrect dimensionality in grid_potential!\n";
        throw "Incorrect dimensionality in grid_potential!";
    }

    // Read data from file
    this->data = new double[read_dimensions * this->grid_size];
    double val;
    unsigned n = 0;
    while(file >> val)
    {
        this->data[n] = val;
        ++n;
        if (n >= read_dimensions * this->grid_size)
        {
            params::error_file 
                << "Warning: grid_potential file has more lines than expected "
                << "please check it is correct.";
            break;
        }
    }
}

double grid_potential :: potential(particle* p)
{
    int coord = 0;
    int stride = 1;

    for (unsigned i=0; i<params::dimensions; ++i)
    {
        // Work out the coordinate I'm at
        double frac = (p->coords[i] + extent)/(2*extent);
        int c = int(grid_size*frac);
        
        // Outside of grid => infinite potential
        if (c < 0 || c >= grid_size) 
            return INFINITY;

        // Add the resulting stride to the coodinate
        coord  += c*stride;
        stride *= grid_size;
    }

    return data[coord];
}

std::string harmonic_well :: one_line_description()
{
    std::stringstream des;
    des << "Harmonic well (omega = " << omega << ")";
    return des.str();
}

double harmonic_well::potential(particle* p)
{
    double r2 = 0;
    for (unsigned i=0; i<params::dimensions; ++i)
        r2 += p->coords[i]*p->coords[i];
    return 0.5*r2*omega*omega;
}

std::string atomic_potential :: one_line_description()
{
    std::stringstream des;
    des << "Atomic potential (charge = " << charge << ") position: ";
    for (unsigned i=0; i<params::dimensions; ++i)
        des << coords[i] << " ";
    return des.str();
}

double atomic_potential :: potential(particle* p)
{
    double r = 0;
    for (unsigned i=0; i<params::dimensions; ++i)
    {
        double dxi = p->coords[i] - this->coords[i];
        r += dxi * dxi;
    }
    r = sqrt(r);
    return coulomb(this->charge, p->charge, r);
}










