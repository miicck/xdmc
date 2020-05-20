/*

    XDMC
    Copyright (C) Michael Hutcheon (email mjh261@cam.ac.uk)

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

#ifndef __CONSTANTS__
#define __CONSTANTS__

// Physical constants (we work in atomic units)
const double PI  = 3.141592653589793;      // You should know this one
const double AMU = 1822.88486193;          // Unified atomic mass unit / electron mass
const double PROTON_MASS  = 1.007276 *AMU; // Proton mass / electron mass
const double NEUTRON_MASS = 1.0086654*AMU; // Neutron mass / electron mass

// Finite difference amounts
const double EPS_X = 0.01;  // Epsilon for distances

// Error codes
const int MPI_ERROR = 1;    // Thrown when an MPI call fails

#endif

