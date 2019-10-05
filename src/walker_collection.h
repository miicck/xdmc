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

#ifndef __WALKER_COLLECTION__
#define __WALKER_COLLECTION__

#include "walker.h"

// A collection of walkers
class walker_collection
{
public:
    walker_collection();
    ~walker_collection();
    walker_collection* copy();

    bool propagate(walker_collection* walkers_last);
    void write_output(bool reverted);

    double average_weight();
    double average_mod_weight();
    double average_potential();
    double sum_mod_weight();

private:
    walker_collection(std::vector<walker*> walkers_in) : walkers(walkers_in) {}

    double diffused_wavefunction(walker* w, double tau);
    double* diffused_wavefunction_signed(walker* w, double tau);
    double* exchange_diffused_wfn_signed(walker* w, double tau);

    void make_exchange_moves();
    void make_diffusive_moves(walker_collection* walkers_last);
    void diffuse_exact_1d();
    void diffuse_max_seperation(walker_collection* walkers_last);
    void diffuse_stochastic_nodes(walker_collection* walkers_last);
    void exchange_diffuse(walker_collection* walkers_last);
    void apply_renormalization();
    void branch();

    std::vector<walker*> walkers;
};

#endif
