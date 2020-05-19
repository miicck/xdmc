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

#define CATCH_CONFIG_RUNNER
#include "catch.h"
#include "params.h"
#include "particle.h"
#include "walker_collection.h"
#include "random.h"
#include "constants.h"

#include <iostream>

// Run the DMC calculation
void run_dmc()
{
    // Our DMC walkers
    params::progress_file << "Initializing walkers\n";
    walker_collection* walkers = new walker_collection();
    
    // Run our DMC iterations
    params::progress_file << "Starting DMC simulation\n";
    params::progress_file << "    Total setup time: " << params::time() << "s\n";
    params::dmc_start_clock = clock();

    for (params::dmc_iteration = 1;
         params::dmc_iteration <= params::dmc_iterations;
         params::dmc_iteration ++)
    {
        // Apply propagation of walkers
        walker_collection* walkers_last = walkers->copy();
        bool revert = !walkers->propagate(walkers_last);

        if (revert)
        {
            // Revert this iteration
            delete walkers;
            walkers = walkers_last;
        }
        else
            // Keep the new walkers
            delete walkers_last;

        // Estimate the new value for tau_nodes
        if (params::tau_nodes_estimator == "min_sep")
            params::tau_nodes = walkers->tau_nodes_min_sep();
        else if (params::tau_nodes_estimator != "none")
            throw "Unkown tau_nodes estimator!";

        walkers->write_output(revert);
    }

    // Output success message
    params::progress_file << "\nDone, total time: " << params::time() << "s.\n";
    
    // Free memory
    delete walkers;
}

// Check if arg is requesting help
bool is_help_arg(std::string arg)
{
    if (arg == "-h") return true;
    if (arg == "--h") return true;
    if (arg == "-help") return true;
    if (arg == "--help") return true;
    return false;
}

// Check if arg is requesting tests
bool is_test_arg(std::string arg)
{
    if (arg == "-t") return true;
    if (arg == "--t") return true;
    if (arg == "-test") return true;
    if (arg == "--test") return true;
    return false;
}

// Program entrypoint
int main(int argc, char** argv)
{
    // Check command line arguments for run mode
    for (int i=0; i<argc; ++i)
    {
        std::string arg = argv[i];
        if (is_help_arg(arg))
            params::run_mode = "help";
        else if (is_test_arg(arg))
            params::run_mode = "test";
    }

    // Initialize MPI etc.
    params::initialize();

    // Print usage info
    if (params::run_mode == "help")
        params::print_usage_info();

    // Run tests
    else if (params::run_mode == "test")
        for (int i=0; i<argc; ++i)
        {
            std::string arg = argv[i];
            if (is_test_arg(arg))
            {
                // Run tests, passing all other command
                // line argumenst to test session.
                char** new_argv = new char*[argc-1];
                for (int j=0; j<argc-1; ++j)
                    new_argv[j] = j < i ? argv[j] : argv[j+1];
                Catch::Session().run(argc-1, new_argv);
                delete[] new_argv;
            }
        }

    // Run the DMC simulation
    else if (params::load(argc, argv)) run_dmc();

    // Free memory used in the simulation specification
    params::free_memory();
}
