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

#ifndef __SIMULATION__
#define __SIMULATION__

#include <vector>
#include <fstream>

#include "particle.h"
#include "potential.h"

// Class to help with I/O. Leaves the file
// closed until it's needed.
class output_file
{
public:
    output_file() { filename = "/dev/null"; }
    output_file(std::string fn) { filename = fn; }

    template<class T>
    output_file& operator<<(T t)
    {
        if (!file.is_open()) file.open(filename);
        file << t;
        if (auto_flush) flush();
        return (*this);
    }

    void open(std::string fn) { filename = fn; }
    void close() { if(file.is_open()) file.close(); }
    void flush() { if(file.is_open()) file.flush(); }

    // Set to true to automatically flush the file after each write
    bool auto_flush = false; 

private:
    std::string filename;
    std::ofstream file;
};

// This namespace represents a complete specification for
// a simulation to run and progress, including ready-to-write
// output files and control parameters. 
namespace params
{
    //%%%%%%%%%%%%%%%%%%%%//
    // PROGRAM PARAMETERS //
    //%%%%%%%%%%%%%%%%%%%%//

    extern int pid;                   // The MPI PID of this process
    extern int np;                    // The number of MPI processes
    extern int dimensions;            // The dimensionality of our system
    extern int target_population;     // The number of DMC walkers per process
    extern double max_pop_ratio;      // The max allowed population as fraction of target
    extern double min_pop_ratio;      // The min allowed population as fraction of target
    extern int dmc_iterations;        // The number of DMC iterations to carry out
    extern double tau;                // The DMC timestep
    extern double tau_c_ratio;        // The ratio for tau_c to tau
    extern double trial_energy;       // Energy used to control the DMC population
    extern double pre_diffusion;      // The amount that we diffuse walkers before run
    extern bool write_wavefunction;   // True if we should write wavefunction files
    extern bool exchange_moves;       // True if we should make exchange moves
    extern double exchange_prob;      // The probability of making an exchange move
    extern std::string cancel_scheme; // Cancellation scheme to use
    extern bool correct_seperations;  // True if walker seperation correction is carried out

    // The external potentials applied to the system (additive)
    extern std::vector<external_potential*> potentials;

    // The system which will be copied to generate walkers
    extern std::vector<particle*> template_system;

    // Exchange information about the system
    extern std::vector<int> exchange_pairs;
    extern std::vector<int> exchange_values;

    // Output files
    extern output_file wavefunction_file;
    extern output_file evolution_file;
    extern output_file progress_file;
    extern output_file error_file;

    //%%%%%%%%%%%//
    // FUNCTIONS //
    //%%%%%%%%%%%//

    // The total charge in the system
    double total_charge();

    // Flush output files so we have information if a run terminates
    void flush();

    // Output details of the simulation to the progress file
    void output_sim_details();

    // Loads system from input, initializes MPI, opens output files etc.
    void load(int argc, char** argv);

    // Closes output files and frees template_system and potentials
    void free_memory();

    // Get the time since startup
    double time();
    
    // Reads/parses the input file
    void read_input();

    // Parse a particle from an input line
    void parse_particle(std::vector<std::string> split);

    // Parse an atomic potential from an input file
    void parse_atomic_potential(std::vector<std::string> split);

};

#endif


