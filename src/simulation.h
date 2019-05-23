#ifndef __SIMULATION__
#define __SIMULATION__

#include <vector>
#include <fstream>

#include "particle.h"
#include "potential.h"

struct simulation_spec
{
        // Program parameters
        int pid               = 0;      // The MPI PID of this process
        int np                = 1;      // The number of MPI processes
        int dimensions        = 3;      // The dimensionality of our system
        int target_population = 500;    // The number of DMC walkers per process
        int dmc_iterations    = 1000;   // The number of DMC iterations to carry out
        double tau            = 0.01;   // The DMC timestep
        double trial_energy   = 0;      // Energy used to control the DMC population

        // The system which will be copied to generate walkers
        std::vector<particle*> template_system;

        // The potentials applied to the system (additive)
        std::vector<external_potential*> potentials;

        // Output files
        std::ofstream wavefunction_file;
        std::ofstream iterations_file;

        void open_output_files(); // Opens output files
        void free_memory();       // Closes output files and frees template_system and potentials
};

extern simulation_spec simulation;

#endif
