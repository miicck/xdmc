#include <mpi.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>

#include "parameters.h"

// The collection of parameters
parameter_collection params;

void parameter_collection::read_input()
{
    std::ifstream input("input");
    for (std::string line; getline(input, line); )
    {
        
    }
}

void parameter_collection::output_sim_details(){}

void test()
{
   params.set<int>("walkers", 69);
   std::cout << params.get<int>("walkers") << "\n";
}

void parameter_collection :: load(int argc, char** argv)
{
    params["walkers"] = new parameter<int>
    (
        // Default value
        1000,  

        // Description
        "The target number of DMC walkers. The actual number of DMC walkers "
        "will fluctuate during propagation, but will be bias towards this value."
    );

    test();
}

/*
void parameter_collection::load(int argc, char** argv)
{
    // Get the start time so we can time stuff
    start_clock = clock();

    // Initialize mpi
    if (MPI_Init(&argc, &argv) != 0) exit(MPI_ERROR);

    // Get the number of processes and my id within them
    if (MPI_Comm_size(MPI_COMM_WORLD, &np.value)  != 0) exit(MPI_ERROR);
    if (MPI_Comm_rank(MPI_COMM_WORLD, &pid.value) != 0) exit(MPI_ERROR);

    // Seed random number generator
    srand(pid*clock());

    // Read our input and setup parameters accordingly 
    // do for each process sequentially to avoid access issues
    for (int pid_read = 0; pid_read < np; ++ pid_read)
    {
            if (pid == pid_read) read_input();
            MPI_Barrier(MPI_COMM_WORLD);
    }

    // Open various output files
    if (pid == 0)
    {
        // Files on the root process
        progress_file.open("progress");
        evolution_file.open("evolution");
    }

    // Files on all processes have their pid appended
    error_file.open("error_"+std::to_string(pid.value));
    error_file.auto_flush = true;
    wavefunction_file.open("wavefunction_"+std::to_string(pid.value));

    // Output parameters to the progress file
    output_sim_details();
}
*/
