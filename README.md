<h1>XDMC</h1>
XDMC is a code that implements eXchange-Diffusion Monte Carlo, as described here: <br/>
https://link.aps.org/doi/10.1103/PhysRevE.102.042105 <br/>
An open access version is available here: <br/>
https://arxiv.org/abs/2001.02215 <br/>

<h2>Compiling</h2>
To build the executable, run build.py (which can be found in the root directory).
This will generate code, compile and link an executable (called xdmc)
in the root directory. This will detect the number of cores on your system and
use them all to speed up compilation. Unit tests can be ran (in a few seconds), either in
serial or parallel, in the following manner:

        $ xdmc -t # serial
        $ mpirun xdmc -t # parallel

<h2>Example usage</h2>
The xdmc executable requires a single input file, called simply "input". The program must
be executed in the same directory as this file. It can be invoked in serial, or in parallel with MPI:

        $ xdmc # serial
        $ mpirun xdmc # parallel with MPI
        
The input file contains a list of keywords, information about which can be obtained by invoking xdmc with the -h option

        $ xdmc -h
        ... info about input parameters

The input file contains the algorithm settings (number of walkers, timestep, iterations to carry out etc.)
as well as the description of the system to simulate. The system to simulate consists (optionally) of the
specification of an external potential and the specification of (at least one) quantum particle, the wavefunction
of which will be sampled. An example input file for a lithium atom is shown below:

        # Lines starting with any of #!// will be treated as comments
        
        # Algorithm settings
        dimensions          3
        walkers             10000
        iterations          10000
        tau                 0.001
        tau_nodes           1.0
        diffusion_scheme    stochastic_nodes_mpi
        pre_diffusion       0.5
        coulomb_softening   0.00001
        max_weight          4

        # Atomic potential with charge 3 at x=0, y=0, z=0
        atomic_potential 3 0 0 0

        # Quantum particles to treat with monte carlo, described by
        # a mass, charge and spin (the name is just a label).
        #        name      mass  charge  spin  Initial position
        particle electron  1     -1       1    0 0 0
        particle electron  1     -1       1    0 0 0
        particle electron  1     -1      -1    0 0 0
        
Running this input file will produce a variety of output files, listed below. Some types of output will be distributed to different files for each process. These have the PID of the process appended (e.g wavefunction_0 is the wavefunction file for the root process). <br>
- **progress** File updated with a high-level report of the progress of the calculation (human readable). <br>
- **evolution** File containing the evaluation of expectation values for each DMC iteration. <br>
- **wavefunction_n** File containing all of the walker configurations for each iteration (large). <br>
- **nodal_surface_n** File containing the configurations of walkers that were killed due to crossing a nodal surface (not written by default). <br>

<h3>Analysis</h3>
Various analysis scripts can be found in the /src/scripts directory. The simplest of these is the plot_evolution.py script, which will plot the expectation values calculated (including the energy/walker population) vs. DMC timestep.
