This program is a simple DMC implementation that reads a list of atoms and DMC parameters from the "input" file and outputs the potential energy at each timestep to the "out" file. The positions of the electrons are also sampled to the "electrons" file.

# Scripts
None of these scripts require any arguments.
## build.sh
Builds the main.cpp into the dmc executable. Only tested with v7.3.0 of the mpiCC compiler (an mpi wrapper of the gcc compiler).
## run.sh
Runs the dmc executable using mpirun then calls plot.py on the results. Only tested with mpirun v2.1.1 (open MPI).
## plot.py
Plots the population and potential against iteration count. Only tested with python 2.7.15rc1.
## fit_hydrogen_electrons.py
Fits the sampled electron wavefunction contained in the "electrons" file to a hydrogen-1s-like wavefunction and compares with the exact hydrogen 1s wavefunction. Only tested with python 2.7.15rc1.
