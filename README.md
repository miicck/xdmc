This program is a simple DMC implementation that reads a list of atoms and DMC parameters from the "input" file and outputs the potential energy at each timestep to the "out" file. The positions of the electrons are also sampled to the "electrons" file.

# Scripts
## build.sh
Builds the main.cpp into the dmc executable
## run.sh
Runs the dmc executable using mpirun then calls plot.py on the results.
## plot.py
Plots the population and potential against iteration count.
## fit_hydrogen_electrons.py
Fits the sampled electron wavefunction contained in the "electrons" file to a hydrogen-1s-like wavefunction and compares with the exact hydrogen 1s wavefunction.
