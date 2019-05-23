#include "simulation.h"

simulation_spec simulation;

void simulation_spec::open_output_files()
{
	iterations_file.open(  "iterations_"   + std::to_string(pid));
	wavefunction_file.open("wavefunction_" + std::to_string(pid));
}

void simulation_spec::free_memory()
{
	iterations_file.close();
	wavefunction_file.close();

	for (int i=0; i<template_system.size(); ++i)
		delete template_system[i];

	for (int i=0; i<potentials.size(); ++i)
		delete potentials[i];
}
