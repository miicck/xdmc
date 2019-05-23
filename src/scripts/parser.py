import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def parse_evolution():
	# Labels for the plots
	y_axes = ["Population", "Trial energy (Hartree)"]
	how_to_combine = [
		lambda(p) : np.sum(p, axis=0), # Sum populations across processes
		lambda(e) : np.mean(e, axis=0) # Average energies across processes
	]

	# Read in the evolution from each process
	all_data = []
	for f in os.listdir("."):
		if not f.startswith("evolution"): continue

		# Read in our data from the evolution file
		data = []
		for line in open(f).read().split("\n"):
			if len(line.strip()) > 0:
				data.append([float(i) for i in line.split(",")])
		data = zip(*data)
		all_data.append(data)

	# Combine the data from each process
	data = []
	for i in range(0, len(all_data[0])):
		data.append(how_to_combine[i]([d[i] for d in all_data]))
	return np.array([y_axes, data])

def parse_wavefunction():

	# Read the number of particles and coordinates
	# from  the first wavefunction file
	f0 = open("wavefunction_0")
	lines = f0.read().split("\n")
	f0.close()
	particle_count = len(lines[0].split(";")) - 1
	coord_count    = len(lines[0].split(";")[0].split(","))

	# Prepare the wavefunction object
	wavefunction = []
	for i in range(0, particle_count): 
		wavefunction.append([])
		for j in range(0, coord_count): 
			wavefunction[i].append([])

	# Read the wavefunction samples from each file
	for f in os.listdir("."):
		if not f.startswith("wavefunction"): continue
		for line in open(f).read().split("\n")[0:-1]:
			for pi, p in enumerate(line.split(";")[0:-1]):
				for ci, x in enumerate([float(i) for i in p.split(",")]):
					wavefunction[pi][ci].append(x)
	
	return np.array(wavefunction)
