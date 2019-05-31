import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def parse_evolution():

	# Labels for the plots
	y_axes = [
		"Population",
		"Trial energy (Hartree)",
		"Average weight",
		"Average weight^2"
	]

	# How to combine axes across processes
	how_to_combine = [
		lambda(p) : np.sum(p, axis=0) , # Sum populations across processes
		lambda(e) : np.mean(e, axis=0), # Average energies across processes
		lambda(e) : np.mean(e, axis=0), # Average energies across processes
		lambda(w) : np.mean(w, axis=0), # Average weights across processes
		lambda(w) : np.mean(w, axis=0), # Average weights across processes
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

def transpose_wavefunction(wfn):
	# Take a wavefunction of the form
	# wavefunction[iteration number][walker number] = [weight, x0, x1...]
	# and return the array
	# [[w1, w2, w3 ... wN], [x0_1, x0_2, x0_3 ... x0_N], [x1_1, x1_2, x1_3 ... x1_N] ...]
	# where 1 -> N runs over all iterations and walker numbers 
	# (i.e walkers from all iterations are combined)
	wfn_combined = []
	for i in range(0, len(wfn)):
		wfn_combined.extend(wfn[i])
	return np.array(wfn_combined).T

def parse_wavefunction(iter_start=0, iter_end=None):
	
	# Combines wavefunctions across processes
	wfs = []
	for f in os.listdir("."):
		if f.startswith("wavefunction_"):
			wf = parse_wavefunction_file(f, iter_start, iter_end)
			wfs.append(wf)
			continue # Only do one for now
	wfn = wfs[0]
	for wf in wfs[1:]:
		for i in range(0, len(wfn)):
			wfn[i].extend(wf[i])
	return wfn
			

def parse_wavefunction_file(filename, iter_start=0, iter_end=None):
	
	# Parse wavefunction into the form
	# wavefunction[iteration number][walker number] = [weight, x0, x1...]
	# where x0, x1 ... are the coordinates of the particles
	f0 = open(filename)
	lines = f0.read().split("\n")
	f0.close()

	# Get the start location of each iteration
	i_iter_start = [0]
	for i, l in enumerate(lines):
		if not "#" in l: continue
		i_iter_start.append(i+1)
	
	# Validate the requested number of iterations
	if iter_end is None or iter_end > len(i_iter_start) - 1:
		iter_end = len(i_iter_start)-1

	wavefunction = []
	iteration    = []

	# Read the requested number of iterations
	for i in range(i_iter_start[iter_start],
		       i_iter_start[iter_end]):
		l = lines[i]
		if "#" in l:
			wavefunction.append(iteration)
			iteration = []
			continue

		weight, coord_data = l.split(":")
		particles = coord_data.split(";")[0:-1]

		dat = [float(weight)]
		for p in particles:
			x = [float(xi) for xi in p.split(",")]
			dat.append(x)

		iteration.append(dat)

	return wavefunction

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
			weight = float(line.split(":")[0])
			line = line.split(":")[1]
			for pi, p in enumerate(line.split(";")[0:-1]):
				for ci, x in enumerate([float(i) for i in p.split(",")]):
					wavefunction[pi][ci].append([weight, x])
	
	return np.array(wavefunction)
