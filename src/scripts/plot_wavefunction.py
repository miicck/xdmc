import os
import numpy as np
import matplotlib.pyplot as plt

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

# Histogram each coordinate of each particle
for pi in range(0, particle_count):
	for ci in range(0, coord_count):

		psi = wavefunction[pi][ci]
		sp_index = 1 + pi + ci*particle_count
		plt.subplot(particle_count, coord_count, sp_index)
		plt.hist(psi, bins=len(psi)/100)
		plt.xlabel("Particle "+str(pi)+" coord "+str(ci))

plt.show()
