import numpy as np
import matplotlib.pyplot as plt

lines = open("wavefunction").read().split("\n")
particle_count = len(lines[0].split(";"))
coord_count    = len(lines[0].split(";")[0].split(","))

wavefunction = []
for i in range(0, particle_count): 
	wavefunction.append([])
	for j in range(0, coord_count): 
		wavefunction[i].append([])

for line in lines:
	if len(line) == 0: continue
	for pi, p in enumerate(line.split(";")):
		for ci, x in enumerate([float(i) for i in p.split(",")]):
			wavefunction[pi][ci].append(x)

for pi in range(0, particle_count):
	for ci in range(0, coord_count):

		psi = wavefunction[pi][ci]
		sp_index = 1 + pi + ci*particle_count
		plt.subplot(particle_count, coord_count, sp_index)
		plt.hist(psi, bins=len(psi)/100)
		plt.xlabel("Particle "+str(pi)+" coord "+str(ci))

plt.show()
