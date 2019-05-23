import matplotlib.pyplot as plt
from parser import parse_wavefunction

wavefunction = parse_wavefunction()

particle_count = len(wavefunction)
coord_count    = len(wavefunction[0])

# Histogram each coordinate of each particle
for pi in range(0, particle_count):
	for ci in range(0, coord_count):

		psi = wavefunction[pi][ci]
		sp_index = 1 + pi + ci*particle_count
		plt.subplot(particle_count, coord_count, sp_index)
		plt.hist(psi, bins=len(psi)/100)
		plt.xlabel("Particle "+str(pi)+" coord "+str(ci))

plt.show()
