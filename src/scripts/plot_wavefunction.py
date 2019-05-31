import matplotlib.pyplot as plt
import numpy as np
import parser
import sys

start = int(sys.argv[1])
end   = int(sys.argv[2])

wavefunction = parser.parse_wavefunction(start, end)
wfn = parser.transpose_wavefunction(wavefunction)

weights = wfn[0]
coords  = wfn[1:]

num_particles = len(coords)
num_coords    = len(coords[0][0])
print "Particles  : ",num_particles
print "Coordinates: ",num_coords

for ip, particle in enumerate(coords):

	positive_walkers = np.array([p for w, p in zip(weights, particle) if w > 0])
	negative_walkers = np.array([p for w, p in zip(weights, particle) if w < 0])

	for ix, x in enumerate(positive_walkers.T):
		plt.subplot(num_particles, num_coords, ip*num_coords+ix+1)
		plt.hist(x, label="Positive", alpha=0.5, bins=len(x)/10)

	for ix, x in enumerate(negative_walkers.T):
		plt.subplot(num_particles, num_coords, ip*num_coords+ix+1)
		plt.hist(x, label="Negative", alpha=0.5, bins=len(x)/10)
		plt.xlabel("Particle {0} coord {1}".format(ip+1, ix+1))
		plt.legend()

plt.show()
