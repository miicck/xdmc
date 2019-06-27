import matplotlib.pyplot as plt
import numpy as np
import parser
import sys

start = int(sys.argv[1])
end   = int(sys.argv[2])

wfn = parser.parse_wavefunction(start, end)

weights = wfn[0]
coords  = wfn[1:]

num_particles = len(coords)
num_coords    = len(coords[0][0])
print "Particles  : ",num_particles
print "Coordinates: ",num_coords

for ip, particle in enumerate(coords):

	positive_walkers = np.array([p for w, p in zip(weights, particle) if w > 0]).T
	negative_walkers = np.array([p for w, p in zip(weights, particle) if w < 0]).T

	bins = len(positive_walkers[0])/10
	if bins > 1000: bins = 1000
	for ix, x in enumerate(positive_walkers):
		plt.subplot(num_particles, num_coords, ip*num_coords+ix+1)
		plt.hist(x, label="Positive", alpha=0.5, bins=bins)

	for ix, x in enumerate(negative_walkers):
		plt.subplot(num_particles, num_coords, ip*num_coords+ix+1)
		plt.hist(x, label="Negative", alpha=0.5, bins=bins)

	for ix in range(0,len(negative_walkers)):
		plt.subplot(num_particles, num_coords, ip*num_coords+ix+1)
		plt.xlabel("Particle {0} coord {1}".format(ip+1, ix+1))
		plt.legend()

plt.show()
