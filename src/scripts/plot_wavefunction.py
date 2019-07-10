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

print("Particles  : ",num_particles)
print("Coordinates: ",num_coords)

for ip, p in enumerate(coords):
        # ip = particle index
        # p  = sampled particle locations

        p = np.array(p).T
        for ic, c in enumerate(p):
                # ic = coordinate index
                # c  = coordinate samples

                min_x = min(c)
                max_x = max(c)
                bins  = min(1000, int(len(c)/10))
                vals  = np.zeros(bins)
                xs    = np.linspace(min_x, max_x, bins)
                for ix, x in enumerate(c):
                        # ix = sample index
                        # x  = sample
                        index = int(bins * (x - min_x)/(max_x - min_x) - 0.5)
                        vals[index] += weights[ix]

                vals /= np.sqrt(np.trapz(vals**2, x=xs))
                plt.subplot(num_particles, num_coords, 1 + ic + ip*num_coords)
                plt.xlabel("Particle {0} coordinate {1}".format(ip+1, ic+1))
                plt.ylabel("Normalized wavefunction")
                plt.plot(xs, vals)
plt.show()
