# 
#     XDMC
#     Copyright (C) 2019 Michael Hutcheon (email mjh261@cam.ac.uk)
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     For a copy of the GNU General Public License see <https://www.gnu.org/licenses/>.
# 
import matplotlib.pyplot as plt
import numpy as np
import parser
import sys

wfn = parser.parse_wavefunction(sys.argv[1:])

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
                        try:
                            index = int(bins * (x - min_x)/(max_x - min_x) - 0.5)
                        except ValueError as ex:
                            print(ex)
                        vals[index] += weights[ix]

                vals /= np.sqrt(np.trapz(vals**2, x=xs))
                plt.subplot(num_particles, num_coords, 1 + ic + ip*num_coords)
                plt.xlabel("Particle {0} coordinate {1}".format(ip+1, ic+1))
                plt.ylabel("Normalized wavefunction")
                plt.plot(xs, vals)
plt.show()

