# 
#     XDMC
#     Copyright (C) Michael Hutcheon (email mjh261@cam.ac.uk)
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
from parser import parse_wavefunction
import numpy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import sys

# Get input variables
j_fixed = int(input("Please enter the fixed particle index (particle indexing starts at 1): "))
fixed_position = [float(w) for w in input("Please enter it's fixed location: ").split()]

# Read the spins of the particles
index_spins = {}
with open("progress") as f:
    
    started = False
    for l in f:
        if "Particles" in l:
            started = True
            continue
        if not started:
            continue

        try:
            index = int(l.split(":")[0])
            spin  = l.split(":")[-1].replace("]\n","")
            index_spins[index+1] = spin
        except:
            break

# Work out the particles that contribute to the
# conditional density
other_particles = []
for j in index_spins:
    if j == j_fixed: continue
    if index_spins[j] != index_spins[j_fixed]: continue
    other_particles.append(j)
all_particles = list(other_particles)
all_particles.append(j_fixed)

# Print the particles we evaluate the condional density for
pinf  = "Calculating conditional density for particles: "
pinf += ", ".join([str(p) for p in other_particles])
print(pinf)

# Read the wavefunction
wfn = parse_wavefunction(sys.argv[1:])

# Bins for accumulated things
BINS         = 20
SCALE        = 4.0
density      = np.zeros((BINS, BINS))
cond_density = np.zeros((BINS, BINS))

# Convert x to a bin coordinate
def ifc(x):
    return int(BINS*(x+SCALE)/(2*SCALE))

# Array version of above
def coords(x):
    return [ifc(xi) for xi in x]

# How much a particle at x counts as being at y
def att(x,y):

    # Have to be in same bin
    for i in range(0, len(x)):
        if ifc(x[i]) != ifc(y[i]): return 0.0
    return 1.0

fixed_position = coords(fixed_position)

# i indexes walkers
for i in range(0, len(wfn[0])):
    
    # Weight of the walker
    w = wfn[0][i]

    # Check if the fixed particle is at the fixed position
    xjf = coords(wfn[j_fixed][i])
    condition_satisfied = (xjf == fixed_position)

    # Accumulate densities
    for j in other_particles:
        xj = coords(wfn[j][i])

        # Has to be at same z-coordinate
        if xj[2] != fixed_position[2]:
            continue

        try:
            density[xj[0]][xj[1]] += 1.0
            if condition_satisfied:
                cond_density[xj[0]][xj[1]] += 1.0
        except IndexError:
            continue

# Normalize a density to n particles
def normalize(den, n):
    integral =  np.trapz(np.trapz(den))
    if abs(integral) < 10e-4: return den
    ret = den * n / integral
    print(n, np.trapz(np.trapz(ret)))
    return ret

# Plot a density
def plot_density(den):
    plt.imshow(den, extent=(-SCALE,SCALE,-SCALE,SCALE))
    plt.axhline(0, color="black")
    plt.axvline(0, color="black")

plt.subplot(221)
plot_density(density)
plt.subplot(222)
plot_density(cond_density)

# Normalize the densities
density      = normalize(density,      len(all_particles)  )
cond_density = normalize(cond_density, len(all_particles)-1)

nxc = cond_density - density
plt.subplot(223)
plot_density(nxc)

plt.show()

