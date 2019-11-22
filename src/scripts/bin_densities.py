from parser import call_on_walkers
import numpy as np
import sys

# Read the number of spatial dimensions 
with open("progress") as f:
    for l in f:
        if "dimensions" in l:
            dims = int(l.split()[-1])
            break

print("Dimensions: {0}".format(dims))

BINS         = 100
SCALE        = 4.0
shape        = [BINS for i in range(0, dims)]
density      = np.zeros(shape)
cond_density = np.zeros(shape)

def ifc(x):
    # Get the bin index from a coordinate
    return int(BINS*(SCALE+x)/(2*SCALE))

def coords(x):
    # Vectro version of ifc
    return [ifc(xi) for xi in x]

# Read the spins of the particles
index_spins = {}
index_particles = []
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
            spin  = l.split(":")[-1].replace("]\n","").strip()
            index_spins[index] = spin
            index_particles.append(l.replace("\n",""))
        except:
            break

for p in index_particles:
    print ("particle "+p)

plot_indicies   = [int(i) for i in input("Particle(s) to accumulate density for: ").split()]
ignore_indicies = []
fixed_positions = {}

if len(plot_indicies) < len(index_particles):

    ignore_indicies = [int(i) for i in input("Particle(s) to ignore: ").split()]

    if len(plot_indicies) + len(ignore_indicies) < len(index_particles):
        for j in range(0, len(index_particles)):
            if not (j in plot_indicies):
                if not (j in ignore_indicies):
                    fs = "Fixed position for particle {0}: ".format(j)
                    fixed_positions[j] = [float(x) for x in input(fs).split()]

print("Plotting particles: ", plot_indicies)
print("Ignoring particles: ", ignore_indicies)
print("Fixing particles  : ")
for i in fixed_positions:
    print("    {0} : {1}".format(i, fixed_positions[i]))

def in_range(c):
    # Returns true if the coordinates c are
    # in binning range
    for ci in c:
        if (ci < 0) or (ci >= BINS):
            return False
    return True

def to_call(i, w, x):
    # This will be called on every walker
    global density, cond_density
    global plot_indicies, ignore_indicies, fixed_positions

    # Evaluate if fixed particles are at the
    # right place
    condition = True
    for i in fixed_positions:
        if coords(x[i]) != coords(fixed_positions[i]):
            condition = False
            break

    # Sum up particle densities
    for i in plot_indicies:
        p = coords(x[i])

        # Check coordinates are in bin range
        if not in_range(p): continue

        # Accumulate the weight
        density[tuple(p)] += w
        if i > 0 and condition:
            cond_density[tuple(p)] += w
            
start = int(input("start iteration: "))
end   = int(input("end iteration  : "))
call_on_walkers(to_call, skip_iter=lambda i: i <= start or i >= end)

def normalized(d, norm):
    n = np.array(d)
    while True:
        try:
            n = np.trapz(n)
        except:
            break
    if abs(n) < 10e-6: return d
    return norm*d/n

with open("wavefunction_0") as f:
    for line in f:
        if "#" in line: continue
        particle_count = len(line.split(";"))
        break

density      = normalized(density**2,      particle_count)
cond_density = normalized(cond_density**2, particle_count-1)

np.save("density.binned",      density)
np.save("cond_density.binned", cond_density)
