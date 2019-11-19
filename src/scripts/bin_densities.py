from parser import call_on_walkers
import numpy as np
import sys

BINS         = 40
SCALE        = 4.0
density      = np.zeros((BINS,BINS))
cond_density = np.zeros((BINS,BINS))

def ifc(x):
    # Get the bin index from a coordinate
    return int(BINS*(SCALE+x)/(2*SCALE))

def coords(x):
    # Vectro version of ifc
    return [ifc(xi) for xi in x]

def on_plane(coords):
    # Returns true if the coordinates 
    # are on the x-y plane
    for c in coords[2:]:
        if c != ifc(0):
            return False
    return True

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
            index_spins[index] = spin
        except:
            break

# Same spin particles
same_spin = [i for i in index_spins if index_spins[i] == index_spins[0]] 
print("Particle indicies: ", same_spin)

def to_call(i, w, x):
    # This will be called on every walker
    global density, cond_density, same_spin

    # Evaluate if condition is true
    condition = coords(x[0]) == coords([1.0, 0, 0])

    # Sum up particle densities
    for i in same_spin:
        p = coords(x[i])
        if not on_plane(p):
            continue

        try:
            density[p[0],p[1]] += 1.0
            if condition and (i > 0):
                cond_density[p[0], p[1]] += 1.0
        except IndexError:
            continue
            
start = int(sys.argv[1])
end   = int(sys.argv[2])
call_on_walkers(to_call, skip_iter=lambda i: i <= start or i >= end)

def normalized(d, norm):
    n = np.trapz(np.trapz(d))
    if abs(n) < 10e-6: return d
    return norm*d/n

with open("wavefunction_0") as f:
    for line in f:
        if "#" in line: continue
        particle_count = len(line.split(";"))
        break

density = normalized(density, particle_count)
cond_density = normalized(cond_density, particle_count-1)

with open("density.binned", "w") as f:
    for x in density:   
        f.write(" ".join([str(xi) for xi in x])+"\n")
    f.close()

with open("cond_density.binned", "w") as f:
    for x in cond_density:   
        f.write(" ".join([str(xi) for xi in x])+"\n")
    f.close()
