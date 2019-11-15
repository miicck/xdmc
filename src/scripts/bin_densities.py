from parser import call_on_walkers
import numpy as np
import sys

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

def on_plane(coords):
    # Returns true if the coordinates 
    # are on the x-y plane
    if len(coords) < 2:
        return True
    for c in coords[2:]:
        if c != ifc(0):
            return False
    return True

fixed_coords = [float(x) for x in input("Eneter the fixed particle coords: ").split()]
print("Fixed coords", fixed_coords)

def to_call(i, w, x):
    # This will be called on every walker
    global density, cond_density, fixed_coords

    # Evaluate if condition is true
    condition = coords(x[0]) == coords(fixed_coords)

    # Sum up particle densities
    for i, p in enumerate(x):
        p = coords(p)
        if not on_plane(p):
            continue

        try:
            density[p] += 1.0
            if i > 0 and condition:
                cond_density[p] += 1.0
        except IndexError:
            continue
            
start = int(sys.argv[1])
end   = int(sys.argv[2])
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

density = normalized(density, particle_count)
cond_density = normalized(cond_density, particle_count-1)

np.save("density.binned", density)
np.save("cond_density.binned", cond_density)
