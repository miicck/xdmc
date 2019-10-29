import numpy as np
import sys
import os

# The fixed coordinates that we plot the wavefunction for
read_fixed = {}
def fixed(i):
    global read_fixed
    if not i in read_fixed:
        s = "Please enter the fixed location of particle {0}: ".format(i)
        read_fixed[i] = [float(x) for x in input(s).split()]
    return read_fixed[i]

# The attenuation applied to weights that
# aren't quite at the fixed coordinates that we want
att_scale = float(input("Please enter the attenuation scale: "))
def attenuation(x,y):
    r2 = 0
    for i in range(0,len(x)):
        r2 += (x[i] - y[i])**2
    r2 /= att_scale**2
    return np.exp(-r2)

# Count the number of particles
with open("wavefunction_0") as f:
    for line in f:
        if "#" in line: continue
        particle_count = len(line.split(";"))
        break

# Get the fixed positions for the other particles
for i in range(1, particle_count):
    fixed(i)

# Some parameters
start_iter = int(input("Please enter the start iteration: "))
end_iter   = int(input("Please enter the end iteration  : "))
RES        = 100
SCALE      = 4.0
bins       = np.zeros((RES, RES))

# Loop over wavefunction files
datapoints     = 0
eff_datapoints = 0
i_wfn          = -1
while True:
    i_wfn += 1

    filename = "wavefunction_{0}".format(i_wfn)
    if not os.path.isfile(filename):
        break

    print("Reading "+filename)
    with open(filename) as f:

        iteration = 0

        # Loop over lines in wfn file
        for line in f:

            # We've reached a new iteration
            if line.startswith("#"):
                iteration += 1
                if iteration % 1000 == 0:
                    print("    Iteration {0}".format(iteration))
                continue

            if iteration < start_iter: continue
            if iteration > end_iter  : break

            # Parse the coordinates from this line
            w, coords = line.split(":")
            w = float(w)
            xs = [[float(x) for x in c.split(",")] for c in coords.split(";")]

            # The attenuation
            att = 1

            # Fix other coordinates of first particle to 0
            for j in range(2, len(xs[0])):
                att *= attenuation([xs[0][j]],[0])

            # Fix other particles to their fixed locations
            for j in range(1, len(xs)):
                xj   = xs[j]
                att *= attenuation(xj, fixed(j))

            # Bin this x, y point
            xi = int(RES*(xs[0][0] + SCALE) / (2*SCALE))
            yi = int(RES*(xs[0][1] + SCALE) / (2*SCALE))

            try:
                bins[xi][yi]   += w * att
                datapoints     += 1
                eff_datapoints += att
            except IndexError:
                continue

print("Datapoints          : {0}".format(datapoints))
print("Effective datapoints: {0} ({1} %)".format(eff_datapoints, 100.0*eff_datapoints/datapoints))

bins_file = open("binned","w")
for i in read_fixed:
    x = read_fixed[i]
    s = ",".join([str(xi) for xi in x])
    s += ";"
    bins_file.write(s)
bins_file.write("\n")

for i in range(0, RES):
    for j in range(0, RES):
        bins_file.write(",{0}".format(bins[i][j]))
    bins_file.write(";")
