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
from parser import parse_wavefunction
import numpy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import sys

# Resolution/range of wavefunction plots
RES   = 100
RANGE = 3

# Get the wavefunction
wfn = parse_wavefunction(sys.argv[1:])

# Fixed locations of particles 2 and 3
x2 = [0.5, 0, 0]
x3 = [-0.5,0, 0]

bins = np.zeros((RES,RES))

def attenuation(x,y,scale=1.0):
    n = 0
    for i in range(0, len(x)):
        n += ((x[i] - y[i])/scale)**2
    return np.exp(-n)

for i in range(0, len(wfn[0])):
    
    # Weight
    w = wfn[0][i]

    # Location of particle 1
    x1 = wfn[1][i]

    # Apply z attenuation
    att = attenuation([0], [x1[2]])

    for j in range(2, len(wfn)):
    
        # Location of particle j
        xj = wfn[j][i]
        
        # Attenuation due to this
        # particle not being at desired location
        if   j == 2: att *= attenuation(xj, x2)
        elif j == 3: att *= attenuation(xj, x3)

    # x,y indicies on plot
    x_index = int(RES*(x1[0] + RANGE)/(2*RANGE))
    y_index = int(RES*(x1[1] + RANGE)/(2*RANGE))

    try:
        bins[x_index][y_index] += w*att
    except IndexError:
        continue

plt.imshow(bins, extent=(-RANGE,RANGE,-RANGE,RANGE))
#plt.scatter(
#    [x[0] for x in [x2,x3]],
#    [x[1] for x in [x2,x3]])

plt.show()


