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

def norm(x1, x2):
    n = 0
    for i in range(0, len(x1)):
        n += (x2[i] - x1[i])**2
    return n**0.5

rs  = []
drs = []
wfn = parse_wavefunction(sys.argv[1:])
for i in range(0, len(wfn[0])):
    
    w = wfn[0][i]
    for j in range(1, len(wfn)):
        xj = wfn[j][i]

        rs.append(la.norm(xj))
        for k in range(1, j):
            xk = wfn[k][i]
            drs.append(norm(xj,xk))

plt.subplot(221)
plt.hist(rs, bins=1000)

plt.subplot(222)
plt.hist(drs, bins=1000)
plt.show()

