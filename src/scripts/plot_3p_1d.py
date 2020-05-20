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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import parser
import sys

# Pick start and end indicies
start = int(sys.argv[1])
end   = int(sys.argv[2])

# Read wavefunction files
ws, xs, ys, zs = parser.parse_wavefunction(start, end)

# Convert [[x], [x], [x]..] to [x,x,x..]
xs = [x[0] for x in xs]
ys = [y[0] for y in ys]
zs = [z[0] for z in zs]

ax = plt.gcf().add_subplot(111, projection="3d")
ax.scatter(xs, ys, zs, c=ws)

plt.show()


