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
