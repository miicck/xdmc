import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import parser
import sys

# Pick start and end indicies
start = int(sys.argv[1])
end   = int(sys.argv[2])

plt.show()

elev = -30
azim = -133

n = -1
while True:

	# View from same angle
	ax = plt.gcf().add_subplot(111, projection="3d")
	ax.elev = elev
	ax.azim = azim

	n = (n+1)%(end-start)

	# Read wavefunction files
	ws, xs, ys, zs = parser.parse_wavefunction(start+n, start+n+1)

	# Convert [[x], [x], [x]..] to [x,x,x..]
	xs = [x[0] for x in xs]
	ys = [y[0] for y in ys]
	zs = [z[0] for z in zs]

	# Scatter plot the data
	plt.suptitle("Iteration {0}".format(start+n))
	ax.scatter(xs, ys, zs, c=ws)

	plt.draw()
	plt.pause(0.001)
	elev = ax.elev
	azim = ax.azim

	plt.clf()
