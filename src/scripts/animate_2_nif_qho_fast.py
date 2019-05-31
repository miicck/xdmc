import parser
import matplotlib.pyplot as plt
import numpy as np
import sys
import nif_plotter as nifp

start = int(sys.argv[1])
end   = int(sys.argv[2])
if "save" in sys.argv: save = True
else: save = False

plt.show()
n_saved = []
n = -1

while True:

	n = (n+1)%(end-start)
	nifp.plot_2nif_fast(n, n+1)
	plt.suptitle("Iteration {0}".format(start+n))
	plt.draw()
	plt.pause(0.001)

	if not n in n_saved:
		if save:
			plt.savefig("iter_{0}".format(start+n))
		n_saved.append(n)

	plt.clf()
