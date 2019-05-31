import parser
import matplotlib.pyplot as plt
import numpy as np
import sys
import nif_plotter as nifp

start = int(sys.argv[1])
end   = int(sys.argv[2])
wavefunction = parser.parse_wavefunction(start, end)

plt.show()
n_saved = []
n = -1
while True:
	n = (n+1)%len(wavefunction)
	nifp.plot_2nif_fast([wavefunction[n]])
	plt.suptitle("Iteration {0}".format(start+n))
	plt.draw()
	plt.pause(0.001)
	if not n in n_saved:
		plt.savefig("iter_{0}".format(start+n))
		n_saved.append(n)
	plt.clf()
