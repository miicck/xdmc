import parser
import matplotlib.pyplot as plt
import numpy as np
import sys
import nif_plotter as nifp

count = int(sys.argv[1])
wavefunction = parser.parse_wavefunction(count)
plt.show()

n_saved = []
n = -1
while True:
	n = (n+1)%len(wavefunction)
	nifp.plot_2nif_fast([wavefunction[n]])
	plt.suptitle("Iteration {0}".format(n+1))
	plt.draw()
	plt.pause(0.001)
	if not n in n_saved:
		plt.savefig("iter_{0}".format(n+1))
		n_saved.append(n)
	plt.clf()
