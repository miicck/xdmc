import parser
import matplotlib.pyplot as plt
import numpy as np
import sys
import nif_plotter as nifp

count = int(sys.argv[1])
wavefunction = parser.parse_wavefunction()
plt.show()

n = -1
while True:
	n = (n+1)%count + len(wavefunction) - count
	nifp.plot_2nif([wavefunction[n]])
	plt.suptitle("Iteration {0}".format(n+1))
	plt.draw()
	plt.pause(0.001)
	plt.savefig("iter_{0}".format(n+1))
	plt.clf()
