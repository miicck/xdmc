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
import parser
import matplotlib.pyplot as plt
import numpy as np
import sys
import nif_plotter as nifp

start = int(sys.argv[1])
end   = int(sys.argv[2])
if "save" in sys.argv[3:]: save = True
else: save = False

plt.show()
n_saved = []
n = -1

while True:

	n = (n+1)%(end-start)
	nifp.plot_2nif(n, n+1)
	plt.suptitle("Iteration {0}".format(n+start))
	plt.draw()
	plt.pause(0.001)
	if not n in n_saved:
		if save:
			plt.savefig("iter_{0}".format(n+start))
		n_saved.append(n)
	plt.clf()


