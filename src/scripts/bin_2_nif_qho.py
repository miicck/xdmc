import matplotlib.pyplot as plt
import nif_plotter
import sys

nif_plotter.bin_2nif(sys.argv[1:])
plt.show()
