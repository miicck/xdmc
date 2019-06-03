import matplotlib.pyplot as plt
import nif_plotter
import sys

start = int(sys.argv[1])
end   = int(sys.argv[2])

nif_plotter.bin_2nif(start, end, RES=20)
plt.show()
