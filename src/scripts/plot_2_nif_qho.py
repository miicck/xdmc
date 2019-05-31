import matplotlib.pyplot as plt
import nif_plotter
import sys

start = int(sys.argv[1])
end   = int(sys.argv[2])
nif_plotter.plot_2nif(start, end)
plt.show()
