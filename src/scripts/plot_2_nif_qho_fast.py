import nif_plotter
import matplotlib.pyplot as plt
import sys

font = {'size' : 16}
plt.rc('font', **font)

start = int(sys.argv[1])
end   = int(sys.argv[2])
nif_plotter.plot_2nif_fast(start, end)
plt.show()
