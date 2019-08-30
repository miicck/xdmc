import matplotlib.pyplot as plt
import nif_plotter
import sys

start = int(sys.argv[1])
end   = int(sys.argv[2]) if len(sys.argv) > 2 else None

nif_plotter.bin_2nif(start, end)
plt.show()
