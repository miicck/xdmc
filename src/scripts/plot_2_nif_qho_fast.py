import nif_plotter
import matplotlib.pyplot as plt
import parser
import sys

wfn = parser.parse_wavefunction(int(sys.argv[1]))
nif_plotter.plot_2nif_fast(wfn)
plt.show()
