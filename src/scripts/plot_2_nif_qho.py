import parser
import matplotlib.pyplot as plt
import nif_plotter as nifp
import sys

count = int(sys.argv[1])
nifp.plot_2nif(parser.parse_wavefunction()[-count:])
plt.show()
