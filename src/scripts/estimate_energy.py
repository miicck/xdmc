import matplotlib.pyplot as plt
import sys
from parser import parse_evolution
from scipy.optimize import curve_fit
import numpy as np

plot = False
if "plot" in sys.argv[1:]: plot = True

y_axes, data = parse_evolution()
for i in range(0, len(y_axes)):
	if not "Trial energy" in y_axes[i]: continue

	# Fit energy vs iteration to exponential decay
	# using initial guess informed by data
	energies = data[i]
	to_fit = lambda n, a, b, c: a + b*np.exp(-n/c)
	p0 = [energies[-1], energies[0] - energies[-1], len(energies)/4]
	bounds = [[min(energies), -np.inf, 1],
		  [max(energies),  np.inf, len(energies)]]
	par, cov = curve_fit(to_fit, range(0, len(energies)), energies, p0, bounds=bounds)
	residuals = [e - to_fit(n, *par) for n, e in enumerate(energies)]
	print par[0] ,"+/-", np.std(residuals)

	if plot:
		
		# Plot the data, fit and residuals
		plt.subplot(211)
		plt.plot(energies, label="DMC")
		xs = np.linspace(0, len(energies), 1000)
		ys = [to_fit(x, *par) for x in xs]
		plt.plot(xs, ys, label="Fit")
		plt.ylabel("Trial energy (Hartree)")
		plt.legend()

		plt.subplot(212)
		plt.plot(residuals)
		plt.ylabel("Residuals (Hartree)")

		plt.show()
