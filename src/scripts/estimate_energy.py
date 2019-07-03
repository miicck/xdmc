import sys
from parser import parse_evolution
import numpy as np

def exp_fit():
	from scipy.optimize import curve_fit
	for name, energies in zip(*parse_evolution()):
		if not "Trial energy" in name: continue

		# Fit energy vs iteration to exponential decay
		# using initial guess informed by data
		to_fit = lambda n, a, b, c: a + b*np.exp(-n/c)
		p0 = [energies[-1], energies[0] - energies[-1], len(energies)/4]
		bounds = [[min(energies), -np.inf, 1],
			  [max(energies),  np.inf, 10*len(energies)]]
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
			plt.axhline(0, color="black")
			plt.ylabel("Residuals (Hartree)")

			plt.show()

def multi_exp_fit(exponents=3):
	from scipy.optimize import curve_fit
	for name, energies in zip(*parse_evolution()):
		if not "Trial energy" in name: continue
		to_fit = lambda n, a, b, c: a + b*np.exp(-n/c)

		residuals  = np.array(energies)
		parameters = []
		for i in range(0, exponents):
			p0 = [residuals[-1], residuals[0] - residuals[-1], len(residuals)/4]
			bounds = [[min(residuals), -np.inf, 1],
				  [max(residuals), np.inf,  10*len(residuals)]]
			par, cov = curve_fit(to_fit, range(0, len(residuals)),
					     residuals, p0, bounds=bounds)
			fit = [to_fit(n, *par) for n in range(0, len(residuals))]
			residuals -= fit
			parameters.append(par)

		model = lambda n : sum(to_fit(n, *p) for p in parameters)
		print model(10e10)

		if plot:
			fit = [model(n) for n in range(0, len(energies))]
			ext = [model(n) for n in range(0, len(energies)*10)]
			plt.subplot(211)
			plt.plot(energies)
			plt.plot(fit)
			plt.plot(ext)
			plt.subplot(212)
			plt.plot(np.array(energies)-np.array(fit))
			plt.axhline(0, color="black")
			plt.show()

def average():
	for name, energies in zip(*parse_evolution()):
		if not "Trial energy" in name: continue
		start = int(sys.argv[2])

		energies = energies[start:]
		print np.mean(energies), " +/- ", np.std(energies)
		if plot:
			plt.plot(energies)
			plt.show()

methods = {
"exp_fit" : exp_fit,
"multi_exp_fit" : multi_exp_fit,
"average" : average
}

if not sys.argv[1] in methods:
	new_argv = [sys.argv[0], "exp_fit"]
	new_argv.extend(sys.argv[1:])
	sys.argv = new_argv

plot = "plot" in sys.argv[2:]
if plot: import matplotlib.pyplot as plt

methods[sys.argv[1]]()
