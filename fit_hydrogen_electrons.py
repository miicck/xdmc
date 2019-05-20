import matplotlib.pyplot as plt;
import sys
import numpy as np

from scipy.optimize import curve_fit

# Attempts to fit a hydrogen 1s ground state to the electron
# wavefunction sampled in the "electrons" file.

# The radial probability of a DMC walker in the
# hydrogen ground state. The walkers will sample
# the distribution \psi_0 = np.exp(-r/a), but a 
# phase space factor of r^2 also appears.
def rad_prob(r,a,b):
	return b*r*r*np.exp(-r/a)

# The maximum radius we consider
MAX_R = 20

# read in the radii from the "electrons" file
to_open = "electrons"
if len(sys.argv) > 1: to_open = sys.argv[1]
rs = []
for line in [l for l in open(to_open).read().split("\n") if len(l) > 0]:
	x = np.array([float(i) for i in line.split(",")])
	r = np.linalg.norm(x)
	if r < MAX_R:
		rs.append(r)

# Histogram the radii
y,x,_ = plt.hist(rs, bins=100)
for i in range(0,len(x)-1):
	x[i] = (x[i] + x[i+1])/2
x = x[0:-1]
rgrid = np.linspace(0, MAX_R, 1000)

# Fit to the hydrogen 1s state
par, covar = curve_fit(rad_prob, x, y, p0=[1, y.max()])
yfit = rad_prob(rgrid, *par)
lab = "Hydrogen 1s fit: a0 = " + str(round(par[0],5)) + " +/- "
lab += str(round(np.sqrt(covar[0][0]), 5))
plt.plot(rgrid, yfit, label=lab)
y1s = rad_prob(rgrid, 1, 1)
y1s *= max(yfit)/max(y1s)
plt.plot(rgrid, y1s,label="Hydrogen true 1s: a0 = 1", linestyle=":")
plt.xlabel("Radius (atomic units)")
plt.ylabel("Counts")
plt.legend(loc="best")
plt.show()
