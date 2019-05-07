import matplotlib.pyplot as plt;
import numpy as np

from scipy.optimize import curve_fit

# Attempts to fit a hydrogen 1s ground state to the electron
# wavefunction sampled in the "electrons" file.

# The radial probability density of a 1s state
# including the r^2 prefactor (parameterised)
def rad_prob(r,a,b):
	return b*r*r*np.exp(-2*r/a)/np.exp(-2)

# The maximum radius we consider
MAX_R = 20

# read in the radii from the "electrons" file
rs = []
for line in [l for l in open("electrons").read().split("\n") if len(l) > 0]:
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
plt.plot(rgrid, y1s,label="Hydrogen true 1s: a0 = 1")
plt.xlabel("Radius (atomic units)")
plt.ylabel("Counts")
plt.legend(loc="best")
plt.show()
