import parser
import matplotlib.pyplot as plt
import numpy as np
import sys

wfn = parser.parse_wavefunction()
particles = len(wfn)
coords    = len(wfn[0])

if particles != 2:
	print "Error ", sys.argv[0], " only works with 2 particle systems!"
	quit()

if coords != 1:
	print "Error ", sys.argv[0], " only works with 1 dimensional systems!"
	quit()

xs = [p[1] for p in wfn[0][0][:]]
ys = [p[1] for p in wfn[1][0][:]]
zs = [p[0] for p in wfn[0][0][:]]

RES = 20
min_lim = min(min(xs), min(ys))
max_lim = max(max(xs), max(ys))
xn  = np.linspace(min_lim, max_lim, RES)
yn  = np.linspace(min_lim, max_lim, RES)
zn  = [] 
tau = 0.01 
for y in yn:
	row = []
	for x in xn:
		val = 0
		for xi, yi, zi in zip(xs, ys, zs):
			dx = xi - x
			dy = yi - y
			val += zi*np.exp(-(dx**2+dy**2)/(2*tau))
		row.append(val)
	zn.append(row)

xn, yn = np.meshgrid(xn, yn)

plt.subplot(211)
plt.contour(xn, yn, zn)
plt.plot([min_lim, max_lim], [min_lim, max_lim], color="black")

def psi0(x):
        return np.exp(-x**2)

def psi1(x):
        return x * psi0(x)

plt.subplot(212)
xs = np.linspace(min_lim,max_lim,100)
xs, ys = np.meshgrid(xs, xs)
zs = [psi0(x)*psi1(y)-psi0(y)*psi1(x) for x, y in zip(xs,ys)]
plt.contour(xs,ys,zs)


plt.show()
