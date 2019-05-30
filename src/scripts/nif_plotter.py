import parser
import matplotlib.pyplot as plt
import numpy as np

def plot_2nif(wavefunction):
	iterations = len(wavefunction)
	wfn = parser.transpose_wavefunction(wavefunction)

	particles = len(wfn)-1
	coords    = len(wfn[1][0])

	if particles != 2:
		print "Error ", sys.argv[0], " only works with 2 particle systems!"
		quit()

	if coords != 1:
		print "Error ", sys.argv[0], " only works with 1 dimensional systems!"
		quit()

	xs = [x[0] for x in wfn[1]]
	ys = [x[0] for x in wfn[2]]
	zs = wfn[0]
	plt.suptitle("{0} walkers from {1} dmc iteration(s)".format(len(zs), iterations))

	RES = 20
	min_lim = min(min(xs), min(ys))
	max_lim = max(max(xs), max(ys))
	xn  = np.linspace(min_lim, max_lim, RES)
	yn  = np.linspace(min_lim, max_lim, RES)

	for iplt, tau in enumerate([0.01, 0.1, 1.0]):

		zn  = [] 
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

		xplot, yplot = np.meshgrid(xn, yn)

		plt.subplot(2,2,1+iplt)
		plt.contour(xplot, yplot, zn)
		plt.plot([min_lim, max_lim], [min_lim, max_lim], color="black")

	def psi0(x):
		return np.exp(-x**2)

	def psi1(x):
		return x * psi0(x)

	xs = np.linspace(min_lim,max_lim,100)
	xs, ys = np.meshgrid(xs, xs)
	zs = [psi0(x)*psi1(y)-psi0(y)*psi1(x) for x, y in zip(xs,ys)]

	plt.subplot(224)
	plt.contour(xs,ys,zs)
