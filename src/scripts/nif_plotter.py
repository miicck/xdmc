import parser
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

def label_plot():
	plt.xlabel("Particle 1 position")
        plt.ylabel("Particle 2 position")
	plt.gca().set_aspect(1.0)

def get_weight_cmap(ws=np.linspace(-1,1,1000), alpha=True):
	cmap = plt.cm.seismic
	if not alpha: return cmap
	my_cmap = cmap(np.arange(cmap.N))
	max_alpha = min(0.5, 500.0/len(ws))
	my_cmap[:,-1] = max_alpha*abs(np.linspace(-1, 1, cmap.N))
	my_cmap = ListedColormap(my_cmap)
	return my_cmap

def plot_analytic_2nif_harmonic(min_lim, max_lim, LEVELS, alpha=1.0):
	def psi0(x):
                return np.exp(-x**2/2)

        def psi1(x):
                return x * psi0(x)

        xs = np.linspace(min_lim,max_lim,100)
        xs, ys = np.meshgrid(xs, xs)
        zs = [psi0(x)*psi1(y)-psi0(y)*psi1(x) for x, y in zip(xs,ys)]
        plt.contour(xs,ys,zs,LEVELS,alpha=alpha, cmap=get_weight_cmap(alpha=False))
	plt.plot([min_lim, max_lim], [min_lim, max_lim], color="black")
	label_plot()
	plt.title("Analytic solution")

def plot_2nif_fast(start_iter, end_iter):
	
	wavefunction = parser.parse_wavefunction(start_iter, end_iter)
	wfn = parser.transpose_wavefunction(wavefunction)
	xs = [x[0] for x in wfn[1]]
	ys = [x[0] for x in wfn[2]]
	zs = wfn[0]

	plot_analytic_2nif_harmonic(-4,4,5)
	plt.scatter(xs,ys,c=zs,cmap=get_weight_cmap(zs))
	plt.title("Walker positions\n(With analytic solution overlaid)")
	plt.xlim([-4,4])
	plt.ylim([-4,4])
	label_plot()

def plot_2nif(start_iter, end_iter):
	wavefunction = parser.parse_wavefunction(start_iter, end_iter)
	iterations = len(wavefunction)
	wfn = parser.transpose_wavefunction(wavefunction)

	xs = [x[0] for x in wfn[1]]
	ys = [x[0] for x in wfn[2]]
	zs = wfn[0]

	fs = "{0} walkers from {1} dmc iteration(s) {2} to {3}"
	plt.suptitle(fs.format(len(zs), iterations, start_iter, end_iter))

	RES    = 40
	LEVELS = 40
	min_lim = -4
	max_lim = 4

	plt.subplot(2,2,3)
	alpha = min(0.2, 200.0/len(zs))

	plt.scatter(xs, ys, c=zs, cmap=get_weight_cmap(zs))
	plt.xlim([min_lim,max_lim])
	plt.ylim([min_lim,max_lim])
	plt.plot([min_lim, max_lim], [min_lim, max_lim])
	plt.title("Walker positions")
	label_plot()

	xn = np.linspace(min_lim, max_lim, RES)
	yn = np.linspace(min_lim, max_lim, RES)
	
	for iplt, tau in enumerate([0.05, 1.0]):

		zn  = np.zeros((len(yn), len(xn)))
		for iy, y in enumerate(yn):
			for ix, x in enumerate(xn):
				dr2 = (xs - x)**2 + (ys - y)**2
				val = np.sum(zs*np.exp(-dr2/(2*tau)))
				zn[iy][ix] = val

		xplot, yplot = np.meshgrid(xn, yn)


		plt.subplot(2,2,1+iplt)
		plt.contour(xplot, yplot, zn, LEVELS, cmap=get_weight_cmap(zs, alpha=False))
		plt.plot([min_lim, max_lim], [min_lim, max_lim], color="black")
		plt.title("Averaged wavefunction\n(Using Guassian kernel with Tau = {0})".format(tau))
		label_plot()

	plt.subplot(224)
	plot_analytic_2nif_harmonic(min_lim, max_lim, LEVELS)
