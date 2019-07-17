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

def analytic(xs,ys):
        psi0 = lambda x : np.exp(-x**2/2)
        psi1 = lambda x : x * psi0(x)
        return np.array([(psi0(x)*psi1(y)-psi0(y)*psi1(x)) for x, y in zip(xs,ys)])

def plot_analytic_2nif_harmonic(min_lim, max_lim, LEVELS, alpha=1.0):
        xs = np.linspace(min_lim,max_lim,100)
        xs, ys = np.meshgrid(xs, xs)
        zs = analytic(xs, ys)
        plt.contour(xs,ys,zs,LEVELS,alpha=alpha, cmap=get_weight_cmap(alpha=False))
        plt.plot([min_lim, max_lim], [min_lim, max_lim], color="black")
        label_plot()
        plt.title("Analytic solution")

def plot_2nif_fast(start_iter, end_iter):
        
        wfn = parser.parse_wavefunction(start_iter, end_iter)
        xs = [x[0] for x in wfn[1]]
        ys = [x[0] for x in wfn[2]]
        zs = wfn[0]

        plot_analytic_2nif_harmonic(-4,4,5)
        plt.scatter(xs,ys,c=zs,cmap=get_weight_cmap(zs))
        plt.title("Walker positions\n(With analytic solution overlaid)\n")
        plt.xlim([-4,4])
        plt.ylim([-4,4])
        label_plot()

def average_r2(xs, ys, bins):
        
        tot_r2 = 0
        tot_w  = 0
        for i in range(0, len(xs)):
                for j in range(0, len(xs[i])):
                        w = bins[i][j]**2
                        tot_r2 += w*(xs[i][j]**2 + ys[i][j]**2)
                        tot_w  += w
        return tot_r2/tot_w
                        

def bin_2nif(start_iter, end_iter, RES=40):
        wfn = parser.parse_wavefunction(start_iter, end_iter)
        xs = np.array([x[0] for x in wfn[1]])
        ys = np.array([x[0] for x in wfn[2]])
        ws = wfn[0]

        MAX  = 4.0
        MIN  = -4.0
        bins = np.zeros((RES,RES))

        for x,y,w in zip(xs,ys,ws):
                xi = int(RES*(x - MIN)/(MAX - MIN))
                yi = int(RES*(y - MIN)/(MAX - MIN))
                if xi < 0:    continue
                if yi < 0:    continue
                if xi >= RES: continue
                if yi >= RES: continue
                bins[xi][yi] += w
        
        xs = np.linspace(MIN, MAX, RES)
        ys = np.linspace(MIN, MAX, RES)
        xs, ys = np.meshgrid(xs, ys)

        plt.subplot(221)
        plt.contour(xs, ys, bins, 40)
        plt.xlabel("Particle 1 position")
        plt.ylabel("Particle 2 position")
        plt.gca().title.set_text("DMC wavefunction\n<r^2> = {0}".format(average_r2(xs,ys,bins)))

        plt.subplot(222)
        anal_sol = analytic(xs,ys)
        plt.contour(xs, ys, anal_sol, 40)
        plt.xlabel("Particle 1 position")
        plt.ylabel("Particle 2 position")
        plt.gca().title.set_text("Analytic wavefunction\n<r^2> = {0}".format(average_r2(xs,ys,anal_sol)))

def plot_2nif(start_iter, end_iter):
        wfn = parser.parse_wavefunction(start_iter, end_iter)
        iterations = end_iter - start_iter
        
        print(len(wfn[1]))
        xs = np.array([x[0] for x in wfn[1]])
        ys = np.array([x[0] for x in wfn[2]])
        zs = wfn[0]
        zs_an = analytic(xs, ys)

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
        
        for iplt, tau in enumerate([0.01, 1.0]):

                zn    = np.zeros((len(yn), len(xn)))
                zn_an = np.zeros((len(yn), len(xn)))
                for iy, y in enumerate(yn):
                        for ix, x in enumerate(xn):
                                dr2 = (xs - x)**2 + (ys - y)**2
                                val = np.sum(zs*np.exp(-dr2/(2*tau)))
                                val_an = np.sum(zs_an*np.exp(-dr2/(2*tau)))
                                zn[iy][ix] = val
                                zn_an[iy][ix] = val_an

                xplot, yplot = np.meshgrid(xn, yn)

                plt.subplot(2,4,1+2*iplt)
                plt.contour(xplot, yplot, zn_an, LEVELS, cmap=get_weight_cmap(zs, alpha=False))
                plt.title("Averaged analytic wavefunction\n(Using Guassian kernel with Tau = {0})".format(tau))
                plt.plot([min_lim, max_lim], [min_lim, max_lim], color="black")
                label_plot()

                plt.subplot(2,4,2+2*iplt)
                plt.contour(xplot, yplot, zn, LEVELS, cmap=get_weight_cmap(zs, alpha=False))
                plt.title("Averaged QMC wavefunction\n(Using Guassian kernel with Tau = {0})".format(tau))
                plt.plot([min_lim, max_lim], [min_lim, max_lim], color="black")
                label_plot()

        plt.subplot(224)
        plot_analytic_2nif_harmonic(min_lim, max_lim, LEVELS)
