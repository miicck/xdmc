import matplotlib.pyplot as plt
import numpy as np
import parser
import sys

BINS = 50
MIN_U  = -4
MAX_U  =  4
MIN_X  = -4
MAX_X  =  4
plt.rc("text", usetex=True)
plt.rc("font", size=14)

# The 2d projection that we will visualize
# (looking along the (111) direction
def U(x,y,z): return (x-y)/np.sqrt(2)
def V(x,y,z): return (2*z-x-y)/np.sqrt(6)

# Plot a symmetry line in the 
# x,y,z cartesian direction
def plot_symmetry_line(x,y,z):
    global MAX_U
    us = np.array([-U(x,y,z), U(x,y,z)])
    vs = np.array([-V(x,y,z), V(x,y,z)])
    scale = max(max(us), max(vs))
    us *= MAX_U/scale
    vs *= MAX_U/scale
    plt.plot(us, vs, color="black", linestyle=":")

def psi0(x): return np.exp(-x**2/2.0)
def psi1(x): return np.sqrt(2.0)*x*np.exp(-x**2/2.0)
def psi2(x): return (2*x**2-1)*np.exp(-x**2/2.0)/np.sqrt(2.0)

def psi(x,y,z):
    r    = [x,y,z]
    psis = [psi0, psi1, psi2]
    mat   = np.zeros((3,3))
    for i in range(0,3):
        for j in range(0,3):
            mat[i][j] = psis[i](r[j])
    return np.linalg.det(mat)

# Returns \psi(u, v) = \int \psi(x,y,z) \delta(u(x,y,z)-u) \delta(v(x,y,z)-v)
def psi_uv(u, v):
    d = np.array([1.0,1.0,1.0])
    y =  MIN_X
    x =  np.sqrt(2.0)*u + y
    z = (np.sqrt(6.0)*v + x + y)/2.0
    p = np.array([x,y,z])
    return np.trapz([psi(*(p+d*l)) for l in np.linspace(0,MAX_X-MIN_X,20)])

def plot_psi():
    bins = np.zeros((BINS, BINS))
    us = np.linspace(MIN_U, MAX_U, BINS)
    vs = np.linspace(MIN_U, MAX_U, BINS)
    for ui, u in enumerate(us):
        for vi, v in enumerate(vs):
            bins[vi][ui] = psi_uv(u, v)
    plt.contour(us, vs, bins, 11)
    plt.xlabel(r"$u = \frac{x - y}{\sqrt{2}}$")
    plt.ylabel(r"$v = \frac{2z - x - y}{\sqrt{6}}$")
    plt.gca().title.set_text(r"$\langle r^2 \rangle = 6.0$")

    # Plot symmetry lines
    plot_symmetry_line(0,1,1)
    plot_symmetry_line(1,0,1)
    plot_symmetry_line(1,1,0)

def project_wavefunction(wfn):

    global BINS, MIN_U, MAX_U

    # Range in U, V space to plot
    bins = np.zeros((BINS,BINS))

    av_r2 = 0
    tot_w = 0

    if (len(wfn) == 0):
        print("Wavefunction has no samples, skipping...")
        return

    # Bin the wavefunction into this range
    print("Projecting wavefunction of length", len(wfn))
    for w, x, y, z in wfn:

        av_r2 += abs(w)*(x[0]*x[0]+y[0]*y[0]+z[0]*z[0])
        tot_w += abs(w)

        # Work out the u coordinate
        u = U(x[0],y[0],z[0])
        ui = int(BINS*(u-MIN_U)/(MAX_U-MIN_U))
        if ui < 0: continue
        if ui >= BINS: continue

        # Work out the v coordinate
        v = V(x[0],y[0],z[0])
        vi = int(BINS*(v-MIN_U)/(MAX_U-MIN_U))
        if vi < 0: continue
        if vi >= BINS: continue

        bins[vi, ui] += w

    if tot_w == 0:
        print("Error, wavefunction has zero weight!")

    av_r2 /= tot_w

    # Plot the resulting contours
    us, vs = np.meshgrid(np.linspace(MIN_U,MAX_U,BINS), np.linspace(MIN_U,MAX_U,BINS))
    plt.contour(us, vs, bins, 11)
    plt.xlabel(r"$u = \frac{x - y}{\sqrt{2}}$")
    plt.ylabel(r"$v = \frac{2z - x - y}{\sqrt{6}}$")
    plt.gca().title.set_text(r"$\langle r^2 \rangle = {0}$".format(av_r2))

    # Plot symmetry lines
    plot_symmetry_line(0,1,1)
    plot_symmetry_line(1,0,1)
    plot_symmetry_line(1,1,0)

def project_nodal_surface(wfn):

    print("Projecting nodal surface of length", len(wfn))
    us = []
    vs = []
    for w, x, y, z in wfn:
        us.append(U(x[0],y[0],z[0]))
        vs.append(V(x[0],y[0],z[0]))

    plt.scatter(us,vs, alpha=0.1)

# Plot DMC samples
plt.subplot(121)
wfn = list(zip(*parser.parse_wavefunction(sys.argv[1:])))
#if "nodal_surface" in sys.argv[1:]:
#    project_nodal_surface(wfn)
#else:
project_wavefunction(wfn)

# Plot analytic
plt.subplot(122)
plot_psi()

plt.show()
