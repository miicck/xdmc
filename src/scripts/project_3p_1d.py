import matplotlib.pyplot as plt
import numpy as np
import parser
import sys

BINS = 50
MIN  = -4
MAX  =  4

# The 2d projection that we will visualize
# (looking along the (111) direction
def U(x,y,z): return (x-y)/np.sqrt(2)
def V(x,y,z): return (2*z-x-y)/np.sqrt(6)

# Plot a symmetry line in the 
# x,y,z cartesian direction
def plot_symmetry_line(x,y,z):
	global MAX
	us = np.array([-U(x,y,z), U(x,y,z)])
	vs = np.array([-V(x,y,z), V(x,y,z)])
	scale = max(max(us), max(vs))
	us *= MAX/scale
	vs *= MAX/scale
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

def project_wavefunction(wfn):

	global BINS, MIN, MAX

	# Range in U, V space to plot
	bins = np.zeros((BINS,BINS))

	# Bin the wavefunction into this range
	print "Projecting wavefunction of length", len(wfn)
	for w, x, y, z in wfn:

		# Work out the u coordinate
		u = U(x[0],y[0],z[0])
		ui = int(BINS*(u-MIN)/(MAX-MIN))
		if ui < 0: continue
		if ui >= BINS: continue

		# Work out the v coordinate
		v = V(x[0],y[0],z[0])
		vi = int(BINS*(v-MIN)/(MAX-MIN))
		if vi < 0: continue
		if vi >= BINS: continue

		bins[vi, ui] += w

	# Plot the resulting contours
	us, vs = np.meshgrid(np.linspace(MIN,MAX,BINS), np.linspace(MIN,MAX,BINS))
	plt.contour(us, vs, bins, 10)
	plt.xlabel("u = (x - y)/sqrt(2)")
	plt.ylabel("v = (2z - x - y)/sqrt(6)")

	# Plot symmetry lines
	plot_symmetry_line(0,1,1)
	plot_symmetry_line(1,0,1)
	plot_symmetry_line(1,1,0)

	plt.legend()

# Pick start and end indicies
start = int(sys.argv[1])
end   = int(sys.argv[2])

plt.subplot(121)
wfn = zip(*parser.parse_wavefunction(start, end))
project_wavefunction(wfn)

# Plot the analytic solution
awfn = []
for n in range(0, len(wfn)):
	x,y,z = np.random.rand(3)*(MAX-MIN) + MIN
	awfn.append([psi(x,y,z),[x],[y],[z]])

plt.subplot(122)
project_wavefunction(awfn)

plt.show()
