# 
#     XDMC
#     Copyright (C) 2019 Michael Hutcheon (email mjh261@cam.ac.uk)
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     For a copy of the GNU General Public License see <https://www.gnu.org/licenses/>.
# 
import matplotlib.pyplot as plt
import numpy as np

from numpy.polynomial.hermite import hermval

def factorial(n):
    if n == 0: return 1
    else: return n*factorial(n-1)

def harmonic_osc(x, n, omega=1.0):
    # Returns the n^th exited state of a 1D harmonic
    # oscillator evaluated at x
    xs = np.sqrt(omega) * x
    c  = [1.0 if i == n else 0.0 for i in range(0, n+1)]
    norm = (2**n * factorial(n))**(-0.5) * (omega/np.pi)**(0.25)
    return norm * np.exp(-0.5*xs**2) * hermval(xs, c) 

def exited_indicies(n, dim):
    # Returns the one-dimensional exitations in the n^th
    # exited state in dim dimensions
    ret = [0 for d in range(0,dim)]
    return ret

def all_indicies(max_index=4, dims=1, set_so_far=[[]]):
    # Return all sets of integers up to size dims, that
    # contain integers up to max_index

    if len(set_so_far[0]) >= dims:
        return set_so_far

    new_set = []
    for l in set_so_far:
        for j in range(0, max_index):
            new_l = list(l)
            new_l.append(j)
            new_set.append(new_l)

    return all_indicies(max_index=max_index, dims=dims, set_so_far=new_set)

def exited_states(max_state, dims):
    # Return the extied state indicies in order of
    # increasing energy
    inds = all_indicies(max_index=max_state, dims=dims)
    inds.sort(key = lambda i:sum(i))
    return inds

def psi(x, n, exitation_indicies):
    # Returns the n^th exited state of a harmonic oscillator
    # in len(x) dimensions, evaluated at x
    return np.prod([harmonic_osc(x[i], exitation_indicies[n][i]) for i in range(0, len(x))])

def slater_det(x, exitation_indicies):
    # Returns the slater determinant of n = len(x) states
    # evaluated at the many-particle configuration x
    mat = np.zeros((len(x),len(x)))
    for i in range(0, len(x)):
        for j in range(0, len(x)): 
            mat[i][j] = psi(x[i], j, exitation_indicies)
    return np.linalg.det(mat)

def plot_2d_nodes(scale=4.0, res=100, other_particles = [[0,0]]):

    # Construct the configuration
    config = [[0,0]]
    config.extend(other_particles)

    # Move the first particle around, and work out wavefunction
    xs = np.linspace(-scale, scale, res)
    ys = np.linspace(-scale, scale, res)
    psis = np.zeros((len(xs), len(ys)))

    min_psi = np.inf
    max_psi = -np.inf

    exitations = exited_states(len(other_particles)+1, 2)
    for xi, x in enumerate(xs):
        for yi, y in enumerate(ys):
            config[0][0] = x
            config[0][1] = y
            psi          = slater_det(config, exitations)
            psis[yi][xi] = psi
            if psi < min_psi: min_psi = psi
            if psi > max_psi: max_psi = psi

    max_abs = max(abs(min_psi), abs(max_psi))
    levels = np.linspace(-max_abs, max_abs, 41)

    #plt.imshow(psis, extent=(scale, -scale, scale, -scale))
    plt.contour(xs, ys, psis, levels=levels)
    plt.scatter([p[0] for p in config[1:]], [p[1] for p in config[1:]])
    plt.gca().set_aspect(1.0)

for n in range(0, 50):
    print("{0}/50".format(n))
    plot_2d_nodes(other_particles=[8*(np.random.rand(2)-0.5) for i in range(0,10)])
    plt.savefig("nodes_{0}.png".format(n))
    plt.clf()

