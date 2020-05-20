# 
#     XDMC
#     Copyright (C) Michael Hutcheon (email mjh261@cam.ac.uk)
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
import numpy as np
import matplotlib.pyplot as plt

# Parameters
PARTICLES    = 1
DIMENSIONS   = 1
TAU          = 0.01
WALKERS      = 100
ITERATIONS   = 1000
TRIAL_ENERGY = 0.0

class particle:
        
    def __init__(self):
        
        # Initialize at the origin
        self.position = np.zeros(DIMENSIONS)

    def copy(self):
        
        # Return an exact copy of this particle
        p = particle()
        p.position = np.array(self.position)
        return p

    def potential(self):
        
        # Return the potential evaluated at my location
        return 0.5*np.linalg.norm(self.position)**2

    def make_move(self, tau=TAU):
        
        # Make a move according to our proposal
        # distribution
        self.position += np.random.normal(size=DIMENSIONS, scale=tau**0.5)

    def sq_distance_to(self, other):
        
        # Return the squared distance |x-y|^2 between
        # x=self y=other
        return np.linalg.norm(self.position-other.position)**2

class walker:
    
    def __init__(self):

        # Weight starts as 1.0
        self.weight = 1.0

        # Initialize the particles
        self.particles = []
        for n in range(0, PARTICLES):
            self.particles.append(particle())

    def copy(self):

        # Return an exact copy of this walker
        w           = walker()
        w.weight    = self.weight
        w.particles = [p.copy() for p in self.particles]
        return w

    def branch_copy(self):
        
        # Return a copy of this walker that has
        # been spawned by branching
        w        = self.copy()
        w.weight = np.sign(self.weight)
        return w

    def make_move(self):
        
        # Make a move according to our proposal distribution
        for p in self.particles: p.make_move()

    def sq_distance_to(self, other):

        # Return the squared distance in configuration space
        # |x-y|^2 between x=self and y=other
        r2 = 0
        for i in range(0, len(self.particles)):
            r2 += self.particles[i].sq_distance_to(other.particles[i])
        return r2


    def potential(self):
        
        # Evaluate the potential
        pot = 0.0
        for p in self.particles: pot += p.potential()
        return pot

    def greens_function(self, other, tau=TAU):
        
        # Evaluate the greens function G(x, y, \tau)
        # between two configurations x=self y=other
        r2 = self.sq_distance_to(other)
        ex  = -0.5*r2/tau                                     # Diffusive part
        ex += -0.5*tau*(self.potential() + other.potential()) # Potential part
        ex += tau * TRIAL_ENERGY                              # Normalization part
        return np.exp(ex)

def propagated_wavefunction(walker, walkers_to_propagate):
    
    # Evaluate the propagated wavefunction of the given
    # walker set at the configuration of the given walker 
    # \psi(x) = \sum_i w_i G(x, x_i, \tau)
    gf = 0
    for w in walkers_to_propagate:
        gf += w.weight * w.greens_function(walker)
    return gf

# Data saved over iterations
data = []

# Initialize walkers
walkers = [walker() for i in range(0, WALKERS)]

# Run DMC iterations
for dmc_iter in range(0, ITERATIONS):

    # Save last iteration walkers
    walkers_last = [w.copy() for w in walkers]

    # Propagate
    accepted_moves = 0
    rejected_moves = 0
    for i, w in enumerate(walkers):

        # Make proposal move
        w_new = w.copy()
        w_new.make_move()

        # Get propagated wavefunction before/after
        wfn_before = propagated_wavefunction(w, walkers_last)
        wfn_after  = propagated_wavefunction(w_new, walkers_last)

        # Move acceptance probability
        accept_prob = abs(wfn_after/wfn_before)

        if np.random.rand() < accept_prob:
            
            # Accept the move
            walkers[i] = w_new
            walkers[i].weight = np.sign(wfn_after)
            accepted_moves += 1

        else:

            # Reject the move
            walkers[i] = w
            walkers[i].weight = np.sign(wfn_before)
            rejected_moves += 1

    # Branch
    walkers_new = []
    for w in walkers:
    
        # Work out how many walkers should survive at 
        # the configuration of w
        surviving = int(np.random.rand() + abs(w.weight))
        for i in range(0, surviving):
            walkers_new.append(w.branch_copy())

    walkers = walkers_new

    # Data for this iteration
    iter_data = {
        "iter" : dmc_iter,
        "pop"  : int(len(walkers)),
        "acc"  : float(accepted_moves)/float(accepted_moves+rejected_moves),
        "et"   : TRIAL_ENERGY,
        "pot"  : np.mean([w.potential() for w in walkers])
    }
    data.append(iter_data)
    
    # Print iteration results
    fs  = "Iteration {iter}\n"
    fs += "  Population       {pop}\n"
    fs += "  Acceptance ratio {acc}\n"
    fs += "  E_T              {et} \n"
    fs += "  <V>              {pot}\n"

    print(fs.format(**iter_data))

# Plot the results
for i, k in enumerate(data[0]):
    if k == "iter": continue

    plt.subplot(len(data[0])-1,1,i)
    plt.plot([d[k] for d in data])
    plt.ylabel(k)
    plt.xlabel("Iteration")

plt.show()


