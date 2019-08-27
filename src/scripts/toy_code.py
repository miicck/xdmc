import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from math import erf

def extern_potential(particle):
        # Harmonic well
        return 0.5*la.norm(particle)**2

def potential(walkers):
        # Return the sum of the single 
        # particle potentials
        ret = []
        for w in walkers:
                pot = 0
                for p in w:
                        pot += extern_potential(p)
                ret.append(pot)
        return np.array(ret)

def cancel_func(w1, w2, tau):
        return erf(la.norm(w1-w2)/(2*np.sqrt(2*tau)))

# DMC algorithm setup
DIM         = 1
PARTICLES   = 2
DMC_WALKERS = 500
DMC_ITERS   = 500
TAU         = 0.01
E_T         = 0

# The particle indicies that 
# obey fermionic exchange
FERMION_GROUPS = [[0,1]]

# Initialize the DMC walkers
data = []
walkers = np.random.normal(scale=1.0,size=(DMC_WALKERS, PARTICLES, DIM))
signs   = np.zeros(DMC_WALKERS)+1.0

for dmc_iter in range(0, DMC_ITERS):

        # Carry out random exchange of walkers
        for i, wi in enumerate(walkers):
                for g in FERMION_GROUPS:

                        # Pick a random exchange from this group
                        i_ex = g[int(np.random.rand()*len(g))]
                        j_ex = g[int(np.random.rand()*len(g))]
                        if i_ex == j_ex: continue

                        # Perform the exchange
                        tmp      = wi[i_ex].copy()
                        wi[i_ex] = wi[j_ex].copy()
                        wi[j_ex] = tmp
                        signs[i] *= -1

        # Carry out diffusion of walkers
        pots_before = potential(walkers)
        walkers    += np.random.normal(scale=np.sqrt(TAU), size=walkers.shape)
        pots_after  = potential(walkers)

        # Caclulate potential renormalization
        weights = np.exp(-TAU*(pots_before + pots_after - 2.0*E_T)/2.0)

        # Apply inter-walker cancellations
        for i, wi in enumerate(walkers):
                for j, wj in enumerate(walkers):
                        if signs[i] == signs[j]: continue
                        weights[i] *= cancel_func(wi, wj, TAU)

        # Apply branching / accumulate averages
        new_walkers = []
        new_signs   = []
        av_pot      = 0

        for i, w in enumerate(walkers):
                
                survive = int(np.random.rand() + weights[i])
                av_pot += pots_after[i] * survive

                for j in range(0, survive):
                        new_walkers.append(w)
                        new_signs.append(signs[i])

        walkers = np.array(new_walkers)
        signs = new_signs

        # Update trial energy to control population
        flw     = float(len(walkers))
        av_pot /= flw
        E_T     = av_pot - np.log(flw/float(DMC_WALKERS))

        data.append((E_T, av_pot, len(walkers)))
        print (dmc_iter, E_T, av_pot, len(walkers))

# Plot the resulting evolution
data = list(zip(*data))
for i, d in enumerate(data):
        plt.subplot(len(data),1,i+1)
        plt.plot(data[i])
AV_ITER = int(DMC_ITERS/2)
print(np.mean(data[0][AV_ITER:]))
plt.show()
