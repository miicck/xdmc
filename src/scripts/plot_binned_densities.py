import matplotlib.pyplot as plt
import numpy as np

density      = np.load("density.binned.npy")**2
cond_density = np.load("cond_density.binned.npy")**2

def plot_density(d):
    if len(d.shape) == 1:
        x = np.linspace(-4,4,len(d))
        plt.plot(x,d)
    if len(d.shape) == 2:
        plt.imshow(d, extent=(-4,4,-4,4))
        plt.axhline(0, color="black")
        plt.axvline(0, color="black")

plt.subplot(221)
plot_density(density)
plt.subplot(222)
plot_density(cond_density)
plt.subplot(223)
plot_density(cond_density/density)

plt.show()
