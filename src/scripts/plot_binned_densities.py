import matplotlib.pyplot as plt
import numpy as np

density      = np.load("density.binned.npy")
cond_density = np.load("cond_density.binned.npy")

def plot_density(d):
    if len(d.shape) == 1:
        x = np.linspace(-4,4,len(d))
        plt.plot(x,d)
    elif len(d.shape) == 2:
        plt.imshow(d, extent=(-4,4,-4,4))
        plt.axhline(0, color="black")
        plt.axvline(0, color="black")
    elif len(d.shape) == 3:
        h = int(len(d[0,0,:])/2)
        plt.imshow(d[:,:,h], extent=(-4,4,-4,4))
        plt.axhline(0, color="black")
        plt.axvline(0, color="black")
    else:
        raise ValueError("Cannot plot data of shape {0}!".format(d.shape))

plt.subplot(221)
plot_density(density)
plt.subplot(222)
plot_density(cond_density)
plt.subplot(223)
plot_density(cond_density-density)

plt.show()
