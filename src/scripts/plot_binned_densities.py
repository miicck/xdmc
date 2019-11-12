import matplotlib.pyplot as plt
import numpy as np

def read_density(filename):
    d = []
    with open(filename) as f:
        for l in f:
            d.append([float(x) for x in l.split()])
    return np.array(d)

density      = read_density("density.binned")
cond_density = read_density("cond_density.binned")

def plot_density(d):
    plt.imshow(d, extent=(-4,4,-4,4))
    plt.axhline(0, color="black")
    plt.axvline(0, color="black")

plt.subplot(221)
plot_density(density)
plt.subplot(222)
plot_density(cond_density)
plt.subplot(223)
plot_density(cond_density - density)

plt.show()
