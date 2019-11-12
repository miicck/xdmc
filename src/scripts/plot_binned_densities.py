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

plt.subplot(221)
plt.imshow(density)

plt.subplot(222)
plt.imshow(cond_density)

plt.subplot(223)
plt.imshow(cond_density - density)

plt.show()
