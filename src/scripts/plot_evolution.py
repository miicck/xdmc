import sys
import numpy as np
import matplotlib.pyplot as plt
from parser import parse_evolution

# Read in the evolution data
y_axes, data = parse_evolution()
y_size = int(len(data)**0.5)
x_size = 0
while x_size * y_size < len(data): x_size += 1

x_size = 1
y_size = len(data)

e_targ = None
for arg in sys.argv:
    if arg.startswith("-e="):
        e_targ = float(arg.split("=")[-1])
        break

v_targ = None
for arg in sys.argv:
    if arg.startswith("-v="):
        v_targ = float(arg.split("=")[-1])
        break

log_scales = ["Cancelled weight"]

# Plot each dataseries on its own subplot
for i, d in enumerate(data):
    plt.subplot(y_size,x_size,i+1)
    plt.plot(d)
    plt.xlabel("Iteration")
    plt.ylabel(y_axes[i])

    # Set the axis scales
    if y_axes[i] in log_scales:
        plt.yscale("log")

    elif "norm" in sys.argv:
        equil_n = int(len(d)/2)
        mean = np.mean(d[equil_n:])
        std  = np.std(d[equil_n:])
        plt.ylim([mean - 4*std, mean + 4*std])


    # Plot additional things
    if "Average weight" in y_axes[i]:
        plt.axhline(0, color="black")

    elif "Trial energy" in y_axes[i]:
        if e_targ != None:
            plt.axhline(e_targ, color="red")

    elif "<V>" in y_axes[i]:
        if v_targ != None:
            plt.axhline(v_targ, color="red")

# Spread the subplots out, and show them
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
