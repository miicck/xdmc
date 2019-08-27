import sys
import numpy as np
from math import floor, log10
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from parser import parse_evolution
    
round_to_n = lambda x, n: round(x, -int(floor(log10(abs(x)))) + (n - 1))

# Parse command line args
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

# Read in the evolution data
y_axes, data = parse_evolution()

log_scales = ["Cancelled weight"]

fig, axes = plt.subplots(len(data),2, gridspec_kw = {"width_ratios":[8,1]})

# Plot each dataseries on its own subplot
for i, d in enumerate(data):

    # Plot dataseries
    axes[i,0].plot(d)
    axes[i,0].set_xlabel("Iteration")

    half = int(len(d)/2)
    mean = round_to_n(np.mean(d[half:]), 4)
    std  = round_to_n(np.std(d[half:]),  4)

    axes[i,0].set_ylabel(y_axes[i]+"\n{0} +/- {1}".format(mean,std))
    axes[i,0].axhline(mean, color="green", linestyle=":")

    # Set the axis scales
    if y_axes[i] in log_scales:
        axes[i,0].set_yscale("log")
        axes[i,1].set_yscale("log")

    elif "norm" in sys.argv:
        axes[i,0].set_ylim([mean - 4*std, mean + 4*std])

    # Plot additional things
    if "Average weight" in y_axes[i]:
        axes[i,0].axhline(0, color="black")

    elif "Trial energy" in y_axes[i]:
        if e_targ != None:
            axes[i,0].axhline(e_targ, color="red")

    elif "<V>" in y_axes[i]:
        if v_targ != None:
            axes[i,0].axhline(v_targ, color="red")

    # Histogram distribution
    d = [di for di in d if di > axes[i,0].get_ylim()[0] and di < axes[i,0].get_ylim()[1]]
    bins = int(min(100, int(len(d))/10))
    axes[i,1].hist(d, bins=bins, orientation="horizontal")
    axes[i,1].set_ylim(axes[i,0].get_ylim())
    axes[i,1].set_yticks([])

# Spread the subplots out, and show them
plt.subplots_adjust(wspace=0, hspace=0)
fig.align_ylabels(axes)
plt.show()
