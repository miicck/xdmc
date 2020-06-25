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
import sys
import numpy as np
from math import floor, log10
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from parser import parse_evolution
    
round_to_n = lambda x, n: round(x, -int(floor(log10(abs(x)))) + (n - 1)) if abs(x) > 10e-9 else 0

# Parse command line args
e_targ = None
for arg in sys.argv:
    if arg.startswith("-e="):
        val = arg.split("=")[-1]
        if   val == "lithium":   e_targ = -7.47807
        elif val == "beryllium": e_targ = -14.667353
        elif val == "boron":     e_targ = −24.65386608
        else: e_targ = float(val)
        break

v_targ = None
for arg in sys.argv:
    if arg.startswith("-v="):
        v_targ = float(arg.split("=")[-1])
        break

# Read in the evolution data
y_axes, data = parse_evolution()

i_include=[]
for a in sys.argv:
    try:
        i_include.append(int(a))
    except:
        pass

if len(i_include) > 0:
    new_data  = []
    new_yaxes = []
    for i in range(0, len(data)):
        if not i in i_include: continue
        new_data.append(data[i])
        new_yaxes.append(y_axes[i])
    y_axes = new_yaxes
    data   = new_data

print("Plotting {0} variables.".format(len(data)))

log_scales = ["Cancelled weight"]

fig, axes = plt.subplots(len(data), 2, gridspec_kw = {"width_ratios":[8,1]})

# Plot each dataseries on its own subplot
for i, d in enumerate(data):

    sigma_range = 3

    # Get statistics for dataseries
    half = int(len(d)/2)
    std  = round_to_n(np.std(d[half:]),  4)
    mean = round_to_n(np.mean(d[half:]), 4)
    err  = round_to_n(std / (float(half)**0.5), 4)

    # Plot dataseries
    axes[i,0].plot(d)
    axes[i,0].set_xlabel("Iteration")
    axes[i,0].set_ylabel(y_axes[i]+"\n{0} +/- {1}".format(mean,err))
    axes[i,0].axhline(mean, color="green", linestyle=":")

    # Set the axis scales
    if y_axes[i] in log_scales:
        axes[i,0].set_yscale("log")
        axes[i,1].set_yscale("log")

    elif "norm" in sys.argv:
        delta = sigma_range*std
        if delta < 10e-5: delta = 10e-5
        axes[i,0].set_ylim([mean - delta, mean + delta])

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
    bins = int(min(100, int(len(d)/10 + 1)))
    axes[i,1].hist(d, bins=bins, orientation="horizontal")
    axes[i,1].set_ylim(axes[i,0].get_ylim())
    axes[i,1].set_yticks([])

# Spread the subplots out, and show them
plt.subplots_adjust(wspace=0, hspace=0)
try:
    fig.align_ylabels(axes)
except:
    print("Could not align y labels!")
plt.show()


