import matplotlib.pyplot as plt
from parser import parse_evolution
import sys

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

# Plot each dataseries on its own subplot
for i, d in enumerate(data):
    plt.subplot(y_size,x_size,i+1)
    plt.plot(d)
    plt.xlabel("Iteration")
    plt.ylabel(y_axes[i])

    if y_axes[i] == "Cancelled weight":
        plt.yscale("log")

    elif "Average weight" in y_axes[i]:
        plt.axhline(0, color="black")

    elif "Trial energy" in y_axes[i]:
        if e_targ != None:
            plt.axhline(e_targ, color="red")

# Spread the subplots out, and show them
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
