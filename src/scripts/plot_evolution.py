import matplotlib.pyplot as plt
from parser import parse_evolution

# Read in the evolution data
y_axes, data = parse_evolution()
sq_size = int(len(data)**0.5)+1

# Plot each dataseries on its own subplot
for i, d in enumerate(data):
	plt.subplot(sq_size,sq_size,i+1)
	plt.plot(d)
	plt.xlabel("Iteration")
	plt.ylabel(y_axes[i])

# Spread the subplots out, and show them
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
