import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Labels for the plots
y_axes = ["Population", "Trial energy (Hartree)"]
how_to_combine = [
	lambda(p) : np.sum(p, axis=0), # Sum populations across processes
	lambda(e) : np.mean(e, axis=0) # Average energies across processes
]

# Take a rolling average with a window of sample_length
def rolling_av(data, sample_length):
	ret = []
	val = data[0]
	r = 1.0/sample_length
	for d in data:
		ret.append(val)
		val = val * (1-r) + d*r; 
	return ret

# Read in the evolution from each process
all_data = []
for f in os.listdir("."):
	if not f.startswith("evolution"): continue

	# Read in our data from the evolution file
	data = []
	for line in open(f).read().split("\n"):
		if len(line.strip()) > 0:
			data.append([float(i) for i in line.split(",")])
	data = zip(*data)
	all_data.append(data)

# Combine the data from each process
data = []
for i in range(0, len(all_data[0])):
	data.append(how_to_combine[i]([d[i] for d in all_data]))

# Plot each dataseries on its own subplot
for i, d in enumerate(data):
	plt.subplot(len(data),1,i+1)
	plt.plot(d)
	plt.plot(rolling_av(d,len(d)/100))
	plt.xlabel("Iteration")
	plt.ylabel(y_axes[i])

# Spread the subplots out, and show them
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
