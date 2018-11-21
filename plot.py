import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Labels for the plots
y_axes = ["Population", "Potential energy (Hartree)"]

# Take a rolling average with a window of sample_length
def rolling_av(data, sample_length):
	ret = []
	val = data[0]
	r = 1.0/sample_length
	for d in data:
		ret.append(val)
		val = val * (1-r) + d*r; 
	return ret

# Read in our data from the "out" file
data = []
for line in open("out").read().split("\n"):
	if len(line.strip()) > 0:
		data.append([float(i) for i in line.split(",")])

# Plot each dataseries on it's own subplot
data = zip(*data)
for i, d in enumerate(data):
	plt.subplot(len(data),1,i+1)
	plt.plot(d)
	plt.plot(rolling_av(d,len(d)/100))
	plt.xlabel("Iteration")
	plt.ylabel(y_axes[i])

# Spread the subplots out, and show them
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
