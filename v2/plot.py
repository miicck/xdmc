import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def rolling_av(data, sample_length):
	ret = []
	val = data[0]
	r = 1.0/sample_length
	for d in data:
		ret.append(val)
		val = val * (1-r) + d*r; 
	return ret

data = []
for line in open("out").read().split("\n"):
	if len(line.strip()) > 0:
		data.append([float(i) for i in line.split(",")])

data = zip(*data)
plot_size = np.ceil(np.sqrt(len(data)))
for i, d in enumerate(data):
	plt.subplot(plot_size,plot_size,i+1)
	plt.plot(d)
	plt.plot(rolling_av(d,1000))

plt.show()
