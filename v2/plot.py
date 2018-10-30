import matplotlib.pyplot as plt
import numpy as np
import sys

data = []
for line in open(sys.argv[1]).read().split("\n"):
	try:
		data.append([float(i) for i in line.split(",")])
	except:
		continue

data = zip(*data)
plot_size = np.ceil(np.sqrt(len(data)))
for i, d in enumerate(data):
	plt.subplot(plot_size,plot_size,i+1)
	plt.plot(d)

plt.show()
