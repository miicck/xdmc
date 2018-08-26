import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def title(n):
	titles = ["Population","Population change","Total energy","Potential energy","Kinetic energy"]
	if n < len(titles):
		return titles[n]
	return "No title"

lines = open(sys.argv[1]).read().split("\n")

data = []
for n in range(2,len(lines)):
	line = lines[n]
	dat = [n]
	try:
		for d in line.split(','):
			dat.append(float(d))
	except:
		continue

	if len(dat) > 1:

		bad_data = False
		for d in dat:
			if math.isnan(d):
				bad_data = True
				break
		if not bad_data:
			data.append(dat)

cols = zip(*data)

plots_x = np.ceil(np.sqrt(len(cols)-1))
plots_y = plots_x

for i in range(1,len(cols)):

	plt.subplot(plots_x,plots_y,i)
	plt.plot(cols[0],cols[i])
	plt.xlabel("Iteration")
	plt.ylabel(title(i-1))

plt.show()
