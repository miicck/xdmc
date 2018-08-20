import sys
import math
import numpy as np
import matplotlib.pyplot as plt

lines = open(sys.argv[1]).read().split("\n")
titles = lines[0].split(',')
equils = lines[1].split(',')

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

	if equils[i-1] == "1":
		second_half = cols[i][len(cols[i])/2:]
		av = np.mean(second_half)
		sd = np.std(second_half)

		av_arr = []
		for j in range(0,len(cols[0])):
			av_arr.append(av)
		av_arr = np.array(av_arr)

		plt.plot(cols[0],av_arr,label="Equilibrated average = " + str(av))
		plt.legend(loc="best")

	plt.plot(cols[0],cols[i])
	plt.ylabel(titles[i-1])
	plt.xlabel("Iteration")

plt.show()
