import matplotlib.pyplot as plt
import sys
import numpy as np

lines = open(sys.argv[1]).read().split('\n')
data = []
for line in lines:
	try:
		data.append(float(line))
	except:
		continue

plt.subplot(211)
plt.xlim([-4,4])
plt.ylim([0,50000])
plt.hist(data,bins=1000)
plt.subplot(212)
plt.xlim([-4,4])
plt.ylim([0,50000])
plt.hist(np.random.normal(size=10000000),bins=1000)
plt.show()
