import matplotlib.pyplot as plt
from parser import parse_evolution

# Take a rolling average with a window of sample_length
def rolling_av(data, sample_length):
	ret = []
	val = data[0]
	r = 1.0/sample_length
	for d in data:
		ret.append(val)
		val = val * (1-r) + d*r;
	return ret

# Read in the evolution data
y_axes, data = parse_evolution()
sq_size = int(len(data)**0.5)+1

# Plot each dataseries on its own subplot
for i, d in enumerate(data):
	plt.subplot(sq_size,sq_size,i+1)
	plt.plot(d)
	plt.plot(rolling_av(d,len(d)/100))
	plt.xlabel("Iteration")
	plt.ylabel(y_axes[i])

# Spread the subplots out, and show them
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
