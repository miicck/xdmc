import numpy as np
import matplotlib.pyplot as plt
import sys

# Read pre-binned data in
with open(sys.argv[1]) as f:
    par, wfn = f.read().split("\n")
    particles = [[float(x) for x in s.split(",") if len(x) > 0] for s in par.split(";") if len(s) > 0]
    data = [d for d in wfn.split(";") if d != ""]
data = [[float(x) for x in d.split(",") if len(x) > 0] for d in data]
data = np.array(data).T
abs_max = max(-np.min(data), np.max(data))

RES   = 100
RANGE = 4.0

#plt.imshow(data, extent=(-RANGE,RANGE,-RANGE,RANGE))
xs = np.linspace(-RANGE, RANGE, 100)
lev = np.linspace(-abs_max, abs_max, 11)
plt.contour(xs,xs,data,levels=lev)
plt.gca().set_aspect(1.0)
plt.scatter([p[0] for p in particles],[p[1] for p in particles],color="red")
plt.show()
