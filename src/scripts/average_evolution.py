from parser import parse_evolution 
import numpy as np

yaxes, data = parse_evolution()
for yn, d in zip(yaxes, data):
        print yn, np.mean(d[len(d)/2:])
