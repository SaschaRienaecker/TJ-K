""" Utility functions shared """
import numpy as np

dt = 1e-6 # sampling time in [s]
R = np.arange(start=-32, step=5, stop=28.1, dtype=int) * 1e-3 # probe radii (r-a) [m]
NR = R.size
Z = np.array([-5, 0, 5]) * 1e-3 # probe vertical distance to midplane [m]
dZ = 5e-3 # vertical dist. between probes [m]

def normalized(a):
    return (a - a.mean() ) / a.std()
