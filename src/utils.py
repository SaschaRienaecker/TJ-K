""" Utility functions shared """
import numpy as np
import scipy.stats as scstats
import scipy.signal as scsignal

dt = 1e-6 # sampling time in [s]
R = np.arange(start=-32, step=5, stop=28.1, dtype=int) * 1e-3 # probe radii (r-a) [m]
NR = R.size
Z = np.array([-5, 0, 5]) * 1e-3 # probe vertical distance to midplane [m]
dZ = 5e-3 # vertical dist. between probes [m]

def normalized(a):
    return (a - a.mean() ) / a.std()

def statistical_properties(array):
    mean_array = np.mean(array)
    kurtosis_array = scstats.kurtosis(array)
    skew_array = scstats.skew(array)
    autocorr_array = scsignal.correlate(array, array)

    return mean_array, kurtosis_array, skew_array, autocorr_array
