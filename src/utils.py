""" Utility functions shared """
import numpy as np
import scipy.stats as scstats
import scipy.signal as scsignal

dt = 1e-6 # sampling time in [s]

# For radial probes
R = np.arange(start=-32, step=5, stop=28.1, dtype=int) * 1e-3 # probe radii (r-a) [m]
NR = R.size
Z = np.array([-5, 0, 5]) * 1e-3 # probe vertical distance to midplane [m]
dZ = 5e-3 # vertical dist. between probes [m]






# For poloidal probes
dX = 8e-3 # distance between two successive probes [m]
l  = -1.5e-2 #distance to the separatrix
Bt = -72e-3 #mean magnetic field for now [T]

# Angle for each probe [rad]
theta_array_OPA = np.array([-3.1235280, -3.0137569, -2.9073485, -2.8051417, 
-2.7071847, -2.6126248, -2.5215145, -2.4330854, -2.3478083, -2.2648591,
-2.1842471, -2.1040693, -2.0234246, -1.9421344, -1.8614906, -1.7821543,
-1.7024737, -1.6212623, -1.5376113, -1.4496697, -1.3572263, -1.2596373,
-1.1571439, -1.0493933, -0.93729199, -0.82152736, -0.70318702, -0.58321263,
-0.46279019, -0.34215081, -0.22190346, -0.10223254, 0.017224783, 0.13673208,
0.25606317, 0.37513834, 0.49413465, 0.61282050, 0.73003303, 0.84526151, 0.95768375,
1.0665708, 1.1713279, 1.2715163, 1.3671756, 1.4586130, 1.5469808, 1.6322630,
1.7158621, 1.7989729, 1.8817460, 1.9647129, 2.0477027, 2.1307102, 2.2136460,
2.2962418, 2.3788249, 2.4615021, 2.5462539, 2.6346722, 2.7287933, 2.8293017,
2.9353781, 3.0454504])

# magnetic field at actual probe tip [T]
Bt_array_OPA    = np.array([81.413653,  81.554947, 81.688886, 81.785425,
81.805867, 81.677173, 81.418237, 80.978830, 80.410903, 79.665974, 78.762967,
77.573484, 76.043523, 74.128822, 72.002351, 69.823963, 67.703469, 65.680512,
63.775266, 61.973155, 60.285373, 58.697925, 57.225119, 55.836454, 54.551215,
53.374512, 52.301156, 51.361428, 50.541115, 49.880512, 49.379216, 49.074743,
48.995043, 49.098658, 49.399573, 49.909240, 49.413465, 51.404843, 52.347126,
53.423241, 54.605184, 55.900020, 57.305778, 58.811793, 60.431098, 62.159817,
64.012493, 65.972022, 68.033278, 70.161485, 72.300145, 74.362368, 76.264515,
77.940651, 79.336531, 80.374150, 81.064533,  81.410054, 81.525579, 81.490495,
81.405933, 81.336058, 81.291581, 81.331378])*10**(-3)






def normalized(a):
    return (a - a.mean() ) / a.std()

<<<<<<< HEAD
def fluctuations(a):
    return a - a.mean()

=======
>>>>>>> 2c3df1a974a6bff37bf385e116df7bddc5c64190
def statistical_properties(array):
    mean_array = np.mean(array)
    kurtosis_array = scstats.kurtosis(array)
    skew_array = scstats.skew(array)
    autocorr_array = scsignal.correlate(array, array)
<<<<<<< HEAD

    return mean_array, kurtosis_array, skew_array, autocorr_array

# def correlation_properties()

=======
    
    return mean_array, kurtosis_array, skew_array, autocorr_array
>>>>>>> 2c3df1a974a6bff37bf385e116df7bddc5c64190
