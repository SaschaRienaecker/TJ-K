""" Utility functions shared """
import numpy as np
from numpy import pi
import scipy.stats as scstats
import scipy.signal as scsignal
import matplotlib.pyplot as plt

dt = 1e-6 # sampling time in [s]

# For radial probes
R = np.arange(start=-32, step=5, stop=28.1, dtype=int) * 1e-3 # probe radii (r-a) [m]
NR = R.size
Z = np.array([-5, 0, 5]) * 1e-3 # probe vertical distance to midplane [m]
dZ = 5e-3 # vertical dist. between probes [m]


# For poloidal probes
dX = 8e-3 # distance between two successive probes [m]
dx_pol = 8e-3 # polidal distance between adjacent probes
X_theta = np.arange(64) * dx_pol # poloidal positions array
Theta = np.linspace(0, 2*np.pi, 64, endpoint=False) - np.pi
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

def fluctuations(a):
    return a - a.mean()

def statistical_properties(array):
    mean_array = np.mean(array)
    kurtosis_array = scstats.kurtosis(array)
    skew_array = scstats.skew(array)
    autocorr_array = scsignal.correlate(array, array)
    return mean_array, kurtosis_array, skew_array, autocorr_array

def find_nearest(a, val):
    imin = np.argmin(np.abs(a - val))
    return imin, a[imin]


def pdf_stat(Dat, shot='radial', itor=0):
    """Histogram of the density fluctuations (probability density) for the different radial/poloidal positions."""

    if shot=='radial':
        Isat = normalized(Dat[1, :, :])
    elif shot=='poloidal':
        iIsat = np.arange(1,64, step=2, dtype=int)
        Isat = Dat[iIsat, itor]

    X = Isat
    N = Isat.shape[0]

    for i in range(N):

        x = normalized(X[i])

        hist, bin_edges = np.histogram(x, bins=30, range=(-5,5))

        if i == 0:
            Hist = np.zeros((N, *hist.shape))
            Kurtosis = np.zeros(N)
            Skew = np.zeros(N)

        binw = bin_edges[1] - bin_edges[0]

        # normalize
        hist = hist / (np.sum(hist) * binw)
        Hist[i] = hist
        bin_centers = bin_edges[1:] - binw/2

        # pdf propoerties
        _, Kurtosis[i], Skew[i], _ = statistical_properties(x)

    return bin_centers, Hist, Skew, Kurtosis

def plot_pdf(bin_centers, Hist, Skew, Kurtosis, shot='radial', axs=None):

    if axs is None:
        fig, axs = plt.subplots(1,2, figsize=(7,3))

    [ax, ax2] = axs

    if shot=='radial':
        x = R * 1e3
        xlab = r'$r -a$ [mm]'
    elif shot=='poloidal':
        iIsat = np.arange(1,64, step=2, dtype=int)
        x = Theta[iIsat]
        xlab = r'poloidal angle $\theta$ [rad]'

    # for hist in Hist:
        # ax.plot(bin_centers, hist, alpha=0.2, color='black')
    imax_kurt = np.argmax(Kurtosis)
    ax.plot(bin_centers, Hist[imax_kurt], alpha=0.5, label='max. Kurt', color='red')
    imin_kurt = np.argmin(Kurtosis)
    ax.plot(bin_centers, Hist[imin_kurt], alpha=0.5, label='min. Kurt', color='blue')

    ax.plot(bin_centers, np.mean(Hist, axis=0), ls='--', color='black', label=r'$\theta$-average')
    ax.set_xlabel(r'$\tilde{I}$ [$\sigma_I$]')
    ax.set_ylabel(r'probab. density $\tilde{I}$')
    ax.legend(handlelength=1.0, loc='upper left')
    ax.set_xlim(left=-12)

    ax2.plot(x, Kurtosis, label=r'Kurt [$\tilde{I}$]')
    ax2.plot(x, Skew, label=r'Skew [$\tilde{I}$]')
    ax2.set_xlabel(xlab)
    ax2.legend(frameon=False)
    ax2.plot(x[imax_kurt], Kurtosis[imax_kurt], 'ro', ms=4)
    ax2.plot(x[imin_kurt], Kurtosis[imin_kurt], 'bo', ms=4)
    if shot=='poloidal':
        annot_poloidal_xaxis(ax2)
        # ax2.set_xticks([0, np.pi, 2 *np.pi])
        # ax2.set_xticklabels(['$0$', '$\pi$', '$2\pi$'])
    # plt.tight_layout()

def fluct_level(Dat, shot='poloidal', itor=0):
    """
    Returns the normalized density and potential fluctuations
    for all radii/poloidal positions.
    NOTE: Not implemented for radial shot yet.
    """
    if shot=='radial':
        Isat = Dat[1, :, :]
        phi  = Dat[0, :, :]

    elif shot=='poloidal':
        iPhi = np.arange(0,64, step=2, dtype=int)
        phi = Dat[iPhi, itor]
        iIsat = np.arange(1,64, step=2, dtype=int)
        Isat = Dat[iIsat, itor]


    Te = 9 # electron temperature in [eV] at probe tip position

    dn = np.std(Isat , axis=-1) / np.mean(Isat) # relative density fluctuations
    dphi = - np.std(phi, axis=-1) / np.mean(phi) # potential fluctuations in [V]

    return dn, dphi

def annot_poloidal_xaxis(ax, label=True, LFS=True, HFS=True):
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])

    if LFS:
        null_lab = '$0$\n(LFS)'
    else:
        null_lab = '$0$'

    if HFS:
        neg_pilab = '$-\pi$\n(HFS)'
        pos_pilab = '$\pi$\n(HFS)'
    else:
        neg_pilab = '$-\pi$'
        pos_pilab = '$\pi$'

    ax.set_xticklabels([neg_pilab, '$-\pi/2$', null_lab, '$\pi/2$', pos_pilab])

    # ax.set_xticklabels(['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'])
    if label:
        ax.set_xlabel(r'$\theta$ [rad]')

def plot_v_pol(CrossCorr_phi, CrossCorr_I, ax=None):

    labels = [r'$\tilde{\phi}$', r'$\tilde{I}$']
    iPhi = np.arange(0,64, step=2, dtype=int)

    if ax is None:
        fig, ax = plt.subplots()

    for i, CCorr in enumerate([CrossCorr_phi, CrossCorr_I]):
        tau_max = CCorr[:,2]
        max_corr = CCorr[:,3]

        # remove outliers
        #print(tau_max)
        itake = (np.abs(tau_max) > 1e-6) & (max_corr > 0.4)
        vfluct = - 2 * dx_pol / tau_max[itake]

        ax.plot(theta_array_OPA[iPhi][itake], vfluct * 1e-3, 'x-', label=labels[i])
        ax.set_ylabel(r'$v_\theta$ [km/s]')
        ax.axhline(0, ls='--', alpha=0.6, color='black')

    ax.legend()
    annot_poloidal_xaxis(ax)

# Blobs functions

# Gives the position of each value beyond alpha*sigma in absolute value
# Warning : normalized array needed
def blobholes(normalized_array, alpha = 2.3):
    indice = np.where(abs(normalized_array) >= alpha)[0]
    indiceblob = []
    indicehole = []
    for i in indice:
        if normalized_array[i]>0:
            indiceblob.append(i)
        else:
            indicehole.append(i)
    return indiceblob, indicehole

# Reduces the number of positions in order to take into account a time windowing
# This means that each blob/hole is referenced by a unique position (its center)
# The list of indices can be either blobs or holes
def blobhole_windowing(indices, window = 100):
    l = [[indices[0]]]
    for i in indices[1:]:
        if i - l[-1][0] < 2*window :
            l[-1].append(i)
        else:
            l.append([i])
    n = len(l)
    windowed_indices = np.zeros(n, dtype = int)
    for j in range(n):
        windowed_indices[j] = int((l[j][-1]+ l[j][0])//2)
    return windowed_indices

# Computes the average over every bolb/hole of the profile over time (mean evolution over time).
# If f_i is the time profile for blob nb i, it returns 1/N \sum f_i
window = 500
averaging_nb = 10

def blobholes_meanprof(normalized_array, windices, window = window):
    #averaging blob/hole
    mean_profile = np.zeros(2*window)

    for i in windices:
        mean_profile += normalized_array[ i- window : i + window]

    mean_profile /= len(windices)
    return mean_profile

# Computes the average over smaller groups of blobs (e.g. 10 by 10) of the profile
# if (i1_1, ... i1_10, i2_1,...i2_10,...,in_10) characterizes each blob, then
# it returns multiple profiles 1/10 \sum_j f_i1_j,..., 1/10 \sum_j f_in_10
def blobholes_local_meanprof(normalized_array, windices, averaging_nb = averaging_nb):
    #averaging blob/hole but only over a few pics
    nb_profiles = len(windices)//averaging_nb

    mean_profiles = np.zeros((nb_profiles, 2*window))

    compteur = 0
    m = 0
    while compteur < len(windices) and m < nb_profiles:
        somme = 0
        while somme <= averaging_nb:
            somme += 1
            i = windices[compteur]
            mean_profiles[m] += normalized_array[ i- window : i + window]/averaging_nb
        compteur += 1
        m +=1
        return mean_profiles
"""
These lines of command work provided you have the set of data datr_norm

datr_norm = normalized(Datr[0, 0])
Isat_norm = normalized(Datr[1, 0])

indiceblob, indicehole = blobholes(datr_norm)

windiceblob = blobhole_windowing(indiceblob)

phi_meanprof = blobholes_meanprof(datr_norm, windiceblob)
Isat_meanprof = blobholes_meanprof(Isat_norm, windiceblob)
phi_mean_profiles = blobholes_local_meanprof(datr_norm, windiceblob)
phi_mean_profiles.shape
"""

