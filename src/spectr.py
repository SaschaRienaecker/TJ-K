""" Spectral analysis """
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch, csd, correlate, coherence
from scipy.signal import correlation_lags # Note: requires a recent version of SciPy
from utils import R, dt, NR, Z, dZ, normalized

dt = 1e-6 # sampling time in [s]
nperseg = 4 * 1024

def get_tau_corr(tau, corr, threshold=np.exp(-1)):

    ip = tau > 0
    tau_p = tau[ip]
    itau_c_p = np.argmin(np.abs(corr[ip] - corr.max() * threshold ))
    tau_c_p = tau_p[itau_c_p]

    i_n = tau <= 0
    tau_n = tau[i_n]
    itau_c_n = np.argmin(np.abs(corr[i_n] - corr.max() * threshold ))
    tau_c_n = tau_n[itau_c_n]

    return tau_c_n, tau_c_p

def Corr_profile(Dat, ip1, ip2, threshold=np.exp(-1)):
    """
    Returns for each radial position the left and right 1/e time delay,
    the time corresponding to the correlation maximum and the correlation value at that maximum.
    """

    Corr = np.zeros((R.size,4))
    for iR in range(R.size):
        dat1 = normalized(Dat[ip1, iR])
        dat2 = normalized(Dat[ip2, iR])

        corr = correlate(dat1, dat2, method='fft') / dat1.size
        lags = correlation_lags(dat1.size, dat2.size)
        tau = lags * dt

        Corr[iR,:2] = get_tau_corr(tau, corr, threshold)
        Corr[iR, 2] = tau[np.argmax(corr)]
        Corr[iR, 3] = np.max(corr)

    return Corr

def plot_spec(Spec, f, ax=None, cbar=True, angle=False):

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    Z = Spec / Spec.max()
    if angle:
        Z = Spec
    im = ax.imshow(Z, aspect='auto', extent=[0, f.max()/1e3, R.min(), R.max()], origin='lower')

    ax.set_xlabel('$f$ [kHz]')
    ax.set_ylabel('$(r-a)$ [mm]')
    ax.set_xlim(right=20)

    if cbar:
        lab = 'PSD [a.u.]' if not angle else 'angle [rad]'
        fig.colorbar(im, ax=ax, label=lab)
