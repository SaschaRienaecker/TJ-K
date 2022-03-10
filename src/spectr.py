""" Spectral analysis """
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch, csd, correlate, coherence
from scipy.signal import correlation_lags # Note: requires a recent version of SciPy
from utils import R, dt, NR, Z, dZ, normalized, theta_array_OPA, annot_poloidal_xaxis

dt = 1e-6 # sampling time in [s]
nperseg = 4 * 1024
dx_pol = 8e-3 # polidal distance between adjacent probes
Theta = np.linspace(0, 2*np.pi, 64, endpoint=False) - np.pi
X_theta = np.arange(64) * dx_pol


def get_tau_corr(tau, corr, threshold=np.exp(-1)):
    """Return tau at 1/e width (left and right value) of the correlation maximum."""
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

def Corr_profile_poloidal(Dat, quantity='phi', itor=0, threshold=np.exp(-1)):
    """
    Returns for each poloidal position the left and right 1/e time delay,
    the time corresponding to the auto-correlation maximum and the correlation value at that maximum.
    """
    iPhi = np.arange(0,64, step=2, dtype=int)
    Phi = Dat[iPhi, itor]
    iIsat = np.arange(1,64, step=2, dtype=int)
    Isat = Dat[iIsat, itor]

    Q = Phi if quantity=='phi' else Isat
    Npol = Q.shape[0]

    Corr = np.zeros((Npol,4))
    for i in range(Npol):
        dat = normalized(Q[i])
        corr = correlate(dat, dat, method='fft') / dat.size
        lags = correlation_lags(dat.size, dat.size)
        tau = lags * dt

        Corr[i,:2] = get_tau_corr(tau, corr, threshold)
        Corr[i, 2] = tau[np.argmax(corr)]
        Corr[i, 3] = np.max(corr)

    return Corr

def adjacent_Corr(Dat, quantity='phi', threshold=np.exp(-1), itor=0):
    """
    Returns for each radial position the left and right 1/e time delay,
    the time corresponding to the correlation maximum and the correlation value at that maximum.
    """

    iPhi = np.arange(0,64, step=2, dtype=int)
    Phi = Dat[iPhi]
    iIsat = np.arange(1,64, step=2, dtype=int)
    Isat = Dat[iIsat]

    E = -(Phi - np.roll(Phi, -1, axis=0)) / (2 * dx_pol)


    if quantity=='phi':
        Q = Phi
    elif quantity=='I':
        Q = Isat
    elif quantity=='E_I':
        Q = E
    elif quantity=='phi_I':
        Q = Phi

    Corr = np.zeros((Q.shape[0],4))

    for i in range(Q.shape[0]):

        dat1 = normalized(Q[i,itor])

        if quantity=='E_I':
            dat2 = normalized(Isat[i, itor])
        elif quantity=='phi_I':
            dat2 = normalized(Isat[i, itor])
        else:
            if i==Q.shape[0]-1:
                i2 = 0
            else:
                i2 = i+1

            dat2 = normalized(Q[i2, itor])

        corr = correlate(dat1, dat2, method='fft') / dat1.size
        lags = correlation_lags(dat1.size, dat2.size)
        tau = lags * dt

        Corr[i,:2] = get_tau_corr(tau, corr, threshold)
        Corr[i, 2] = tau[np.argmax(corr)]
        Corr[i, 3] = np.max(corr)

    return Corr

def plot_spec(Spec, f, ax=None, cbar=True, angle=False, shot='radial', lognorm=False, vmin=1e-3, vmax=1):
    """Just a convenience function for plotting"""
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    Z = Spec / Spec.max()
    if angle:
        Z = Spec

    if shot=='radial':
        ext = [0, f.max()/1e3, R.min(), R.max()]
        ax.set_ylabel('$(r-a)$ [mm]')
        ax.set_xlabel('$f$ [kHz]')
        ax.set_xlim(right=20)
    elif shot=='poloidal':
        ext = [theta_array_OPA.min(), theta_array_OPA.max(), 0, f.max()/1e3]
        annot_poloidal_xaxis(ax)
        ax.set_ylabel('$f$ [kHz]')
        ax.set_ylim(top=20)
        Z = Z.T

    if lognorm:
        from matplotlib.colors import LogNorm
        norm=LogNorm(vmin, vmax)
        im = ax.imshow(Z, aspect='auto', extent=ext, origin='lower', norm=norm)
    else:
        im = ax.imshow(Z, aspect='auto', extent=ext, origin='lower')
    return im


    if cbar:
        lab = 'PSD [a.u.]' if not angle else 'angle [rad]'
        fig.colorbar(im, ax=ax, label=lab)

def get_kspec(Dat, quantity='phi', itor=0, it_step=10):
    """

    """
    from numpy.fft import fft, fftfreq, fftshift
    dx_pol = 8e-3

    iPhi = np.arange(0,64, step=2, dtype=int)
    Phi = Dat[iPhi]
    iIsat = np.arange(1,64, step=2, dtype=int)
    Isat = Dat[iIsat]

    Isat_m = np.mean(Isat[:,itor,:], axis=-1)
    itake = Isat_m > 1
    Isat[itake] = 0


    Q = Phi if quantity=='phi' else Isat
    # Q = normalized(Q)



    # select some time frames
    it_select = np.arange(start=0, stop=Q.shape[-1], step=it_step)

    for i, it in enumerate(it_select):
        kspec = fft(Q[:, itor, it])

        if i == 0:
            kth = 2 * np.pi * np.fft.fftfreq(Q.shape[0], d=2 * dx_pol)
            kth =  np.fft.fftshift(kth)
            Kspec = np.zeros((it_select.size, *kspec.shape), dtype='complex')

        kspec = np.fft.fftshift(kspec)
        Kspec[i] = kspec
    return kth, Kspec
