
import numpy as np
import pandas as pd
from pathlib import Path
import scipy.stats as scstats
import scipy.signal as scsignal

datap = Path('../data')

def read_rad_prof(rad_position, probe_nr):
    """
    rad_position : index for the radial position of the triple probe. Must be an integer between  0 and 12.
    probe_nr: index of the probe in the triple probe setup. Must be either, 0,1 or 2.
    """
    if probe_nr == 0 or probe_nr == 2:
        ext = ".ufl"
    if probe_nr == 1:
        ext = ".isa"
    if probe_nr != 0 and probe_nr != 1 and probe_nr !=2:
        print("Probe number must be 0,1 or 2 depending on the probe we use : 0 and 2 are measuring potential and 1 is measuring ion saturation current.")
        ext = 'NaN'

    if ext != 'NaN':
        # path= datap / str("20100216#006709/TJ-K20100216#006709pos00"+"%02d" % (rad_position,)+"_0"+str(probe_nr)+ext )
        path = datap / "20100216#006709/TJ-K20100216#006709pos00{:02d}_0{}{}".format(rad_position,probe_nr,ext)

        #Load measurements in a panda dataframe.
        data_pandas=pd.read_csv(path,skiprows=20,engine='python',header=None,delim_whitespace=True,skipfooter=2)
        #Convert it to numpy array
        data_array = data_pandas.values
        #Reshape it from (174763, 6) to  (1048576)
        data_array=np.reshape(data_array,[np.shape(data_array)[0]*np.shape(data_array)[1]])
        #Remove nans (they come from number of values in each data columns:not the same length)
        data_nan=~np.isnan(data_array)
        data_array=data_array[data_nan]
        return(data_array)

def extract_to_binary(shot='radial'):

    if shot=='radial':
        Np = 3    # number of probes
        NR = 13   # number or radii


        for ip in np.arange(Np):
            for iR in np.arange(NR):

                dat = read_rad_prof(iR,ip)

                if ip==0 and iR==0:
                    Dat = np.zeros((Np, NR, dat.size)) # placeholder

                Dat[ip, iR] = dat

        p_binary = datap / '20100216#006709/dat.npy'
        np.save(p_binary, Dat)
    elif shot == 'poloidal':
        Np = 64    # number of probes
        NR = 2   # number or radii


        for ip in np.arange(Np):
            for iR in np.arange(NR):

                dat = read_rad_prof(iR,ip)

                if ip==0 and iR==0:
                    Dat = np.zeros((Np, NR, dat.size)) # placeholder

                Dat[ip, iR] = dat

        p_binary = datap / '20100920#007192/dat.npy'
        np.save(p_binary, Dat)

def statistical_properties(array):
    mean_array = np.mean(array)
    kurtosis_array = scstats.kurtosis(array)
    skew_array = scstats.skew(array)
    autocorr_array = scsignal.correlate(array, array)
    
    return mean_array, kurtosis_array, skew_array, autocorr_array

# def correlation_properties()
    






















