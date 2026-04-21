import os
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

from kuibit.simdir import SimDir, load_SimDir
from kuibit.grid_data import UniformGrid, UniformGridData
from kuibit.timeseries import TimeSeries, remove_duplicated_iters
from kuibit import grid_data as gd
from kuibit import grid_data_utils as gdu
from kuibit import series


N_increase_f_resolution=1
savgol_window_size=11
savgol_window_order=5

# def spectrogram(strain):
# take from Sam's code
def get_spectrogram(strain_in, 
    tminus, tplus, 
    fminus=0.6, fplus=5, tmerge=None,
    f22=None, f21=None, fmax=None, f3=None, f3times=None,
    axs=None, LabelSize=LabelSize, NFFT=80, fgwmax_t=None, 
    vmin = -40, vmax = 5, **kwargs):
   
    from scipy import signal
    from scipy.interpolate import interp1d

    strain = strain_in.copy()

    tmerge = strain.time_at_maximum() if tmerge is None else tmerge
    strain.align_at_maximum()
    
    hp = strain.y  * 10**22. / 40. / MpcToLMsol
    ht = strain.t

    # Sampling rate
    ttot = ht[-1] - ht[0]
    #Fs = len(ht) / ttot
    Fs = 1. / (ht[1] - ht[0])
    lht = len(ht)
    # Samples per window
#    NFFT= 20#int(len(ht)/)
    #NFFT= 2**10
    # Number of overlapping samples
    #Nover = NFFT / 2
    Nover = int(NFFT*0.98)
    f, t, Sxx = signal.spectrogram(
        hp, 
        fs=Fs, 
        nperseg=NFFT, 
        noverlap = Nover, 
        scaling='spectrum', 
        window=signal.get_window("blackman", NFFT), 
        detrend='linear')
    tfinal = t * 1e3 - tmerge - tminus / 2.
    trange = np.where((tfinal >= tminus) & (tfinal <= tplus))
    frange = np.where((f/1e3 >= fminus) & (f/1e3 <= fplus))

    #print("fspec dt: {}".format(tfinal[-1] - tfinal[0]), tfinal[0], tfinal[-1], tmerge/1e3)
    Sxx = Sxx[frange,:][0]
    f = f[frange]
    
    # here we limit the spectrogram to the relevant times
    newSxx = list()
    for val in Sxx:
        newSxx.append(val[trange])
    Sxx = np.array(newSxx)
    logSxx = 10.*np.log10(Sxx)
    tfinal = tfinal[trange]

    dt = tminus - tfinal[0]


    return f, logSxx


# effective strain as done in PostCactus gw_utils
def h_eff(strain, polarization="both", truncate_inspiral=False, end=None):
    from scipy import signal
    from kuibit.frequencyseries import FrequencySeries
    ts = strain.regular_resampled()
    if truncate_inspiral:
    # sets t=0 to strain max    
        ts.align_at_maximum()
    # clips all data up to t=0
        ts.clip(init=0)

    if end != None:
        ts.clip(end=end)

    # We extract the plus and cross components of the strain
    # Note: for an unaltered Kuibit strain, this will in reality be
    # r * hp, r * hc respectively
    hp = ts.real()
    hc = -ts.imag()
  
    # Then, we take the Fourier transform.
    # By extracting hp and hc seperately we obtain
    # the real signal only in the positive frequency space
    hp_fft = hp.to_FrequencySeries()
    hc_fft = hc.to_FrequencySeries()

    # Finally, we can compute the effective strain amplitude 
    # power spectral density akin to [Eq (8-9) in 1604.00246].
    # For a full discussion, see https://github.com/Sbozzolo/kuibit/pull/27
    h_eff = hp_fft
    
    print("Doing polarization: ", polarization)

    if(polarization=="+"):
        h_eff.fft = h_eff.fft = h_eff.f * np.sqrt(
                 (hp_fft.amp**2. ) / 2.0
                )

    elif(polarization=="x"):
        h_eff.fft = h_eff.fft = h_eff.f * np.sqrt(
                 (hc_fft.amp**2. ) / 2.0
                )
    elif(polarization=="both"):
        h_eff.fft = h_eff.f * np.sqrt(
        (hp_fft.amp**2. + hc_fft.amp**2.) / 2.0
          )
    
    else:
        print("YOU MUST CHOOSE POLARIZATION!")

    return h_eff  

######### GW call
CU_time=4.925553195445955e-06
           for sim in sim_to_strain.keys():
            sim_to_strain[sim] = sim_to_strain[sim].time_unit_changed(CU_time, inverse=True).copy()
            fGW=sim_to_strain[sim].phase_frequency().abs()

                    # find the f_mer frequency:
        print("Find f_mer:")
        index_t_mer=np.where(fGW.t>=0)[0][0]
        print(fGW.t[index_t_mer], fGW.t[0], fGW.t[-1])
        f_mer=fGW.y[index_t_mer]/1000
        print("f_mer: ", f_mer, " kHz" )
     ```