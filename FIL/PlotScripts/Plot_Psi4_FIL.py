
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

dir_name = "/mnt/raarchive/cassing/FIL_runs/"
gomega=0.009 # <- FUKA orbital omega from .info file 

# orbital period
T=2.*np.pi/gomega
# distance in Msun of detector of GW extraction
dist=500
tMsolTomsec=0.00493
fig, ax = plt.subplots(figsize=(10,5), dpi=128)


for sim in ["BNS_SFHo_FUKA"]:
    sim_dir= SimDir(dir_name+sim)
    print(sim)
    
    # loading strain data of dominant l=m=2 mode

    psi4=sim_dir.gws[dist].get_psi4_lm(2,2)#,0.1,window_function='tukey',trim_ends=True)
    
    # extracting ++ and xx polarization
    psi4pp=psi4.real().y
    psi4xx=-psi4.imag().y
    
    # locating maximum of GW amplitude
    amp=list(np.sqrt(psi4pp[:]**2+psi4xx[:]**2))
    # max_amp=max(amp)
    # i_max=amp.index(max_amp)
    time=psi4.real().t

    

    ax.plot(psi4.t*tMsolTomsec, psi4.y)

ax.set_ylabel(r"$\Psi_{4}^{22}$", fontsize=15)

ax.set_xlabel(r"$t~[\rm{ms}]$",fontsize=15)
ax.set_xlim(0,50)
ax.set_xlim(0,2)
#ax.set_ylim(-1e-6,1e-6)
ax.set_ylim(-1e-6,1e-6)

#plt.subplots_adjust(wspace=0.02)
outpath = dir_name+sim+"/Plots/Psi4"
fig.savefig(outpath+".pdf", bbox_inches="tight", dpi=180)
fig.savefig(outpath+".png", bbox_inches="tight", dpi=180)