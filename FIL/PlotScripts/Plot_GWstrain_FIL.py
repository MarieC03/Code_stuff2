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

gomega=0.005

T=2.*np.pi/gomega
dist=500

tMsolTomsec=0.00493
#fig_aux, ax_aux = plt.subplots(figsize=(10,5),dpi=128)
# orbital velocity from initial data
iters_to_merger_largest=-20
iters_to_merger_smallest=100000

sim_to_color={    "BNS_DD2_MNS1_2.7_MNS2_1.2_CHINS1_0.585_CHINS2_0.7_ECC_RED_0_dx_0.144_RK3_CFL0p2" : "sienna"}

sim_to_strain={ "BNS_DD2_MNS1_2.7_MNS2_1.2_CHINS1_0.585_CHINS2_0.7_ECC_RED_0_dx_0.144_RK3_CFL0p2"  : None          }
fig, ax = plt.subplots(figsize=(10,5), dpi=128)


counter=0
for sim in ["bhns_dd2_d50"]:
    sim_dir= SimDir(dir_name+sim)
    
    # loading strain data of dominant l=m=2 mode
    strain=sim_dir.gws[dist].get_strain_lm(2,2,T,0.001,window_function='tukey',trim_ends=False)
    sim_to_tmer_ms[sim]=strain.x_at_abs_maximum_y()*tMsolTomsec
    print(sim, "t_mer = ", sim_to_tmer_ms[sim], "ms")
    
    # strain.align_at_maximum()
    # extracting ++ and xx polarization
    # rescaling to 40 Mpc
    
    # by default, the strain in kuibit is given as r*(h_+ - i h_x),
    # rescaling of the strain is : r_1 / r_2 
    
    hpp=strain.real().y  /  (40. * MpcToLMsol) * 1e21
    hxx=-strain.imag().y /  (40. * MpcToLMsol) * 1e21

    
    # hpp=strain.real().y  /  (40. * MpcToLMsol) 
    # hxx=-strain.imag().y /  (40. * MpcToLMsol) 

    
    # locating maximum of GW amplitude
    amp=list(np.sqrt(hpp[:]**2+hxx[:]**2))
    max_amp=max(amp)
    i_max=amp.index(max_amp)
    time=strain.real().t
    tmerge=time[i_max]
    sim_to_strain[sim]=strain
    
    ax.plot(strain.t*tMsolTomsec, hpp)


    counter+=1
    #ax3.set_xlim(-10,1)
plt.show()
#plt.plot(strain)