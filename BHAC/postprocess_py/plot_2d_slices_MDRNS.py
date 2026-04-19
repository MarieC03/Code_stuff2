import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as signal
import read, amrplot
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.default']

d=read.load(0,file='./MRNS_2d_poly_64x32_lv3_oldmetinit',type='vtu')
d2=read.load(10,file='./MRNS_2d_poly_64x32_lv3_oldmetinit',type='vtu')
time_of_slice1 = d.time
time_of_slice2 = d2.time
print(time_of_slice1, 'time_of_slice1')
print(time_of_slice2, 'time_of_slice2')

fig, axes = plt.subplots(nrows=1, ncols=4, sharex=True)
ax1, ax2, ax3, ax4= axes.flatten()
p1 = amrplot.polyplot(np.log10(d.rho),d,xrange=[0,20],yrange=[0,20], clear=False, axis=ax1, fig=fig) 
p2 = amrplot.polyplot(d.u3/d.lfac,d,xrange=[0,20],yrange=[0,20], clear=False, axis=ax2, fig=fig) 
p3 = amrplot.polyplot(np.log10(d.rho),d2,xrange=[0,20],yrange=[0,20], clear=False, axis=ax3, fig=fig) 
p4 = amrplot.polyplot(d.u3/d.lfac,d2,xrange=[0,20],yrange=[0,20], clear=False, axis=ax4, fig=fig) 

ax1.set(xlabel='$x~[{\mathrm{km}}]$')
ax2.set(xlabel='$x~[{\mathrm{km}}]$')
ax3.set(xlabel='$x~[{\mathrm{km}}]$')
ax4.set(xlabel='$x~[{\mathrm{km}}]$')
ax1.set(ylabel='$z~[{\mathrm{km}}]$')
ax2.set(ylabel='$z~[{\mathrm{km}}]$')
ax3.set(ylabel='$z~[{\mathrm{km}}]$')
ax4.set(ylabel='$z~[{\mathrm{km}}]$')
plt.subplots_adjust(wspace=0.1,hspace=0)


plt.show()
#plt.savefig("./MDRNS_t0_2dslice_logrho_t0.pdf", bbox_inches="tight")
#p2.savefig("../Figures/MDRNS_t0_2dslice_v3_t10.pdf", bbox_inches="tight")
#p3.savefig("../Figures/MDRNS_t0_2dslice_logrho_t0.pdf", bbox_inches="tight")
#p4.savefig("../Figures/MDRNS_t0_2dslice_v3_t10.pdf", bbox_inches="tight")

