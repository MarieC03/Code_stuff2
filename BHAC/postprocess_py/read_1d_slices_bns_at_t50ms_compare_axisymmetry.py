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


d=read.load(0,file='./cylind_dd2bns_50ms_230125_128x128lv5do400_prim_weno5',type='vtu')

# t0 time slice
file =  np.loadtxt('./cylind_dd2bns_20ms_230125_128x128lv5do400_prim_weno5_d2_x+0.00D+00_n0000.csv', delimiter=',', skiprows=1)

y = file[:,1] 
x = file[:,0]
rho = file[:,14]
Bphi = file[:,10]
u3   = file[:,17]
lfac = file[:,20]
vphi = u3/lfac
ye   = file[:,13]
temp   = file[:,12]
psi = file[:,34]
r = x

# change unit
r  = r/0.677140812

plt.plot(r,np.log10(rho),'.', markerfacecolor='none', markersize=5, color='blue',label='$\mathtt{BHAC}$ at $\\bar{t}=50~\mathrm{ms}$')
plt.xlabel('$r~[{\mathrm{km}}]$')
plt.ylabel('$\mathrm{log}_{10}(\\rho)$')
plt.ylim(-8, -2.5)
plt.xlim(0,50) #for RNS
plt.legend(loc='best')
#plt.plot(r4,Bphi4/Bphi, linestyle='dashed')
plt.show()
#plt.savefig("../Figures/DD2BNS_from20+10ms_lv5_128do400_1dslices.pdf")

