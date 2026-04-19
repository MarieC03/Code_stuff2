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


d=read.load(0,file='./cylind_dd2bns_20ms_230125_128x128lv5do400_prim_weno5',type='vtu')
d2=read.load(15,file='./cylind_dd2bns_20ms_230125_128x128lv5do400_prim_weno5',type='vtu')
d3=read.load(30,file='./cylind_dd2bns_20ms_230125_128x128lv5do400_prim_weno5',type='vtu')
d4=read.load(707,file='./cylind_dd2bns_50ms_230125_128x128lv5do400_prim',type='vtu')

time_of_slice1 = d.time
time_of_slice2 = d2.time
time_of_slice3 = d3.time
time_of_slice4 = d4.time
print(time_of_slice1, 'time_of_slice1')
print(time_of_slice2, 'time_of_slice2')
print(time_of_slice3, 'time_of_slice3')
print(time_of_slice4, 'time_of_slice4')

# t0 time slice
file =  np.loadtxt('./cylind_dd2bns_20ms_230125_128x128lv5do400_prim_weno5_d2_x+0.00D+00_n0000.csv', delimiter=',', skiprows=1)
file2 = np.loadtxt('./cylind_dd2bns_20ms_230125_128x128lv5do400_prim_weno5_d2_x+0.00D+00_n0000.csv', delimiter=',', skiprows=1)

# Next time slice
file4 = np.loadtxt('./cylind_dd2bns_20ms_230125_128x128lv5do400_prim_weno5_d2_x+0.00D+00_n0015.csv', delimiter=',', skiprows=1)

#file5 = np.loadtxt('./cylind_dd2bns_50ms_230125_128x128lv5do400_prim_weno5_d2_x+0.00D+00_n0010.csv', delimiter=',', skiprows=1)
#file6 = np.loadtxt('./cylind_dd2bns_50ms_230125_256x256lv5do400_prim_weno5_d2_x+0.00D+00_n0016.csv', delimiter=',', skiprows=1)

file5 = np.loadtxt('./cylind_dd2bns_20ms_230125_128x128lv5do400_prim_weno5_d2_x+0.00D+00_n0030.csv', delimiter=',', skiprows=1)
file6 = np.loadtxt('./cylind_dd2bns_50ms_230125_128x128lv5do400_prim_d2_x+0.00D+00_n0707.csv', delimiter=',', skiprows=1)
#column, 1:x, 2:y, bphi: 11, 15: rho,   18: u3, 21:lfac
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

y2 = file2[:,1] 
x2 = file2[:,0]
rho2 = file2[:,14]
Bphi2 = file2[:,10]
u3_2   = file2[:,17]
lfac2 = file2[:,20]
vphi2 = u3_2/lfac2
ye2   = file2[:,13]
temp2   = file2[:,12]
psi2 = file2[:,34]
r2 = x2


y4 = file4[:,1] 
x4 = file4[:,0]
rho4 = file4[:,14]
Bphi4 = file4[:,10]
u3_4   = file4[:,17]
lfac4 = file4[:,20]
vphi4 = u3_4/lfac4
psi4 = file4[:,34]
ye4   = file4[:,13]
temp4   = file4[:,12]
r4 = x4

y5 = file5[:,1] 
x5 = file5[:,0]
rho5 = file5[:,14]
Bphi5 = file5[:,10]
u3_5   = file5[:,17]
lfac5 = file5[:,20]
vphi5 = u3_5/lfac5
ye5   = file5[:,13]
psi5 = file5[:,34]
temp5   = file5[:,12]
r5 = x5

y6 = file6[:,1]
x6 = file6[:,0]
rho6 = file6[:,14]
Bphi6 = file6[:,10]
u3_6   = file6[:,17]
lfac6 = file6[:,20]
vphi6 = u3_6/lfac6
ye6   = file6[:,13]
psi6 = file6[:,34]
temp6   = file6[:,12]
r6 = x6


# change unit
r  = r/0.677140812
r2 = r2/0.677140812
r5 = r5/0.677140812
r6 = r6/0.677140812
r4 = r4/0.677140812

#plt.ylabel('$\\frac{\\rho}{\\rho_{0}}, \\frac{v^{\\phi}}{v^{\\phi}_{0}}, \\frac{B^{\\phi}}{B^{\\phi}_{0}}$')

max_rho_t0_1stpost    = 1
max_logrho_t0_1stpost = 1
max_bphi_t0_1stpost   = 1 
max_vphi_t0_1stpost   = 1 

max_rho_t0_2ndpost    = 1
max_logrho_t0_2ndpost = 1
max_bphi_t0_2ndpost   = 1
max_vphi_t0_2ndpost   = 1

fig, axes = plt.subplots(nrows=3, ncols=3, sharex=True)
#fig, axes = plt.subplots(nrows=4, ncols=3, sharex=True, sharey='row')
#fig, axes = plt.subplots(nrows=4, ncols=3, constrained_layout=True, sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.1,hspace=0)
ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9= axes.flatten()

ax1.plot(r,rho, '.', markerfacecolor='none', markersize=5, color='blue')
#ax1.plot(r,rho, color='black', linestyle='dashed')
ax1.set_yscale('log')
ax1.set_title('$\\rho$')
#ax1.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

ax2.plot(r,temp, '.', markerfacecolor='none', markersize=5, color='blue')
#ax2.plot(r,temp, color='black', linestyle='dashed')
ax2.set_title('$T~[{{\mathrm{MeV}}}]$')

ax3.plot(r,psi, '.', markerfacecolor='none', label='$\\mathtt{BHAC}$', markersize=5, color='blue')
#ax3.plot(r,psi, color='black', linestyle='dashed')
ax3.legend(loc='upper right',prop={'size': 12})
ax3.set_title('$\\psi$')

ax4.plot(r4,rho4, '.', markerfacecolor='none', markersize=5, color='blue')
#ax4.plot(r,rho,  label='$t=0$',color='black', linestyle='dashed')
ax4.set_yscale('log')

ax5.plot(r4,temp4, '.', markerfacecolor='none', markersize=5, color='blue')
#ax5.plot(r,temp, color='black', linestyle='dashed')

ax6.plot(r4,psi4, '.', markerfacecolor='none', markersize=5, color='blue')
#ax6.plot(r,psi, color='black', linestyle='dashed')

ax7.plot(r5,rho5, '.', markerfacecolor='none', markersize=5, color='blue')
#ax7.plot(r6,rho6, '.', markerfacecolor='none', markersize=5, color='orange')
ax7.set_yscale('log')
#ax7.plot(r,rho, color='black', linestyle='dashed')

ax8.plot(r5,temp5, '.', markerfacecolor='none', markersize=5, color='blue')
#ax8.plot(r6,temp6, '.', markerfacecolor='none', markersize=5, color='orange')
#ax8.plot(r,temp, color='black', linestyle='dashed')

ax9.plot(r5,psi5, '.', markerfacecolor='none', markersize=5, color='blue')
#ax9.plot(r6,psi6, '.', markerfacecolor='none', markersize=5, color='orange')
#ax9.plot(r,psi, color='black', linestyle='dashed')

ax7.set(xlabel='$r~[{\mathrm{km}}]$')
ax8.set(xlabel='$r~[{\mathrm{km}}]$')
ax9.set(xlabel='$r~[{\mathrm{km}}]$')

ax1.set(ylabel='$t = t_{\\mathrm{HO}}~[{\mathrm{ms}}]$')
ax4.set(ylabel='$t = t_{\\mathrm{HO}}+10~[{\mathrm{ms}}]$')
ax7.set(ylabel='$t = t_{\\mathrm{HO}}+30~[{\mathrm{ms}}]$')
# Set the ticks and ticklabels for all axes
plt.setp(axes, xticks=[0, 5, 10,15, 20,25], xticklabels=[0,5, 10, 15,20,25])

ax1.set_ylim(2e-7, 2e-3)
ax4.set_ylim(2e-7, 2e-3)
ax7.set_ylim(2e-7, 2e-3)

#ax2.set_ylim(0,35)
#ax5.set_ylim(0,35)
#ax8.set_ylim(0,35)

ax3.set_ylim(1,1.39)
ax6.set_ylim(1,1.39)
ax9.set_ylim(1,1.39)
[ax1.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]
[ax2.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]
[ax3.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]
[ax4.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]
[ax5.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]
[ax6.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]
[ax7.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]
[ax8.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]
[ax9.axvline(x=i, color='lightgrey', alpha=0.3) for i in [5,10,15,20,25]]

[ax1.axhline(y=i, color='lightgrey', alpha=0.3) for i in [1e-7,1e-6,1e-5,1e-4,1e-3]]
[ax4.axhline(y=i, color='lightgrey', alpha=0.3) for i in [1e-7,1e-6,1e-5,1e-4,1e-3]]
[ax7.axhline(y=i, color='lightgrey', alpha=0.3) for i in [1e-7,1e-6,1e-5,1e-4,1e-3]]

[ax2.axhline(y=i, color='lightgrey', alpha=0.3) for i in [10,20,30]]
[ax5.axhline(y=i, color='lightgrey', alpha=0.3) for i in [10,20,30]]
[ax8.axhline(y=i, color='lightgrey', alpha=0.3) for i in [10,20,30]]

[ax3.axhline(y=i, color='lightgrey', alpha=0.3) for i in [1.1,1.2,1.3]]
[ax6.axhline(y=i, color='lightgrey', alpha=0.3) for i in [1.1,1.2,1.3]]
[ax9.axhline(y=i, color='lightgrey', alpha=0.3) for i in [1.1,1.2,1.3]]



plt.xlim(0,30) #for RNS
#plt.plot(r4,Bphi4/Bphi, linestyle='dashed')
plt.show()
#plt.savefig("../Figures/DD2BNS_from20+10ms_lv5_128do400_1dslices.pdf")

