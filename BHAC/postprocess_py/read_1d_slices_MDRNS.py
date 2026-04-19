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

#file  = np.loadtxt('./code_paper_data/MRNS_2d_poly_d2_x+0.00D+00_n0000.csv', skiprows=1)
file =  np.loadtxt('./MRNS_2d_poly_64x32_lv3_oldmetinit_stag_d2_x+0.79D+00_n0000.csv', skiprows=1)
file2 = np.loadtxt('./MRNS_2d_poly_64x32_lv3_oldmetinit_stag_d2_x+0.16D+01_n0000.csv', skiprows=1)
#file4 = np.loadtxt('./code_paper_data/MRNS_2d_poly_d2_x+0.00D+00_n0840.csv', skiprows=1)
file4 = np.loadtxt('./MRNS_2d_poly_64x32_lv3_oldmetinit_stag_d2_x+0.79D+00_n0010.csv', skiprows=1)
file5 = np.loadtxt('./MRNS_2d_poly_64x32_lv3_oldmetinit_stag_d2_x+0.16D+01_n0010.csv', skiprows=1)

#column, 1:x, 2:y, bphi: 11, 15: rho,   18: u3, 21:lfac
y = file[:,1] 
x = file[:,0]
rho = file[:,14]
Bphi = file[:,10]
u3   = file[:,17]
lfac = file[:,20]
vphi = u3/lfac
r = np.sqrt(x**2 + y**2)

y2 = file2[:,1] 
x2 = file2[:,0]
rho2 = file2[:,14]
Bphi2 = file2[:,10]
u3_2   = file2[:,17]
lfac2 = file2[:,20]
vphi2 = u3_2/lfac2
r2 = np.sqrt(x2**2 + y2**2)

y4 = file4[:,1] 
x4 = file4[:,0]
rho4 = file4[:,14]
Bphi4 = file4[:,10]
u3_4   = file4[:,17]
lfac4 = file4[:,20]
vphi4 = u3_4/lfac4
r4 = np.sqrt(x4**2 + y4**2)

y5 = file5[:,1] 
x5 = file5[:,0]
rho5 = file5[:,14]
Bphi5 = file5[:,10]
u3_5   = file5[:,17]
lfac5 = file5[:,20]
vphi5 = u3_5/lfac5
r5 = np.sqrt(x5**2 + y5**2)


# change unit
r = r/0.677140812
r2 = r2/0.677140812
r5 = r5/0.677140812
r4 = r4/0.677140812


#plt.ylabel('$\\frac{\\rho}{\\rho_{0}}, \\frac{v^{\\phi}}{v^{\\phi}_{0}}, \\frac{B^{\\phi}}{B^{\\phi}_{0}}$')


max_rho_t0_45deg = max(rho)
max_logrho_t0_45deg = np.log10(max(rho))
max_bphi_t0_45deg = max(Bphi) 
max_vphi_t0_45deg = max(vphi) 

max_rho_t0_90deg = max(rho2)
max_logrho_t0_90deg = np.log10(max(rho2))
max_bphi_t0_90deg = max(Bphi2) 
max_vphi_t0_90deg = max(vphi2) 


#plt.axhline(1.346e-3/8e-3 - 1)
#plt.plot(r4,rho4/max_rho_t0_45deg, '.', markerfacecolor='none', color='blue', label='$\\frac{\\rho}{\\rho_{\mathrm{max}}(0)}$ at $\\theta=\\pi/4$')
#plt.plot(r,rho/max_rho_t0_45deg, color='black', linestyle='dashed')
#plt.plot(r4,Bphi4/max_bphi_t0_45deg, '.', markerfacecolor='none', color='green', label='$ \\frac{B^{\\phi}}{B^{\\phi}_{\mathrm{max}}(0)}$ at $\\theta=\\pi/4$')
#plt.plot(r,Bphi/max_bphi_t0_45deg, color='black', linestyle='dashed')
#plt.plot(r4,vphi4/max_vphi_t0_45deg, '.', markerfacecolor='none', color='red', label='$\\frac{v^{\\phi}}{v^{\\phi}_{\mathrm{max}}(0)}$ at $\\theta=\\pi/4$')
#plt.plot(r,vphi/max_vphi_t0_45deg, color='black', linestyle='dashed')

#plt.plot(r5,rho5/max_rho_t0_90deg, '^', markerfacecolor='none', color='blue', label='$\\frac{\\rho}{\\rho_{\mathrm{max}}(0)}$ at $\\theta=\\pi/2$')
#plt.plot(r2,rho2/max_rho_t0_90deg, color='black', linestyle='dashdot')
#plt.plot(r5,Bphi5/max_bphi_t0_90deg, '^', markerfacecolor='none', color='green', label='$ \\frac{B^{\\phi}}{B^{\\phi}_{\mathrm{max}}(0)}$ at $\\theta=\\pi/2$')
#plt.plot(r2,Bphi2/max_bphi_t0_90deg, color='black', linestyle='dashdot')
#plt.plot(r5,vphi5/max_vphi_t0_90deg, '^', markerfacecolor='none', color='red', label='$\\frac{v^{\\phi}}{v^{\\phi}_{\mathrm{max}}(0)}$ at $\\theta=\\pi/2$')
#plt.plot(r2,vphi2/max_vphi_t0_90deg, color='black', linestyle='dashdot')

fig, (ax1, ax2) = plt.subplots(2,constrained_layout=True, sharex=True)
#fig.suptitle('Vertically stacked subplots')

#ax1.plot(x2,y2, label='BU0 Dynamical')
#ax1.set_ylabel('$\\rho_{c}(t)/ \\rho_{c}(0) - 1$')

ax1.plot(r4,rho4/max_rho_t0_45deg, '.', markerfacecolor='none', color='blue', label='$\\frac{\\rho}{\\rho_{\mathrm{max}}(0)}$ at $\\theta=\\pi/4$')
ax1.plot(r,rho/max_rho_t0_45deg, color='black', linestyle='dashed')
ax1.plot(r4,Bphi4/max_bphi_t0_45deg, '.', markerfacecolor='none', color='green', label='$ \\frac{B^{\\phi}}{B^{\\phi}_{\mathrm{max}}(0)}$ at $\\theta=\\pi/4$')
ax1.plot(r,Bphi/max_bphi_t0_45deg, color='black', linestyle='dashed')
ax1.plot(r4,vphi4/max_vphi_t0_45deg, '.', markerfacecolor='none', color='red', label='$\\frac{v^{\\phi}}{v^{\\phi}_{\mathrm{max}}(0)}$ at $\\theta=\\pi/4$')
ax1.plot(r,vphi/max_vphi_t0_45deg, color='black', linestyle='dashed')

ax1.legend()
#ax2.plot(x2,alp2)
#ax2.set_ylabel('$\\alpha_{c}(t)/ \\alpha_{c}(0) - 1$')

ax2.plot(r5,rho5/max_rho_t0_90deg, '^', markerfacecolor='none', color='blue', label='$\\frac{\\rho}{\\rho_{\mathrm{max}}(0)}$ at $\\theta=\\pi/2$')
ax2.plot(r2,rho2/max_rho_t0_90deg, color='black', linestyle='dashdot')
ax2.plot(r5,Bphi5/max_bphi_t0_90deg, '^', markerfacecolor='none', color='green', label='$ \\frac{B^{\\phi}}{B^{\\phi}_{\mathrm{max}}(0)}$ at $\\theta=\\pi/2$')
ax2.plot(r2,Bphi2/max_bphi_t0_90deg, color='black', linestyle='dashdot')
ax2.plot(r5,vphi5/max_vphi_t0_90deg, '^', markerfacecolor='none', color='red', label='$\\frac{v^{\\phi}}{v^{\\phi}_{\mathrm{max}}(0)}$ at $\\theta=\\pi/2$')
ax2.plot(r2,vphi2/max_vphi_t0_90deg, color='black', linestyle='dashdot')

ax2.legend()
#ax2.plt.xlim(0,0.01)
plt.xlabel('$r~[{\mathrm{km}}]$')
plt.xlim(0,20)
plt.ylim(0,1)

#plt.plot(r4,Bphi4/Bphi, linestyle='dashed')
#plt.show()
plt.savefig("../Figures/MDRNS_poly_64x32_lv3_both.pdf", bbox_inches="tight")

