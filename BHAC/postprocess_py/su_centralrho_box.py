import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as signal
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.default']

filepath = "../"

simname = "./output/single_perturb.log"

#x, y = np.loadtxt(open(filepath+simname,'rt').readlines()[:-1], unpack=True ,dtype='float',usecols=(1,3),skiprows=1)
#file = np.loadtxt('./fluidID.csv', skiprows=3)
#file2 = np.loadtxt('./code_paper_data/2d_ls220_specialrefine4_R_15_r_core.log', skiprows=2)

file  = np.loadtxt('./su2d_pert.log', skiprows=2)
#file2 = np.loadtxt('./code_paper_data/tov_ls220.log', skiprows=2)
file2 = np.loadtxt('./SU2D_pert_rho.scalars.asc')
#file3 = np.loadtxt('./code_paper_data/mirgation.txt')
#PM_noinit_readmetricfully_128x24_lv3_do100
y = file[:,2]

print(max(y))


#file2 = file
#file2 = np.loadtxt('./output/similarres_FIL_SU2d_idealgas.log', skiprows=2)
y = file[288:,3]
x = file[288:,1]

y4 = file2[:218,3]
x4 = file2[:218,1]

y2 = file2[:,3]
x2 = file2[:,1]

#y3 = file3[:,1]
#x3 = file3[:,0]

#print(y)
#print('#############')
y = y/file[1,3] 
x = x / 2.03001708E05 + 4.58e-5# shift
#for i in range(len(x)):
#      if x[i-1] < 0.000162:
#         y[i-1] = 0
#x = x / 2.03001708E05 # code unit to sec 

y4 = y4/file2[1,3]
x4 = x4 / 2.03001708E05 # code unit to sec 
y2 = y2/file2[1,3]
x2 = x2 / 2.03001708E05 # code unit to sec 

#y3 = y3/file3[1,1]
#x3 = x3 / 2.03001708E05 # code unit to sec 

x_detail = x
y_detail = y
x2_detail = x2
y2_detail = y2


plt.xlabel('$t~[{\mathrm{s}}]$')
plt.ylabel('$\\rho_{c}(t)/ \\rho_{c}(0)$')
#plt.xlim(0,0.05)
plt.xlim(0,0.01)  #10ms
#plt.ylim(-0.005,0.005)
plt.axhline(y = 1.346e-3/8e-3, color = 'black', linestyle = 'dotted')

#plt.axhline(1.346e-3/8e-3 - 1)
plt.plot(x,y, label='$\mathtt{BHAC}$', color='blue')
plt.plot(x4,y4, color='blue')
plt.plot(x2,y2, label='$\mathtt{FIL}$',linestyle='dashed', color='green')
#plt.plot(x3,y3, label='CC2009',linestyle='dashed')
plt.legend()
x1_pt, y1_pt = 0.001110, 0.411
x2_pt, y2_pt = 0.00677, 0.522
plt.plot([x1_pt, x2_pt], [y1_pt, y2_pt], color='grey', linestyle='--')
#plt.annotate("",
#            xy=(x1_pt, y1_pt), xycoords='data', 
#            xytext=(x2_pt, y2_pt), textcoords='data',
#            arrowprops=dict(arrowstyle="-"
#                            #shrinkA=5, shrinkB=5,
#                            #patchA=None, patchB=None,
#                            #connectionstyle="arc3,rad=0.",
#                            ),
#            )

x3_pt, y3_pt = 0.00288, 1.004
plt.plot([x1_pt, x3_pt], [y1_pt, y3_pt], color='grey', linestyle='--')
#plt.annotate("",
#            xy=(x1_pt, y1_pt), xycoords='data', 
#            xytext=(x3_pt, y3_pt), textcoords='data', 
#            arrowprops=dict(arrowstyle="-"
#                            #shrinkA=5, shrinkB=5,
#                            #patchA=None, patchB=None,
#                            #connectionstyle="arc3,rad=0.",
#                            ),
#            )


 # location for the zoomed portion 
sub_axes = plt.axes([.35, .50, .30, .35]) 

# plot the zoomed portion
sub_axes.plot(x_detail, y_detail, c = 'blue')
sub_axes.plot(x2_detail, y2_detail, c = 'green')
sub_axes.set_xlim(0.001,0.0013)
sub_axes.set_ylim(0.3,0.44)
sub_axes.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)


#plt.show()
plt.savefig('../Figures/spherical_su2d_IG__thr1d-12_fac1_do60_640x32.pdf',bbox_inches="tight")

