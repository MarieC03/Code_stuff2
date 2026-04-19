import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as signal
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.default']


filepath = "../"

simname = "./output/single_perturb.log"

#x, y = np.loadtxt(open(filepath+simname,'rt').readlines()[:-1], unpack=True ,dtype='float',usecols=(1,3),skiprows=1)
file = np.loadtxt('./cowling_IG_640_1d_do50.log', skiprows=2)
#file = np.loadtxt('./output/1-12_ppm_640.log', skiprows=2)

file2 = np.loadtxt('./cfc_IG_640_1d_do50.log', skiprows=2)

y = file[:,3]
x = file[:,1]
alp = file[:,4]

y2 = file2[:,3]
x2 = file2[:,1]
alp2 = file2[:,4]


y_fft = file2[:,3]
x_fft = file2[:,1]

y_fft = y_fft/file2[1,3] - 1
x_fft = x_fft/ 2.03001708E05

#print(y)
#print('#############')
y = y/file[1,3]-1
#y = y/ max(y)
x = x / 2.03001708E05 # code unit to sec 
alp = alp/file[1,4] -1

y2 = y2/file2[1,3]-1
x2 = x2 / 2.03001708E05 # code unit to sec 
alp2 = alp2/file2[1,4] -1


###fft
df_TD = pd.DataFrame({'t' : x_fft, 'rho' : y_fft})
df_TD = df_TD[['t','rho']]
df_TD.to_csv("rho_c_TD.csv", index = False, sep = '\t')

t_min = 0.0
t_max = float(df_TD.nlargest(1,'t')['t'])
#t_max = 1.0E-2
print (t_max)
time_step = 1.0/2.0**12 * t_max
print (time_step)
time = np.arange(0.0, float(t_max), float(time_step))

window = signal.tukey(len(time), alpha = 0.2)
#window = signal.tukey(len(time), alpha = 0.5)
#window = np.hanning(len(time))
#window = np.blackman(len(time))
#window = np.kaiser(len(time), beta = 14)

print(np.array(df_TD['t']))
print(df_TD['rho'])
data_t = interp1d(df_TD['t'], df_TD['rho'], kind = 'cubic')
#data_t = interp1d(np.array(df_TD['t']), np.array(df_TD['rho']), kind = 'cubic')
data_temp = data_t(time)# * 1.0e5
data_max = df_TD['rho'][1]
print (data_max)
data_temp -= data_max
data_temp -= np.mean(data_temp)
data_temp = data_temp * window

data_f = fft(data_temp) * 2.0/time.size
freq = np.linspace(start=0.0, stop=1.0/(2.0 * time_step), num=int(time.size/2) )
data_psd = np.abs(data_f)**2

df_FD = pd.DataFrame( {'f': freq, 'fft': np.abs(data_f[:time.size//2]), 'psd': data_psd[:time.size//2]} )
df_FD = df_FD[['f','fft','psd']]
df_FD.to_csv('rho_c_FD.csv', index = False, sep = '\t')

#fig, (ax1, ax2) = plt.subplots(2,constrained_layout=True,sharex=True)
##fig, (ax1, ax2, ax3) = plt.subplots(3,constrained_layout=True)
#
#ax1.plot(x,y/1e-3, label='BU0-CW')
#ax1.plot(x2,y2/1e-3, label='BU0-DY', color='red')
#ax1.set_ylabel('$\\rho_{c}(t)/ \\rho_{c}(0) - 1$')
#ax1.set_title('$\\times 10^{-3}$', x=0.05, y=0.97, fontsize=12)
#ax1.legend()
##ax2.plot(x2,alp2)
#ax2.set_ylabel('$\\alpha_{c}(t)/ \\alpha_{c}(0) - 1$')
#ax2.set_xlabel('$t~[{\mathrm{ms}}]$')
#ax2.plot(x2,alp2/1e-4, color='red')
#ax2.set_title('$\\times 10^{-4}$', x=0.05, y=0.97, fontsize=12)
#
#ax1.get_xaxis().set_visible(False)
#ax1.sharex(ax2)
#
##ax2.plt.xlim(0,0.01)
#plt.xlabel('$t$ $(s)$')
#plt.xlim(0,0.01)  #10ms
#plt.ylim(-0.005,0.005)
#plt.xlabel('$f$')
#plt.ylabel('$psd$')
plt.xlabel('$f~[{\mathrm{Hz}}]$')
plt.ylabel('$\mathrm{log}_{10}\mathrm{(PSD)}$')
plt.axvline(x=1458, color='red',label='F (H. Dimmelmeier 2006)',linestyle='dashed')
plt.axvline(x=3971, color='blue',label='H1 (H. Dimmelmeier 2006)',linestyle='dashed')
#plt.axvline(x=, color='purple',label='H2')
#plt.axvline(x=8140, color='green', label='H3')
#plt.axvline(x=2690, color='red',label='F')
#plt.axvline(x=4547, color='blue',label='H1')
#plt.axvline(x=6340, color='purple',label='H2')
#plt.axvline(x=8140, color='green', label='H3')
#plt.text(x, .5, 'ls220beta_radialmode')
plt.xlim(500,7000)
plt.ylim(-15,-8)
#plt.plot(x,y, label='x1y1')
#plt.plot(x2,y2, label='x2y2')
plt.legend(loc='upper right', prop={'size':12})
plt.plot(freq,np.log10(data_psd[:time.size//2]))
#plt.plot(freq,data_psd[:time.size//2])
#plt.show()
plt.savefig("../Figures/new_IG_1D_cfc_vs_cowling_fft_10ms_50domain.pdf", bbox_inches="tight")

