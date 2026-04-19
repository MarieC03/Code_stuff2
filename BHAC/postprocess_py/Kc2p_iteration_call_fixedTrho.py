import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as signal
import matplotlib as mpl
#from colorspacious import cspace_converter
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d
from matplotlib import rc
from matplotlib.cm import ScalarMappable
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.default']

filepath = "../"

file   = np.loadtxt('./c2p_iteration_fixedrhoT.log',skiprows=1)
#file   = np.loadtxt('./c2p_rootfinding_call.log',skiprows=1)
#file   = np.loadtxt('./c2p_EOS_call.log',skiprows=1)
#file  = np.loadtxt('./relative_diff_yevs_rhoT_lowestYe.log',skiprows=1)
file2 = np.loadtxt('./iW_of_table.log')
y1 = file[:,2]  #n pmag to p
x1 = file[:,1]  #n W
z1 = file[:,0]
z1_per_cell = z1/64
W_array = file2[:]

nW = 200
n_pmag_to_p = 200
#Z = [[0 for x in range(ntemp)] for y in range(nW)] 

Z = np.zeros((nW,n_pmag_to_p))

count = 0 
for i in range(nW):
    for j in range(n_pmag_to_p):
        count = count+1
#        Z[i][j] = z1[count-1] - 105.66
        Z[j][i] = z1_per_cell[count-1]
#        if Z[j][i] > 1000:
#           Z[j][i] = 0
#        if Z[i][j] < 100:
#           Z[i][j] = 0.0
#        if Z[i][j] > 0.5:
#            Z[i][j] = 0.0
#@        print(Z[i][j],x1[count-1],y1[j], count)
#@        print(x1[count], y1[j], count);


#print('max of Z')
#print(max(map(max, Z))) 

#print(Z[255][162])  ! Z is correct

Y = np.zeros(nW)
X = np.zeros(n_pmag_to_p) 

for i in range(n_pmag_to_p):
    X[i] = np.log10(y1[i])
    #if X[i] > 11.5:
    #   X[i] = 11.5
    #Y[i] = y1[i]
#    print(Y[i])
        
#Y = np.log10(W_array)
Y = np.log10(W_array-1)
#X = W_array

#print(max_z)
#Z = np.log10(Z)
#print(X)

# maximum value of Z array
print('max')
print(max(map(max, Z)) )
print('min')
print(min(map(min, Z)) )
print('average')
print(sum(z1_per_cell)/len(z1_per_cell))


vmax = 50
vmin = 0
#vmin = min(map(min, Z))
#vmax = max(map(max, Z)) 

#Z = np.log10(Z)
levels = np.linspace(vmin, vmax, 51)

fig,ax=plt.subplots(1,1)
cp = ax.contourf(Y, X, Z, levels=levels)
#cp = ax.contourf(Y, X, Z, levels=levels, extend='max')

#cp.cmap.set_over('white')
#cp.cmap.set_over('b')
#set_ylim(bottom=8.1, top=11.5)
#bounds=[0,1]

#fig.colorbar(
#   ScalarMappable(norm=cp.norm, cmap=cp.cmap),
#   ticks=range(vmin, vmax, 1)
#)
#fig.colorbar(cp) # Add a colorbar to a plot
cbar = plt.colorbar(cp) # Add a colorbar to a plot
#cbar.ax.set_yticklabels(['0','5','10','15','20','25','30','35','40','45'])
cbar.set_label('# of iterations', rotation=90)
#cp.cmap.set_over('black')


#plt.colorbar(boundaries=np.linspace(0,0.1,1)) 
ax.set_title('')
#ax.set_title('Kastaun\'s primitive recovery with tabulated EOS')
ax.set_xlabel('$\mathrm{log}_{10}(W - 1)$')
ax.set_ylabel('$\mathrm{log}_{10}(b^2/2p)$')
plt.show()
