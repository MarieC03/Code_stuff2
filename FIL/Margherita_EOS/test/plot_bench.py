import numpy as np
import matplotlib.pyplot as plt
from numpy import log10
import sys

input_string =sys.argv[1]

name_string =input_string
if name_string.endswith('.dat'):
        name_string = name_string[:-4]

x, y, z, r = np.loadtxt(input_string, delimiter='\t', usecols=(0,1,2,3), unpack=True)

jet = plt.cm.get_cmap('jet') 

fig, ax1 = plt.subplots(1, 1)

ax1.set_title(name_string)

# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
CS=ax1.tripcolor(np.log10(x),np.log10(y),np.log10(z), vmin=-16, vmax=-1, cmap=jet,shading='gouraud')

ax1.set_xlabel(r'$\log_{10}\ \rho$')
ax1.set_ylabel(r'$\log_{10}\ T$')
ax1.set_autoscaley_on(True)
ax1.set_autoscalex_on(True)

cb= plt.colorbar(CS, extend='both')
cb.set_label(label=r'$\mathrm{L1\ error\ on\ primitives}$')

xmin = np.amin(x)
ymin = np.amin(y)
xmax = np.amax(x)
ymax = np.amax(y)

ax1.set_xlim([log10(xmin),log10(xmax)])
ax1.set_ylim([log10(ymin),log10(ymax)])

#fig.show()
plt.savefig(name_string + '.pdf')


### Lets plot the residual


jet = plt.cm.get_cmap('jet') 

fig, ax1 = plt.subplots(1, 1)

ax1.set_title(name_string+ '_residual')

# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
CS=ax1.tripcolor(np.log10(x),np.log10(y),np.log10(r), vmin=-16, vmax=-1, cmap=jet,shading='gouraud')

ax1.set_xlabel(r'$\log_{10}\ \rho$')
ax1.set_ylabel(r'$\log_{10}\ T$')
ax1.set_autoscaley_on(True)
ax1.set_autoscalex_on(True)

cb= plt.colorbar(CS, extend='both')
cb.set_label(label=r'$\mathrm{residual}$')

xmin = np.amin(x)
ymin = np.amin(y)
xmax = np.amax(x)
ymax = np.amax(y)

ax1.set_xlim([log10(xmin),log10(xmax)])
ax1.set_ylim([log10(ymin),log10(ymax)])

#fig.show()
plt.savefig(name_string + '_residual.pdf')

