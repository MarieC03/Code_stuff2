import read, amrplot
import numpy as np
import matplotlib as mpl  
import matplotlib.pyplot as plt
from matplotlib import rcParams  
from matplotlib import rc 
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib import rc 
import matplotlib.colors as mcolors
import imageio
import imageio.v2 as imageio  # Fix deprecation warning
from PIL import Image
import os
import scipy as scp
from scipy.integrate import trapezoid
from scipy.integrate import simpson
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from scipy.interpolate import Rbf
import argparse



rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif' 
rcParams['font.serif'] = ['Computer Modern']

rcParams['axes.labelsize']=12#8
#rcParams['axes.labelweight']=600
rcParams['axes.linewidth']=1
rcParams['lines.markersize']=6

rcParams['xtick.labelsize']=14#10
rcParams['ytick.labelsize']=14#10

plt.rcParams.update({
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "pgf.preamble": "\n".join([
        # r"\usepackage{url}",            # load additional packages
         r"\usepackage{unicode-math}",   # unicode math setup
         r"\setmainfont{Computer Modern}",  # serif font via preamble
    ])
})

FilenameTxt = "rho_b.maximum.asc"
data = []
data.append(np.genfromtxt(FilenameTxt))

#colorJ=["#5b0390", "#8632b9", "#b75ea9", "#c56785", "#e87f7f"]
colorJ=["#5b0390", "tab:blue"]

lw1=1.5
figw=6
figh=5.5

tmerg = 13.87 # 12.0

t_ms = data[0][:,1]*4.925490947 * 1.e-3 - tmerg
rho = data[0][:,2]

fig = plt.figure(figsize=(figw, figh), dpi=300)


plt.plot(t_ms,rho, label="BNS $M_{tot}=2.9~M_{sun}$",color=colorJ[1], linewidth=lw1)

#plt.xlabel(r'$t~[{M_{sun}}]$', fontsize=18)
plt.xlabel(r'$t~[\rm ms]$', fontsize=18)
plt.ylabel(r'$ \rho_{\rm max}$', fontsize=18)

#plt.tight_layout(rect=[0, 0, 1, 0.95])
filename = "Rho_b_max.png"
plt.savefig(filename,dpi=300,bbox_inches='tight',pad_inches=0.015,transparent=False)
