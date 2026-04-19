import os
import read
import numpy as np
import matplotlib.pyplot as plt
import amrplot
import matplotlib as mpl
mpl.use('Agg')

#######################

# The aim of this script is comparing two different slices of the BHAC 'vtu' result 'with the same data structure'

#######################


# parameters you have to set to adapt to your problems
data_dir = "/lustre/hpe/ws10/ws10.3/ws/xfpjiang-jiangws/Simulation/sim_bhac/GRHD/BU8_nstag/nstag/" # directory that store the data
out_fig_dir = "../../result/fig_buff/" # directory that store the output figures
slicedir = 3                           # slice direction, the same as that in the amrvac.par file
name_tag = "output_d{}_x+0.00D+00_n".format(slicedir) # name of the output files
slice1, slice2 = 0, 180                # slices you want to compare, result would be slice2-slice1
interested_pars = ['rho', ]      # variables you want to compare, for each variable we genereate a figure, need more variables just need to add into this list
log_pars = ['nothing']                 # if a variable is in this list, then log scale is used
x_cut = 1.e-6                          # not needed unless you want logscale of some variables, if so, you need to read carefully the function strange_log
# end of parameter setting, if you are familliar with the amrplot, you can set its parameters youself below


import kuibit.unitconv as kbunit
myunit = kbunit.geom_umass_msun(1)

# Abs Values below x_cut are all neglected, read this function carefully then you will know how to interpret the colorbar of the output figure
def strange_log(x, x_cut=1.e-5):
	sign_x, abs_x = np.sign(x), np.abs(x)
	abs_x[abs_x<=x_cut] = x_cut
	return sign_x*(np.log10(abs_x)-np.log10(x_cut))

if slicedir==1:
	x_label, y_label, rotateX, rotateY, rotateZ = r'$Y\,[M_{\odot}]$', r'$Z\,[M_{\odot}]$', -90, 0, -90
elif slicedir==2:
	x_label, y_label, rotateX, rotateY, rotateZ = r'$X\,[M_{\odot}]$', r'$Z\,[M_{\odot}]$', -90, 0, 0
else:
	x_label, y_label, rotateX, rotateY, rotateZ = r'$X\,[M_{\odot}]$', r'$Y\,[M_{\odot}]$', 0, 0, 0
bhac_fname = os.path.join(data_dir, name_tag)
da1 = read.load(slice1, file=bhac_fname, type='vtu', silent=True, rotateX=rotateX, rotateY=rotateY, rotateZ=rotateZ)
da2 = read.load(slice2, file=bhac_fname, type='vtu', silent=True, rotateX=rotateX, rotateY=rotateY, rotateZ=rotateZ)
cell_centers = da1.getCenterPoints()

for var in interested_pars:
	print("ploting "+var)
	out_fig_fname = os.path.join(out_fig_dir, "{}.png".format(var))
	aim_var = np.zeros(da1.ncells)
	for i in range(da1.ncells):
		idx_for_da2 = da2.getIcellByPoint(*cell_centers[i])
		aim_var[i] = da2.getVar(var)[idx_for_da2]-da1.getVar(var)[i] 
	if var in log_pars:
		aim_var = strange_log(aim_var, x_cut=x_cut)
	amrplot.polyplot(aim_var, da1, cmap='jet')
	plt.title(var+r" diff at $t={:6.3f}\,ms$".format(da2.time*myunit.time*1000))
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.savefig(out_fig_fname, dpi=500)