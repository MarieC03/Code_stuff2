# -*- coding: utf-8 -*-
"""
Created on Sun May  7 22:20:09 2023

@author: Marie
"""

import kuibit 
#print(kuibit.__version__)
from kuibit import grid_data as gd
from kuibit import grid_data_utils as gdu
from kuibit.grid_data import UniformGrid
from kuibit import series
from kuibit import simdir as sd
from kuibit.simdir import SimDir

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import h5py
import math
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import scipy
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os

nFontSize=15
majorLocator = MultipleLocator(0.2)
majorFormatter = FormatStrFormatter('%d')
majorLocatorx = MultipleLocator(0.5)
minorLocatorx = MultipleLocator(0.05)


mpl.rcParams['axes.linewidth']    = 1.0
mpl.rcParams['xtick.major.size']  = 6
mpl.rcParams['xtick.major.width'] = 0.8
mpl.rcParams['xtick.minor.size']  = 4.5
mpl.rcParams['xtick.minor.width'] = 0.6

mpl.rcParams['ytick.major.size']  = 6
mpl.rcParams['ytick.major.width'] = 0.8
mpl.rcParams['ytick.minor.size']  = 4.5
mpl.rcParams['ytick.minor.width'] = 0.6


    #plot settings
params = {
        'figure.figsize'    : [8, 6],
        'legend.fontsize'   : 14,
        'text.usetex'       : True,
        'axes.titlesize' : 18,
        'axes.labelsize' : 18,  
        'xtick.labelsize' : 18 ,
        'ytick.labelsize' : 18 
  }
matplotlib.rcParams.update(params)

###################################
EoS_name = "_gold_" #"_TNTYST_" 

Flag_Save_Plots = 1

present_2D = False
present_3D = True

    # ---> make a uniform grid : ----------------------
NGRID= 400

lowbound = 0#-20#-15
upbound = 20#20# 15

upbound2D = upbound #15
lowbound2D = lowbound #-15
upboundplot = 20#20#15
lowboundplot = 0#-20#-15

######################################################
Quantity='rho' 
grid_xy = 0
grid_xz = 0
grid_xyz = 1

foldername = "Plots_"+Quantity

if not os.path.exists(foldername):
    os.mkdir(foldername)  # Creates the folder

####################################
s = SimDir(".")
    #print(s)

gf = s.gridfunctions

print(gf)

gf_xy = gf.xy
gf_xz = gf.xz
gf_xyz = gf.xyz
    #print(gf_xy)
if(grid_xy==1):
    Quantity_xy = gf_xy[Quantity]
if(grid_xz==1):
    Quantity_xy = gf_xz[Quantity]
if(present_3D):
    Quantity_xyz  = gf_xyz[Quantity]

grid_2D = UniformGrid([NGRID, NGRID], x0=[lowbound,lowbound], x1=[upbound,upbound])
grid_3D = UniformGrid([NGRID, NGRID, NGRID], x0=[lowbound,lowbound, lowbound], x1=[upbound,upbound, upbound])
    #---> print available iterations : ---------------------

if(present_3D):
    print(Quantity_xyz.available_iterations)
if(present_2D):
    print(Quantity_xy.available_iterations)
print("------------------------------------------------")

imax = 10000  
it_merger = 200704 # 0 #309248 #305152 #301056
I_start = it_merger  # set some number to start iter from
I_end = 0  # set some number to end iter

#++++++++++++++++++++++++++++++++++++++++++++++++++
itera =np.empty([imax])  
tms_name_vec=[]  

j=0

if(present_3D):
    for ii in Quantity_xyz.available_iterations:
        itera[j] = ii #rho_xy.available_iterations[j]
        #print(itera[j])
        if(itera[j]==it_merger):
            jmerger = j  
            print("j-merger=",jmerger,"with iter=",itera[j])
        tmerger = Quantity_xyz.time_at_iteration(it_merger)*0.0049267398258
        t_curr = Quantity_xyz.time_at_iteration(itera[j])*0.0049267398258-tmerger
        t_curr = round(t_curr,3)
        t_curr_name = str(t_curr)
        tms_name_vec.append(t_curr_name)
        #print(t_curr_name)
        j=j+1
if(present_2D):
    for ii in Quantity_xy.available_iterations:
        itera[j] = ii #rho_xy.available_iterations[j]
        #print(itera[j])
        if(itera[j]==it_merger):
            jmerger = j  
            print("j-merger=",jmerger,"with iter=",itera[j])
        tmerger = Quantity_xy.time_at_iteration(it_merger)*0.0049267398258
        t_curr = Quantity_xy.time_at_iteration(itera[j])*0.0049267398258-tmerger
        t_curr = round(t_curr,3)
        t_curr_name = str(t_curr)
        tms_name_vec.append(t_curr_name)
        #print(t_curr_name)
        j=j+1

jmax = j-1    
print("NUMBER of iterations",jmax,"with maximal iteration:",itera[jmax])
if(itera[jmax]<it_merger):
    print("------------------------------------------------")
    print("WARNING: available iterations don't reach merger" )
    print("------------------------------------------------")
 
if(I_start <= 0):
    i=0
elif(I_start == it_merger):
    i=jmerger
    if(it_merger < itera[0]):
        i=0
else:
    i= I_start

if(I_end ==0):
    imax = jmax
else:
    imax = I_end

print("---------------------------------")
#####################################################################

#for i in itera: 
while i<=imax:
    #print(itera)
    it = itera[i]
    tms_name = tms_name_vec[i] 

    if(present_3D):
        #--- > time at a given iteration : -----------------------
        print("At iteration",it," Time in Msun:",Quantity_xyz.time_at_iteration(it))
        #print(rho_xy.time_at_iteration(it))
    
        #-----> CORRECT : TIME IN ms
        #tmerger = 8
        #t_curr = rho_xy.time_at_iteration(it)*0.0049267398258-tmerger
        tmerger = Quantity_xyz.time_at_iteration(it_merger)*0.0049267398258
        t_curr = Quantity_xyz.time_at_iteration(it)*0.0049267398258-tmerger
        print("Post merger time in ms  :",t_curr)
        print(t_curr)
        t_curr_int = int(t_curr)
        #print(t_curr_int)
    if(present_2D):
        #--- > time at a given iteration : -----------------------
        print("At iteration",it," Time in Msun:",Quantity_xy.time_at_iteration(it))
        #print(rho_xy.time_at_iteration(it))
    
        #-----> CORRECT : TIME IN ms
        #tmerger = 8
        #t_curr = rho_xy.time_at_iteration(it)*0.0049267398258-tmerger
        tmerger = Quantity_xy.time_at_iteration(it_merger)*0.0049267398258
        t_curr = Quantity_xy.time_at_iteration(it)*0.0049267398258-tmerger
        print("Post merger time in ms  :",t_curr)
        print(t_curr)
        t_curr_int = int(t_curr)
        #print(t_curr_int)

    t_curr_int1 = 0 #20

    NDIM2 = NGRID*NGRID

    # ---> make 2D or 3D grid object rho : --------------------
    #rho_3D = rho.read_on_grid(it, grid_3D, resample=True)
    #rho_2D_xz = rho_xz.read_on_grid(it, grid_2D, resample=True)
    
    if(present_3D):
        Quantity_3D_xyz = Quantity_xyz.read_on_grid(it, grid_3D, resample=True)
        QuantityG_np = Quantity_3D_xyz.data_xyz # same as above but raw numpy 
        QuantityG_np_2d =QuantityG_np[:,:,0] 
        Quantity_2D_xy = gd.UniformGridData(grid_2D, QuantityG_np_2d)
    if(present_2D):
        Quantity_2D_xy = Quantity_xy.read_on_grid(it, grid_2D, resample=True)

    xgrid = grid_2D.coordinates(as_meshgrid=True)[0]
    ygrid = grid_2D.coordinates(as_meshgrid=True)[1]
    zgrid = grid_2D.coordinates(as_meshgrid=True)[1]

    Quantity_np = np.array(Quantity_2D_xy.data_xyz)
    Quantity_np1 = Quantity_np.reshape([NDIM2,1])
    Quantity_df1 = pd.DataFrame(Quantity_np1)
    Filename = "./Data_np/"+ Quantity +"_2D_np_saved" + EoS_name + tms_name + ".txt"
    Quantity_df1.to_csv(Filename,header = False, index= False)

    ################## 2D plot ################################
    QuantityShow  = Quantity_2D_xy.data_xyz
    
    fig = plt.figure(figsize=(6,6), dpi=128)
    ax = fig.add_subplot()
     #ax = fig.add_subplot(projection='3d', computed_zorder=True)
     
    #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
    pic1 = ax.imshow(QuantityShow, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,1
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
    plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
    if(grid_xy==1):
        gridname = "_xy_it_"
        plt.xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        plt.ylabel(r'$y \, [\rm km]$',fontsize=22, style = 'italic')
    elif(grid_xz==1):
        gridname = "_xz_it_"
        plt.xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        plt.ylabel(r'$z \, [\rm km]$',fontsize=22, style = 'italic')
    elif(grid_xyz==1):
        gridname = "_xy_it_"
        plt.xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        plt.ylabel(r'$y \, [\rm km]$',fontsize=22, style = 'italic')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('top', size="3%", pad=0.3)
    cbar = fig.colorbar(pic1,orientation='horizontal',cax=cax)   
    #cbar_ticks=np.arange(-0.088,-0.05,0.01) #-0.091,-0.051,0.01) 
    #cbar.set_ticks(cbar_ticks)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.xaxis.set_label_position('top')
    #cbar.set_label(r'$W_{\rm lorentz}$',fontsize=18, style = 'italic',fontweight='bold',color="black")
    cbarname = Quantity
    cbar.set_label(Quantity,fontsize=18, style = 'italic',fontweight='bold',color="black")
    plt.plot()
    Filename = "./"+foldername+"/" + Quantity + gridname + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    #plt.show()

    ################################ 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0
    quantity1Dx = Quantity_2D_xy.sliced([None,yslice])

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)

    quantity1 = np.array(quantity1Dx.data_xyz)
    #print(rho1.shape)

    # ------------------------------
    Quantity1x_df = pd.DataFrame(quantity1)
    Filename1 = "./Data_extract/ZSaved_" + Quantity + "_x" + EoS_name + tms_name + ".txt"
    #Rho1x_df.to_csv(Filename1,header = False, index= False)

    ax1.plot(xcoord, quantity1, color='blue')

    ax1.set_xlim(0.,upboundplot)#upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(Quantity,fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./"+foldername+"/" +Quantity+"_x_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    #plt.show()

    ################################ 1D in y #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0
    quantity1Dy = Quantity_2D_xy.sliced([yslice,None])

    #print(rho1Dx.data_xyz)
    #print(rho1Dx.shape)

    # ---> make 1D grid data to a array which can be plotted ----
    ycoord = np.linspace(lowbound,upbound,NGRID)
    print(ycoord.shape)

    quantity1 = np.array(quantity1Dy.data_xyz)
    #print(rho1.shape)

    # ------------------------------
    Quantity1y_df = pd.DataFrame(quantity1)
    Filename1 = "./Data_extract/ZSaved"+Quantity+"_y" + EoS_name + tms_name + ".txt"
    #Rho1y_df.to_csv(Filename1,header = False, index= False)


    ax1.plot(ycoord, quantity1, color='blue')

    ax1.set_xlim(0.,upboundplot)#upbound)
    if(grid_xz==1):
        gridname = "_y_it_"
        ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
    elif(grid_xz==1):
        gridname = "_z_it_"
        ax1.set_xlabel(r'$z \,[\rm km]$',fontsize=22, style = 'italic')        
    elif(grid_xyz==1):
        gridname = "_y_it_"
        ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')

    ax1.set_ylabel(Quantity,fontsize=22, style = 'italic')

    plt.plot()
    Filename = "./"+foldername+"/" +Quantity+ gridname + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    #plt.show()
    ###################################################
    #############################################################
    i = i+1
    continue