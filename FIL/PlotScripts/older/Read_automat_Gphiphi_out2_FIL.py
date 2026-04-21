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

#################################################
#-------  OPTIONS  FOR USERS : -----------------

EoS_name = "_DD2_" #"_TNTYST_" 

Flag_Save_Plots = 1
Flag_Dont_Plot = 0
Flag_no_2D_Plots = 0

Flag_press = 1
Flag_temp = 1
Flag_eps = 0
Flag_entropy = 0#1
Flag_Ye = 1

magfields = 1 # 0 for no mag.
M1_scheme = 1

Flag_M1_Enu = 1
Flag_M1_eps_nu = 1
Flag_M1_eta_nu = 1
Flag_M1_kappa_a_nu = 1
Flag_M1_kappa_s_nu = 1
Flag_M1_Q_nu = 1
Flag_M1_R_nu = 1

it_merger = 297984 #309248 #305152 #301056
I_start = it_merger  # set some number to start iter from
I_end = 0  # set some number to end iter

    # ---> make a uniform grid : ----------------------
NGRID= 400

lowbound = -25#-20#-15
upbound = 25#20# 15

upbound2D = upbound #15
lowbound2D = lowbound #-15
upboundplot = 25#20#15
lowboundplot = -25#-20#-15

#----------------------------------
######################################################


s = SimDir(".")
    #print(s)

gf = s.gridfunctions

print(gf)

gf_xy = gf.xy
gf_xz = gf.xz
gf_xyz = gf.xyz

    #print(gf_xy)

rho_xy = gf_xy['rho']
if(Flag_press==1):
    #press_xy = gf_xy['press']
    press_xy = gf_xy['P']
if(Flag_temp==1):
    #temp_xy = gf_xy['temperature']
    temp_xy = gf_xy['temp']
if(Flag_eps ==1):
    eps_xy = gf_xy['eps']
if(Flag_entropy==1):
    entropy_xy = gf_xy['entropy']

#--------------------------
# quantities needed for rotation:
alpha_xy = gf_xy['alp']
betax_xy = gf_xy['betax']
betay_xy = gf_xy['betay']
betaz_xy = gf_xy['betaz'] 

gxx_xy = gf_xy['gxx']
gxy_xy = gf_xy['gxy']
gxz_xy = gf_xy['gxz'] 
     
gyx_xy = gf_xy['gxy']
gyy_xy = gf_xy['gyy']  
gyz_xy = gf_xy['gyz'] 

gzx_xy = gf_xy['gxz'] 
gzy_xy = gf_xy['gyz']
gzz_xy = gf_xy['gzz']  

vel0_xy = gf_xy['vel[0]']
vel1_xy = gf_xy['vel[1]']
vel2_xy = gf_xy['vel[2]']
#--------------------------

if(Flag_Ye==1):
    Ye_xy = gf_xy['Y_e']

if(magfields == 1):
    Bvec0_xy = gf_xy['Bvec[0]']
    Bvec1_xy = gf_xy['Bvec[1]']
    Bvec2_xy = gf_xy['Bvec[2]']

# M1 quantities
if(M1_scheme == 1):
    if(Flag_M1_Enu==1):
        Enua_xy = gf_xy['Enue_bar']
        Enue_xy = gf_xy['Enue']
        Enux_xy = gf_xy['Enux']
    if(Flag_M1_eps_nu==1):
        eps_nua_xy = gf_xy['eps_nue_bar']
        eps_nue_xy = gf_xy['eps_nue']
        eps_nux_xy = gf_xy['eps_nux']
    if(Flag_M1_eta_nu==1):       
        eta_nua_xy = gf_xy['eta_nue_bar']
        eta_nue_xy = gf_xy['eta_nue']
        eta_nux_xy = gf_xy['eta_nux']
    if(Flag_M1_kappa_a_nu==1):       
        kappa_a_nua_xy = gf_xy['kappa_nue_bar_a']
        kappa_a_nue_xy = gf_xy['kappa_nue_a']
        kappa_a_nux_xy = gf_xy['kappa_nux_a']
    if(Flag_M1_kappa_s_nu==1):       
        kappa_s_nua_xy = gf_xy['kappa_nue_bar_s']
        kappa_s_nue_xy = gf_xy['kappa_nue_s']
        kappa_s_nux_xy = gf_xy['kappa_nux_s']
    if(Flag_M1_Q_nu==1):       
        Q_nua_xy = gf_xy['Qnue_bar']
        Q_nue_xy = gf_xy['Qnue']
        Q_nux_xy = gf_xy['Qnux']
    if(Flag_M1_R_nu==1):       
        R_nua_xy = gf_xy['Qnue_bar']
        R_nue_xy = gf_xy['Qnue']
        R_nux_xy = gf_xy['Qnux']

    # ---> prints available times : ---------------------
    #print(rho_xy.available_times)
    #---> values of rho at some timepoint time1
    #print(rho_xy.get_time(time1))

    # ---> make a uniform grid : ----------------------

grid_2D = UniformGrid([NGRID, NGRID], x0=[lowbound,lowbound], x1=[upbound,upbound])

    #grid_3D = UniformGrid([NGRID, NGRID, NGRID], x0=[-10,-10, 0], x1=[10,10,10])

    #rhoT = SimDir('.').timeseries.maximum['rho']
    #print("Rho time max")
    #print(rhoT)

    #--- > print coords of the grid with 0=x, 1=y ---------
    #print(grid_2D.coordinates()[0])

    #---> print available iterations : ---------------------
print(rho_xy.available_iterations)
print("------------------------------------------------")

imax = 10000    
#++++++++++++++++++++++++++++++++++++++++++++++++++
itera =np.empty([imax])  
tms_name_vec=[]  

j=0
for ii in rho_xy.available_iterations:
    itera[j] = ii #rho_xy.available_iterations[j]
    #print(itera[j])
    if(itera[j]==it_merger):
        jmerger = j  
        print("j-merger=",jmerger,"with iter=",itera[j])
    tmerger = rho_xy.time_at_iteration(it_merger)*0.0049267398258
    t_curr = rho_xy.time_at_iteration(itera[j])*0.0049267398258-tmerger
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
    
#exit(1) 
 
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

    #--- > time at a given iteration : -----------------------
    print("At iteration",it," Time in Msun:",rho_xy.time_at_iteration(it))
    #print(rho_xy.time_at_iteration(it))

    #-----> CORRECT : TIME IN ms
    #tmerger = 8
    #t_curr = rho_xy.time_at_iteration(it)*0.0049267398258-tmerger
    tmerger = rho_xy.time_at_iteration(it_merger)*0.0049267398258
    t_curr = rho_xy.time_at_iteration(it)*0.0049267398258-tmerger
    print("Post merger time in ms  :",t_curr)
    print(t_curr)
    t_curr_int = int(t_curr)
    #print(t_curr_int)

    t_curr_int1 = 0 #20

    # ---> make 2D or 3D grid object rho : --------------------
    #rho_3D = rho.read_on_grid(it, grid_3D, resample=True)
    #rho_2D_xz = rho_xz.read_on_grid(it, grid_2D, resample=True)

    rho_2D_xy = rho_xy.read_on_grid(it, grid_2D, resample=True)
    if(Flag_press == 1):
        press_2D_xy = press_xy.read_on_grid(it, grid_2D, resample=True)
    if(Flag_temp == 1):
        temp_2D_xy = temp_xy.read_on_grid(it, grid_2D, resample=True)
    #eps_2D_xy = eps_xy.read_on_grid(it, grid_2D, resample=True)
    if(Flag_entropy == 1):
        entropy_2D_xy = entropy_xy.read_on_grid(it, grid_2D, resample=True)

    #rho_it = rho_xy[it]
    alpha_2D_xy = alpha_xy.read_on_grid(it, grid_2D, resample=True)
    betax_2D_xy = betax_xy.read_on_grid(it, grid_2D, resample=True)
    betay_2D_xy = betay_xy.read_on_grid(it, grid_2D, resample=True)
    betaz_2D_xy = betaz_xy.read_on_grid(it, grid_2D, resample=True)

    gxx_2D_xy = gxx_xy.read_on_grid(it, grid_2D, resample=True)
    gxy_2D_xy = gxy_xy.read_on_grid(it, grid_2D, resample=True)
    gxz_2D_xy = gxz_xy.read_on_grid(it, grid_2D, resample=True)

    gyx_2D_xy = gyx_xy.read_on_grid(it, grid_2D, resample=True)
    gyy_2D_xy = gyy_xy.read_on_grid(it, grid_2D, resample=True)
    gyz_2D_xy = gyz_xy.read_on_grid(it, grid_2D, resample=True)

    gzx_2D_xy = gzx_xy.read_on_grid(it, grid_2D, resample=True)
    gzy_2D_xy = gzy_xy.read_on_grid(it, grid_2D, resample=True)
    gzz_2D_xy = gzz_xy.read_on_grid(it, grid_2D, resample=True)

    vel0_2D_xy = vel0_xy.read_on_grid(it, grid_2D, resample=True)
    vel1_2D_xy = vel1_xy.read_on_grid(it, grid_2D, resample=True)
    vel2_2D_xy = vel2_xy.read_on_grid(it, grid_2D, resample=True)
    
    Ye_2D_xy = Ye_xy.read_on_grid(it, grid_2D, resample=True)
    
    if(magfields == 1):
        Bvec0_2D_xy = Bvec0_xy.read_on_grid(it, grid_2D, resample=True)
        Bvec1_2D_xy = Bvec0_xy.read_on_grid(it, grid_2D, resample=True)
        Bvec2_2D_xy = Bvec0_xy.read_on_grid(it, grid_2D, resample=True)

    if(M1_scheme):
        if(Flag_M1_Enu==1):
            Enua_2D_xy = Enua_xy.read_on_grid(it, grid_2D, resample=True)
            Enue_2D_xy = Enue_xy.read_on_grid(it, grid_2D, resample=True)
            Enux_2D_xy = Enux_xy.read_on_grid(it, grid_2D, resample=True)
        if(Flag_M1_eps_nu==1):
            eps_nua_2D_xy = eps_nua_xy.read_on_grid(it, grid_2D, resample=True)
            eps_nue_2D_xy = eps_nue_xy.read_on_grid(it, grid_2D, resample=True)
            eps_nux_2D_xy = eps_nux_xy.read_on_grid(it, grid_2D, resample=True)
        if(Flag_M1_eta_nu==1):
            eta_nua_2D_xy = eta_nua_xy.read_on_grid(it, grid_2D, resample=True)
            eta_nue_2D_xy = eta_nue_xy.read_on_grid(it, grid_2D, resample=True)
            eta_nux_2D_xy = eta_nux_xy.read_on_grid(it, grid_2D, resample=True)
        if(Flag_M1_kappa_a_nu==1):
            kappa_a_nua_2D_xy = kappa_a_nua_xy.read_on_grid(it, grid_2D, resample=True)
            kappa_a_nue_2D_xy = kappa_a_nue_xy.read_on_grid(it, grid_2D, resample=True)
            kappa_a_nux_2D_xy = kappa_a_nux_xy.read_on_grid(it, grid_2D, resample=True)
        if(Flag_M1_kappa_s_nu==1):
            kappa_s_nua_2D_xy = kappa_a_nua_xy.read_on_grid(it, grid_2D, resample=True)
            kappa_s_nue_2D_xy = kappa_a_nue_xy.read_on_grid(it, grid_2D, resample=True)
            kappa_s_nux_2D_xy = kappa_a_nux_xy.read_on_grid(it, grid_2D, resample=True)
        if(Flag_M1_Q_nu==1):
            R_nua_2D_xy = Q_nua_xy.read_on_grid(it, grid_2D, resample=True)
            R_nue_2D_xy = Q_nue_xy.read_on_grid(it, grid_2D, resample=True)
            R_nux_2D_xy = Q_nux_xy.read_on_grid(it, grid_2D, resample=True)

    #print(vel2_2D_xy.type)
     
    #xgrid = grid_2D.coordinates()[0]
    #ygrid = grid_2D.coordinates()[1]
    #zgrid = grid_2D.coordinates()[1]

    #xgrid = grid_2D.coordinates(as_same_shape=True)[0]
    #ygrid = grid_2D.coordinates(as_same_shape=True)[1]
    #zgrid = grid_2D.coordinates(as_same_shape=True)[1]

    xgrid = grid_2D.coordinates(as_meshgrid=True)[0]
    ygrid = grid_2D.coordinates(as_meshgrid=True)[1]
    zgrid = grid_2D.coordinates(as_meshgrid=True)[1]

    # Correct Lorentz-factor : 
    Wlor_2D_xy = 1/((1-(gxx_2D_xy*vel0_2D_xy*vel0_2D_xy + gxy_2D_xy*vel0_2D_xy*vel1_2D_xy
                        +gxz_2D_xy*vel1_2D_xy*vel2_2D_xy + gyx_2D_xy*vel1_2D_xy*vel0_2D_xy
                        +gyy_2D_xy*vel1_2D_xy*vel1_2D_xy + gyz_2D_xy*vel1_2D_xy*vel2_2D_xy
                        +gzx_2D_xy*vel2_2D_xy*vel0_2D_xy + gzy_2D_xy*vel2_2D_xy*vel1_2D_xy
                        +gzz_2D_xy*vel2_2D_xy*vel0_2D_xy))**(1./2.))

    # less correct : -----------------------------------------------
    Wlor_2D_xy_new = np.multiply(1,np.reciprocal((1-(gxx_2D_xy.data_xyz*vel0_2D_xy.data_xyz*vel0_2D_xy.data_xyz 
                       + gxy_2D_xy.data_xyz*vel0_2D_xy.data_xyz*vel1_2D_xy.data_xyz
                        +gxz_2D_xy.data_xyz*vel1_2D_xy.data_xyz*vel2_2D_xy.data_xyz 
                        + gyx_2D_xy.data_xyz*vel1_2D_xy.data_xyz*vel0_2D_xy.data_xyz
                        +gyy_2D_xy.data_xyz*vel1_2D_xy.data_xyz*vel1_2D_xy.data_xyz 
                        + gyz_2D_xy.data_xyz*vel1_2D_xy.data_xyz*vel2_2D_xy.data_xyz
                        +gzx_2D_xy.data_xyz*vel2_2D_xy.data_xyz*vel0_2D_xy.data_xyz 
                        + gzy_2D_xy.data_xyz*vel2_2D_xy.data_xyz*vel1_2D_xy.data_xyz
                        +gzz_2D_xy.data_xyz*vel2_2D_xy.data_xyz*vel0_2D_xy.data_xyz))**(1./2.)) )
    Wlor_2D_xy_griddat = gd.UniformGridData(grid_2D, Wlor_2D_xy_new)
    # ----------------------------------------------------------------------
    ut_2D_xy = Wlor_2D_xy/alpha_2D_xy


    #------------------ useless part: ------------------------------------
    rows, cols = (NGRID, NGRID)
    #vphi_2D_xy = [[0]*cols]*rows
    #betaphi_2D_xy  = [[0]*cols]*rows
    #Omega_2D_xy_np = [[0]*cols]*rows
    alpha_np = np.array(alpha_2D_xy.data_xyz)
    vel0_np = np.array(vel0_2D_xy.data_xyz)
    vel1_np = np.array(vel1_2D_xy.data_xyz)
    betax_np = np.array(betax_2D_xy.data_xyz)
    betay_np = np.array(betay_2D_xy.data_xyz)
    xgrid_np = np.array(xgrid)
    ygrid_np = np.array(ygrid)
    W_np = np.array(Wlor_2D_xy.data_xyz)

    NDIM2 = NGRID*NGRID

    # nope : makes 200 x 200 file ...
    #Wlor_df = pd.DataFrame(Wlor_2D_xy.data_xyz)
    #Wlor_df.to_csv('Wlor_2D_data_xyz_saved.txt',header = False, index= False)

    #---> save 2D arrays in one column
    W_np1 = W_np.reshape([NDIM2,1])
    Wlor_df1 = pd.DataFrame(W_np1)
    Wlor_df1.to_csv('./Data_extract/Wlor_2D_np_saved.txt',header = False, index= False)
    
    Vx_np1 = vel0_np.reshape([NDIM2,1])
    Vx_df1 = pd.DataFrame(Vx_np1)
    Vx_df1.to_csv('./Data_extract/Vx_2D_np_saved.txt',header = False, index= False)

    Vy_np1 = vel1_np.reshape([NDIM2,1])
    Vy_df1 = pd.DataFrame(Vy_np1)
    Vy_df1.to_csv('./Data_extract/Vy_2D_np_saved.txt',header = False, index= False)

    Alpha_np1 = alpha_np.reshape([NDIM2,1])
    Alpha_np1_df1 = pd.DataFrame(Alpha_np1)
    Alpha_np1_df1.to_csv('./Data_extract/Alpha_2D_np_saved.txt',header = False, index= False)

    Betax_np1 = betax_np.reshape([NDIM2,1])
    Betax_np1_df1 = pd.DataFrame(Betax_np1)
    Betax_np1_df1.to_csv('./Data_extract/Betax_2D_np_saved.txt',header = False, index= False)

    Betay_np1 = betay_np.reshape([NDIM2,1])
    Betay_np1_df1 = pd.DataFrame(Betay_np1)
    Betay_np1_df1.to_csv('./Data_extract/Betay_2D_np_saved.txt',header = False, index= False)


    xgrid_df = pd.DataFrame(xgrid)
    xgrid_df.to_csv('./Data_extract/xgrid_saved.txt',header = False, index= False)

    ygrid_df = pd.DataFrame(ygrid)
    ygrid_df.to_csv('./Data_extract/ygrid_saved.txt',header = False, index= False)

    #------------------ useles upto here -------------------------------

    #vphi_2D_xy = (xgrid*vel1_2D_xy - ygrid*vel0_2D_xy)/(xgrid**2 + ygrid**2 + zgrid**2)

    #vphi_2D_xy = np.multiply((xgrid*vel1_2D_xy.data_xyz - ygrid*vel0_2D_xy.data_xyz), np.reciprocal((xgrid**2 + ygrid**2 + zgrid**2)))
    vphi_2D_xy = np.multiply((vel1_2D_xy.data_xyz*xgrid - vel0_2D_xy.data_xyz*ygrid), np.reciprocal((xgrid**2 + ygrid**2)))# + zgrid**2)))

    betaphi_2D_xy = np.multiply((xgrid*betay_2D_xy.data_xyz - ygrid*betax_2D_xy.data_xyz), np.reciprocal(xgrid**2 + ygrid**2))# + zgrid**2))
    #Omega_2D_xy = np.multiply(alpha_2D_xy, vphi_2D_xy) - betaphi_2D_xy
    #Omega_2D_xy = vphi_2D_xy*alpha_2D_xy - betaphi_2D_xy

    #print(vphi_2D_xy)

    vphi_2D_xy_griddat = gd.UniformGridData(grid_2D, vphi_2D_xy)
    betaphi_2D_xy_griddat = gd.UniformGridData(grid_2D, betaphi_2D_xy)

    Omega_2D_xy = vphi_2D_xy_griddat.data_xyz*alpha_2D_xy.data_xyz - betaphi_2D_xy_griddat.data_xyz

    Omega_2D_xy_griddat = gd.UniformGridData(grid_2D, Omega_2D_xy)

    #gphiphi_2D_xy = np.multiply((gyz_2D_xy.data_xyz*xgrid - gxz_2D_xy.data_xyz*ygrid), np.reciprocal((xgrid**2 + ygrid**2)))# + zgrid**2)))
    # works :
    #gphiphi_2D_xy = xgrid**2*gxx_2D_xy.data_xyz + ygrid**2*gyy_2D_xy.data_xyz #+ zgrid**2 
    
    # gphiphi if gyx = gxy :
    ###gphiphi_2D_xy = xgrid**2*gxx_2D_xy.data_xyz + ygrid**2*gyy_2D_xy.data_xyz -2.0*xgrid*ygrid*gxy_2D_xy.data_xyz
    gphiphi_2D_xy = xgrid**2*gxx_2D_xy.data_xyz + ygrid**2*gyy_2D_xy.data_xyz -xgrid*ygrid*gxy_2D_xy.data_xyz -xgrid*ygrid*gyx_2D_xy.data_xyz
    gphiphi_2D_xy_griddat = gd.UniformGridData(grid_2D, gphiphi_2D_xy)

    vr_2D_xy = np.multiply((vel0_2D_xy.data_xyz*xgrid + vel1_2D_xy.data_xyz*ygrid), np.reciprocal((xgrid**2 + ygrid**2)**(1./2.)))# + zgrid**2)))
    vr_2D_xy_griddat = gd.UniformGridData(grid_2D, vr_2D_xy)

    gphir_2D_xy = ( ygrid*xgrid*np.reciprocal((xgrid**2 + ygrid**2)**(1./2.))*gxx_2D_xy.data_xyz + ygrid**2*np.reciprocal((xgrid**2 + ygrid**2)**(1./2.))*gxy_2D_xy.data_xyz 
                 - xgrid**2 *np.reciprocal((xgrid**2 + ygrid**2)**(1./2.)) *gyx_2D_xy.data_xyz - xgrid*ygrid*np.reciprocal((xgrid**2 + ygrid**2)**(1./2.))*gyy_2D_xy.data_xyz)
    gphir_2D_xy_griddat = gd.UniformGridData(grid_2D, gphir_2D_xy)
    
    gphiz_2D_xy = ygrid*gxz_2D_xy.data_xyz - xgrid*gyz_2D_xy.data_xyz
    gphiz_2D_xy_griddat = gd.UniformGridData(grid_2D, gphiz_2D_xy)
    
    #Jmom_2D_xy = np.multiply((Wlor_2D_xy.data_xyz*Wlor_2D_xy.data_xyz), np.reciprocal(alpha_2D_xy.data_xyz)) *vphi_2D_xy_griddat.data_xyz * gphiphi_2D_xy_griddat.data_xyz
    #Jmom_2D_xy = np.multiply((Wlor_2D_xy.data_xyz*Wlor_2D_xy.data_xyz), np.reciprocal(alpha_2D_xy.data_xyz)) *vphi_2D_xy_griddat* gphiphi_2D_xy_griddat
    #Jmom_2D_xy = Wlor_2D_xy*Wlor_2D_xy/ alpha_2D_xy *vphi_2D_xy_griddat* gphiphi_2D_xy_griddat
    
    # Jmom only g_phiphi
    ####Jmom_2D_xy = np.multiply((Wlor_2D_xy.data_xyz*Wlor_2D_xy.data_xyz), np.reciprocal(alpha_2D_xy.data_xyz)) *vphi_2D_xy_griddat.data_xyz * gphiphi_2D_xy_griddat.data_xyz
    
    # Jmom correct
    Jmom_2D_xy = np.multiply((Wlor_2D_xy.data_xyz*Wlor_2D_xy.data_xyz), np.reciprocal(alpha_2D_xy.data_xyz)) *(vphi_2D_xy_griddat.data_xyz * gphiphi_2D_xy_griddat.data_xyz 
                    + vr_2D_xy_griddat.data_xyz*gphir_2D_xy_griddat.data_xyz + vel2_2D_xy.data_xyz*gphiz_2D_xy_griddat.data_xyz)
    
    Jmom_2D_xy_griddat = gd.UniformGridData(grid_2D, Jmom_2D_xy)
    
    #---> save 2D arrays in one column
    Omega_np = np.array(Omega_2D_xy_griddat.data_xyz)
    Om_np1 = Omega_np.reshape([NDIM2,1])
    Om_df1 = pd.DataFrame(Om_np1)
    Filename = "./Data_np/Omega_2D_np_saved" + EoS_name + tms_name + ".txt"
    Om_df1.to_csv(Filename,header = False, index= False)

    Jmom_np = np.array(Jmom_2D_xy_griddat.data_xyz)
    Jmom_np1 = Jmom_np.reshape([NDIM2,1])
    Jmom_df1 = pd.DataFrame(Jmom_np1)
    Filename = "./Data_np/Jmom_2D_np_saved" + EoS_name + tms_name + ".txt"
    Jmom_df1.to_csv(Filename,header = False, index= False)

    Rho_np = np.array(rho_2D_xy.data_xyz)
    Rho_np1 = Rho_np.reshape([NDIM2,1])
    Rho_df1 = pd.DataFrame(Rho_np1)
    Filename = "./Data_np/Rho_2D_np_saved" + EoS_name + tms_name + ".txt"
    Rho_df1.to_csv(Filename,header = False, index= False)

    if(Flag_temp==1):
        Temp_np = np.array(temp_2D_xy.data_xyz)
        Temp_np1 = Temp_np.reshape([NDIM2,1])
        Temp_df1 = pd.DataFrame(Temp_np1)
        Filename = "./Data_np/Temp_2D_np_saved" + EoS_name + tms_name + ".txt"
        Temp_df1.to_csv(Filename,header = False, index= False)

    if(Flag_press ==1):
        Press_np = np.array(press_2D_xy.data_xyz)
        Press_np1 = Press_np.reshape([NDIM2,1])
        Press_df1 = pd.DataFrame(Press_np1)
        Filename = "./Data_np/Press_2D_np_saved" + EoS_name + tms_name + ".txt"
        Press_df1.to_csv(Filename,header = False, index= False)
    
    if(Flag_entropy == 1):
        Entropy_np = np.array(entropy_2D_xy.data_xyz)
        Entropy_np1 = Entropy_np.reshape([NDIM2,1])
        Entropy_df1 = pd.DataFrame(Entropy_np1)
        Filename = "./Data_np/Entropy_2D_np_saved" + EoS_name + tms_name + ".txt"
        Entropy_df1.to_csv(Filename,header = False, index= False)

    if(Flag_Ye == 1):
        Ye_np = np.array(Ye_2D_xy.data_xyz)
        Ye_np1 = Ye_np.reshape([NDIM2,1])
        Ye_df1 = pd.DataFrame(Ye_np1)
        Filename = "./Data_np/Ye_2D_np_saved" + EoS_name + tms_name + ".txt"
        Ye_df1.to_csv(Filename,header = False, index= False)
    
      
        
    Gphiphi_np = np.array(gphiphi_2D_xy_griddat.data_xyz)
    Gphiphi_np1 = Gphiphi_np.reshape([NDIM2,1])
    Gphiphi_df1 = pd.DataFrame(Gphiphi_np1)
    Filename = "./Data_np/Gphi2_2D_np_saved" + EoS_name + tms_name + ".txt"
    Gphiphi_df1.to_csv(Filename,header = False, index= False)

    Gphir_np = np.array(gphir_2D_xy_griddat.data_xyz)
    Gphir_np1 = Gphir_np.reshape([NDIM2,1])
    Gphir_df1 = pd.DataFrame(Gphir_np1)
    Filename = "./Data_np/Gphir_2D_np_saved" + EoS_name + tms_name + ".txt"
    Gphir_df1.to_csv(Filename,header = False, index= False)
    
    Gphiz_np = np.array(gphiz_2D_xy_griddat.data_xyz)
    Gphiz_np1 = Gphiz_np.reshape([NDIM2,1])
    Gphiz_df1 = pd.DataFrame(Gphiz_np1)
    Filename = "./Data_np/Gphiz_2D_np_saved" + EoS_name + tms_name + ".txt"
    Gphiz_df1.to_csv(Filename,header = False, index= False)

    Vphi_np = np.array(vphi_2D_xy_griddat.data_xyz)
    Vphi_np1 = Vphi_np.reshape([NDIM2,1])
    Vphi_df1 = pd.DataFrame(Vphi_np1)
    Filename = "./Data_np/Vphi_2D_np_saved" + EoS_name + tms_name + ".txt"
    Vphi_df1.to_csv(Filename,header = False, index= False)

    Vr_np = np.array(vr_2D_xy_griddat.data_xyz)
    Vr_np1 = Vr_np.reshape([NDIM2,1])
    Vr_df1 = pd.DataFrame(Vr_np1)
    Filename = "./Data_np/Vr_2D_np_saved" + EoS_name + tms_name + ".txt"
    Vr_df1.to_csv(Filename,header = False, index= False)
   
    Alpha_np = np.array(alpha_2D_xy.data_xyz)
    Alpha_np1 = Alpha_np.reshape([NDIM2,1])
    Alpha_df1 = pd.DataFrame(Alpha_np1)
    Filename = "./Data_np/Alpha_2D_np_saved" + EoS_name + tms_name + ".txt"
    Alpha_df1.to_csv(Filename,header = False, index= False)
   
    Betaphi_np = np.array(betaphi_2D_xy_griddat.data_xyz)
    Betaphi_np1 = Betaphi_np.reshape([NDIM2,1])
    Betaphi_df1 = pd.DataFrame(Betaphi_np1)
    Filename = "./Data_np/Betaphi_2D_np_saved" + EoS_name + tms_name + ".txt"
    Betaphi_df1.to_csv(Filename,header = False, index= False)
    
    if(magfields == 1):
        Bvec0_np = np.array(Bvec0_2D_xy.data_xyz)
        Bvec0_np1 = Bvec0_np.reshape([NDIM2,1])
        Bvec0_df1 = pd.DataFrame(Bvec0_np1)
        Filename = "./Data_np/Bvec0_2D_np_saved" + EoS_name + tms_name + ".txt"
        Bvec0_df1.to_csv(Filename,header = False, index= False)
 
        Bvec1_np = np.array(Bvec1_2D_xy.data_xyz)
        Bvec1_np1 = Bvec1_np.reshape([NDIM2,1])
        Bvec1_df1 = pd.DataFrame(Bvec1_np1)
        Filename = "./Data_np/Bvec1_2D_np_saved" + EoS_name + tms_name + ".txt"
        Bvec1_df1.to_csv(Filename,header = False, index= False)
   
        Bvec2_np = np.array(Bvec2_2D_xy.data_xyz)
        Bvec2_np1 = Bvec2_np.reshape([NDIM2,1])
        Bvec2_df1 = pd.DataFrame(Bvec2_np1)
        Filename = "./Data_np/Bvec2_2D_np_saved" + EoS_name + tms_name + ".txt"
        Bvec2_df1.to_csv(Filename,header = False, index= False)
   
    # M1:
    if(M1_scheme == 1):   
        if(Flag_M1_Enu ==1): 
            Enua_np = np.array(Enua_2D_xy.data_xyz)
            Enua_np1 = Enua_np.reshape([NDIM2,1])
            Enua_df1 = pd.DataFrame(Enua_np1)
            Filename = "./Data_np/Enua_2D_np_saved" + EoS_name + tms_name + ".txt"
            Enua_df1.to_csv(Filename,header = False, index= False)
     
            Enue_np = np.array(Enue_2D_xy.data_xyz)
            Enue_np1 = Enue_np.reshape([NDIM2,1])
            Enue_df1 = pd.DataFrame(Enue_np1)
            Filename = "./Data_np/Enue_2D_np_saved" + EoS_name + tms_name + ".txt"
            Enue_df1.to_csv(Filename,header = False, index= False)
    
            Enux_np = np.array(Enux_2D_xy.data_xyz)
            Enux_np1 = Enux_np.reshape([NDIM2,1])
            Enux_df1 = pd.DataFrame(Enux_np1)
            Filename = "./Data_np/Enux_2D_np_saved" + EoS_name + tms_name + ".txt"
            Enux_df1.to_csv(Filename,header = False, index= False)

        if(Flag_M1_eps_nu ==1):
            eps_nua_np = np.array(eps_nua_2D_xy.data_xyz)
            eps_nua_np1 = eps_nua_np.reshape([NDIM2,1])
            eps_nua_df1 = pd.DataFrame(eps_nua_np1)
            Filename = "./Data_np/Eps_nua_2D_np_saved" + EoS_name + tms_name + ".txt"
            eps_nua_df1.to_csv(Filename,header = False, index= False)
    
            eps_nue_np = np.array(eps_nue_2D_xy.data_xyz)
            eps_nue_np1 = eps_nue_np.reshape([NDIM2,1])
            eps_nue_df1 = pd.DataFrame(eps_nue_np1)
            Filename = "./Data_np/Eps_nue_2D_np_saved" + EoS_name + tms_name + ".txt"
            eps_nue_df1.to_csv(Filename,header = False, index= False)
                                              
            eps_nux_np = np.array(eps_nux_2D_xy.data_xyz)
            eps_nux_np1 = eps_nux_np.reshape([NDIM2,1])
            eps_nux_df1 = pd.DataFrame(eps_nux_np1)
            Filename = "./Data_np/Eps_nux_2D_np_saved" + EoS_name + tms_name + ".txt"
            eps_nux_df1.to_csv(Filename,header = False, index= False)

        if(Flag_M1_eta_nu ==1):
            eta_nua_np = np.array(eta_nua_2D_xy.data_xyz)
            eta_nua_np1 = eta_nua_np.reshape([NDIM2,1])
            eta_nua_df1 = pd.DataFrame(eta_nua_np1)
            Filename = "./Data_np/Eta_nua_2D_np_saved" + EoS_name + tms_name + ".txt"
            eta_nua_df1.to_csv(Filename,header = False, index= False)
    
            eta_nue_np = np.array(eta_nue_2D_xy.data_xyz)
            eta_nue_np1 = eta_nue_np.reshape([NDIM2,1])
            eta_nue_df1 = pd.DataFrame(eta_nue_np1)
            Filename = "./Data_np/Eta_nue_2D_np_saved" + EoS_name + tms_name + ".txt"
            eta_nue_df1.to_csv(Filename,header = False, index= False)
                                           
            eta_nux_np = np.array(eta_nux_2D_xy.data_xyz)
            eta_nux_np1 = eta_nux_np.reshape([NDIM2,1])
            eta_nux_df1 = pd.DataFrame(eta_nux_np1)
            Filename = "./Data_np/Eta_nux_2D_np_saved" + EoS_name + tms_name + ".txt"
            eta_nux_df1.to_csv(Filename,header = False, index= False)

        if(Flag_M1_kappa_a_nu ==1):
            kappa_a_nua_np = np.array(kappa_a_nua_2D_xy.data_xyz)
            kappa_a_nua_np1 = kappa_a_nua_np.reshape([NDIM2,1])
            kappa_a_nua_df1 = pd.DataFrame(kappa_a_nua_np1)
            Filename = "./Data_np/Kappa_a_nua_2D_np_saved" + EoS_name + tms_name + ".txt"
            kappa_a_nua_df1.to_csv(Filename,header = False, index= False)
    
            kappa_a_nue_np = np.array(kappa_a_nue_2D_xy.data_xyz)
            kappa_a_nue_np1 = eta_nue_np.reshape([NDIM2,1])
            kappa_a_nue_df1 = pd.DataFrame(kappa_a_nue_np1)
            Filename = "./Data_np/Kappa_a_nue_2D_np_saved" + EoS_name + tms_name + ".txt"
            kappa_a_nue_df1.to_csv(Filename,header = False, index= False)
                                           
            kappa_a_nux_np = np.array(kappa_a_nux_2D_xy.data_xyz)
            kappa_a_nux_np1 = kappa_a_nux_np.reshape([NDIM2,1])
            kappa_a_nux_df1 = pd.DataFrame(eta_nux_np1)
            Filename = "./Data_np/Kappa_a_nux_2D_np_saved" + EoS_name + tms_name + ".txt"
            kappa_a_nux_df1.to_csv(Filename,header = False, index= False)

    if(Flag_Dont_Plot==1):
       exit(1)   

    #####################  TEEEST ##########################################################
    ######################################## RHO  2D ###################################
    fig = plt.figure(figsize=(6,6), dpi=128)
    ax = fig.add_subplot()
     #ax = fig.add_subplot(projection='3d', computed_zorder=True)
     
    #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
    pic1 = ax.imshow(rho_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
    plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
    cbar.set_label(r'$\rho $',fontsize=18, style = 'italic',fontweight='bold',color="black")

    plt.plot()
    Filename = "./Plots/Rho_xy_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()


              
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #############################################################################################

    ################################ RHO 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0
    rho1Dx = rho_2D_xy.sliced([None,yslice])

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)

    rho1 = np.array(rho1Dx.data_xyz)
    #print(rho1.shape)

    # ------------------------------
    Rho1x_df = pd.DataFrame(rho1)
    Filename1 = "./Data_extract/ZSaved_Rho_x" + EoS_name + tms_name + ".txt"
    #Rho1x_df.to_csv(Filename1,header = False, index= False)

    ax1.plot(xcoord, rho1, color='blue')

    ax1.set_xlim(0.,upboundplot)#upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$\rho $',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Rho_x_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ RHO 1D in y #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0
    rho1Dy = rho_2D_xy.sliced([yslice,None])

    #print(rho1Dx.data_xyz)
    #print(rho1Dx.shape)

    # ---> make 1D grid data to a array which can be plotted ----
    ycoord = np.linspace(lowbound,upbound,NGRID)
    print(ycoord.shape)

    rho1 = np.array(rho1Dy.data_xyz)
    #print(rho1.shape)

    # ------------------------------
    Rho1y_df = pd.DataFrame(rho1)
    Filename1 = "./Data_extract/ZSaved_Rho_y" + EoS_name + tms_name + ".txt"
    #Rho1y_df.to_csv(Filename1,header = False, index= False)


    ax1.plot(ycoord, rho1, color='blue')

    ax1.set_xlim(0.,upboundplot)#upbound)

    ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$\rho $',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Rho_y_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

            ################################ Entropy 1D in x #################################
    if(Flag_entropy==1):
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax1 = fig.add_subplot()

            yslice = 0.0
            entropy1Dx = entropy_2D_xy.sliced([None,yslice])


            # ---> make 1D grid data to a array which can be plotted ----
            xcoord = np.linspace(lowbound,upbound,NGRID)
            print(xcoord.shape)

            entropy1 = np.array(entropy1Dx.data_xyz)

            # ------------------------------
            Entropy1x_df = pd.DataFrame(entropy1)
            Filename1 = "./Data_extract/ZSaved_Entropy_x" + EoS_name + tms_name + ".txt"
            #Entropy1x_df.to_csv(Filename1,header = False, index= False)

            ax1.plot(xcoord, entropy1, color='blue')

            ax1.set_xlim(0.,upboundplot)#upbound)

            ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
            ax1.set_ylabel(r'$H $',fontsize=22, style = 'italic')


            plt.plot()
            Filename = "./Plots/Entropy_x_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()
            
        ################################ TEMP 1D in x #################################
    if(Flag_temp ==1):
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        yslice = 0.0
        temp1Dx = temp_2D_xy.sliced([None,yslice])


        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        temp1 = np.array(temp1Dx.data_xyz)

        # ------------------------------
        Temp1x_df = pd.DataFrame(temp1)
        Filename1 = "./Data_extract/ZSaved_Temp_x" + EoS_name + tms_name + ".txt"
        #Temp1x_df.to_csv(Filename1,header = False, index= False)

        ax1.plot(xcoord, temp1, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$T \, [MeV] $',fontsize=22, style = 'italic')


        plt.plot()
        Filename = "./Plots/Temp_x_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()

        ################################ TEMP 1D in y #################################
    if(Flag_temp ==1):
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        xslice = 0.0
        temp1Dy = temp_2D_xy.sliced([xslice,None])


        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        temp1 = np.array(temp1Dy.data_xyz)

        # ------------------------------
        Temp1x_df = pd.DataFrame(temp1)
        Filename1 = "./Data_extract/ZSaved_Temp_y" + EoS_name + tms_name + ".txt"
        #Temp1x_df.to_csv(Filename1,header = False, index= False)

        ax1.plot(xcoord, temp1, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$T \, [MeV] $',fontsize=22, style = 'italic')


        plt.plot()
        Filename = "./Plots/Temp_y_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()


        ################################ Press 1D in x #################################
    if(Flag_press==1):
        
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        yslice = 0.0
        press1Dx = press_2D_xy.sliced([None,yslice])

        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        press1 = np.array(press1Dx.data_xyz)

        # ------------------------------
        Press1x_df = pd.DataFrame(press1)
        Filename1 = "./Data_extract/ZSaved_Press_x" + EoS_name + tms_name + ".txt"

        ax1.plot(xcoord, press1, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$ P $',fontsize=22, style = 'italic')


        plt.plot()
        Filename = "./Plots/Press_x_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()

        ################################ Press 1D in y #################################
    if(Flag_press==1):
        
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        xslice = 0.0
        press1Dy = press_2D_xy.sliced([xslice,None])

        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        press1 = np.array(press1Dy.data_xyz)

        # ------------------------------
        Press1x_df = pd.DataFrame(press1)
        Filename1 = "./Data_extract/ZSaved_Press_y" + EoS_name + tms_name + ".txt"

        ax1.plot(xcoord, press1, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$ P $',fontsize=22, style = 'italic')


        plt.plot()
        Filename = "./Plots/Press_y_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()

        ################################ Ye 1D in x #################################
    if(Flag_Ye==1):
        
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        yslice = 0.0
        Ye1Dx = Ye_2D_xy.sliced([None,yslice])

        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        Ye1 = np.array(Ye1Dx.data_xyz)

        # ------------------------------
        Ye1x_df = pd.DataFrame(Ye1)
        Filename1 = "./Data_extract/ZSaved_Ye_x" + EoS_name + tms_name + ".txt"

        ax1.plot(xcoord, Ye1, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$ Y_e $',fontsize=22, style = 'italic')

        plt.plot()
        Filename = "./Plots/Ye_x_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()

        ################################ Ye 1D in y #################################
    if(Flag_Ye==1):
        
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        xslice = 0.0
        Ye1Dy = Ye_2D_xy.sliced([xslice,None])

        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        Ye1 = np.array(Ye1Dy.data_xyz)

        # ------------------------------
        Ye1x_df = pd.DataFrame(Ye1)
        Filename1 = "./Data_extract/ZSaved_Ye_y" + EoS_name + tms_name + ".txt"

        ax1.plot(xcoord, Ye1, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$ Y_e $',fontsize=22, style = 'italic')

        plt.plot()
        Filename = "./Plots/Ye_y_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()


        ################################ Enue 1D in x #################################
    if(M1_scheme ==1):
        if(Flag_M1_Enu ==1):
        
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax1 = fig.add_subplot()
    
            yslice = 0.0
            Enue1Dx = Enue_2D_xy.sliced([None,yslice])
    
            # ---> make 1D grid data to a array which can be plotted ----
            xcoord = np.linspace(lowbound,upbound,NGRID)
            print(xcoord.shape)
    
            Enue1 = np.array(Enue1Dx.data_xyz)
    
            # ------------------------------
            Enue1x_df = pd.DataFrame(Enue1)
            Filename1 = "./Data_extract/ZSaved_Enue_x" + EoS_name + tms_name + ".txt"
            #Press1x_df.to_csv(Filename1,header = False, index= False)
    
            ax1.plot(xcoord, Enue1, color='blue')
    
            ax1.set_xlim(0.,upboundplot)#upbound)
    
            ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
            ax1.set_ylabel(r'$ E_{\nu,e} $',fontsize=22, style = 'italic')
    
            plt.plot()
            Filename = "./Plots/Enue_x_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()

        ################################ Enue 1D in y #################################
    if(M1_scheme ==1):
        if(Flag_M1_Enu ==1):
        
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax1 = fig.add_subplot()
    
            xslice = 0.0
            Enue1Dy = Enue_2D_xy.sliced([xslice,None])
    
            # ---> make 1D grid data to a array which can be plotted ----
            xcoord = np.linspace(lowbound,upbound,NGRID)
            print(xcoord.shape)
    
            Enue1 = np.array(Enue1Dy.data_xyz)
    
            # ------------------------------
            Enue1x_df = pd.DataFrame(Enue1)
            Filename1 = "./Data_extract/ZSaved_Enue_y" + EoS_name + tms_name + ".txt"
            #Press1x_df.to_csv(Filename1,header = False, index= False)
    
            ax1.plot(xcoord, Enue1, color='blue')
    
            ax1.set_xlim(0.,upboundplot)#upbound)
    
            ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
            ax1.set_ylabel(r'$ E_{\nu,e} $',fontsize=22, style = 'italic')
    
            plt.plot()
            Filename = "./Plots/Enue_y_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()
        
        ################################ Enua 1D in x #################################
    if(M1_scheme ==1):        
        if(Flag_M1_Enu ==1):
        
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax1 = fig.add_subplot()
    
            yslice = 0.0
            Enua1Dx = Enua_2D_xy.sliced([None,yslice])
    
            # ---> make 1D grid data to a array which can be plotted ----
            xcoord = np.linspace(lowbound,upbound,NGRID)
            print(xcoord.shape)
    
            Enua1 = np.array(Enua1Dx.data_xyz)
    
            # ------------------------------
            Enue1x_df = pd.DataFrame(Enue1)
            Filename1 = "./Data_extract/ZSaved_Enua_x" + EoS_name + tms_name + ".txt"
            #Press1x_df.to_csv(Filename1,header = False, index= False)
    
            ax1.plot(xcoord, Enua1, color='blue')
    
            ax1.set_xlim(0.,upboundplot)#upbound)
    
            ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
            ax1.set_ylabel(r'$ E_{\nu,a} $',fontsize=22, style = 'italic')
    
            plt.plot()
            Filename = "./Plots/Enua_x_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()


        ################################ eps_nue 1D in x #################################
    if(M1_scheme ==1):         
        if(Flag_M1_eps_nu ==1):
        
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax1 = fig.add_subplot()
    
            yslice = 0.0
            Eps_nue1Dx = eps_nue_2D_xy.sliced([None,yslice])
    
            # ---> make 1D grid data to a array which can be plotted ----
            xcoord = np.linspace(lowbound,upbound,NGRID)
            print(xcoord.shape)
    
            Eps_nue1 = np.array(Eps_nue1Dx.data_xyz)
    
            # ------------------------------
            Eps_nue1x_df = pd.DataFrame(Eps_nue1)
            Filename1 = "./Data_extract/ZSaved_Eps_nue_x" + EoS_name + tms_name + ".txt"
            #Press1x_df.to_csv(Filename1,header = False, index= False)
    
            ax1.plot(xcoord, Eps_nue1, color='blue')
    
            ax1.set_xlim(0.,upboundplot)#upbound)
    
            ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
            ax1.set_ylabel(r'$ \epsilon_{\nu,e} $',fontsize=22, style = 'italic')
    
            plt.plot()
            Filename = "./Plots/Eps_nue_x_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()

        ################################ eps_nua 1D in x #################################
    if(M1_scheme ==1):
        if(Flag_M1_eps_nu ==1):
        
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax1 = fig.add_subplot()
    
            yslice = 0.0
            Eps_nua1Dx = eps_nua_2D_xy.sliced([None,yslice])
    
            # ---> make 1D grid data to a array which can be plotted ----
            xcoord = np.linspace(lowbound,upbound,NGRID)
            print(xcoord.shape)
    
            Eps_nua1 = np.array(Eps_nua1Dx.data_xyz)
    
            # ------------------------------
            Eps_nua1x_df = pd.DataFrame(Eps_nua1)
            Filename1 = "./Data_extract/ZSaved_Eps_nua_x" + EoS_name + tms_name + ".txt"
            #Press1x_df.to_csv(Filename1,header = False, index= False)
    
            ax1.plot(xcoord, Eps_nua1, color='blue')
    
            ax1.set_xlim(0.,upboundplot)#upbound)
    
            ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
            ax1.set_ylabel(r'$ \epsilon_{\nu,a} $',fontsize=22, style = 'italic')
    
            plt.plot()
            Filename = "./Plots/Eps_nua_x_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()
        ################################ eps_nux 1D in x #################################
    if(M1_scheme ==1):
        if(Flag_M1_eps_nu ==1):
        
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax1 = fig.add_subplot()
    
            yslice = 0.0
            Eps_nux1Dx = eps_nux_2D_xy.sliced([None,yslice])
    
            # ---> make 1D grid data to a array which can be plotted ----
            xcoord = np.linspace(lowbound,upbound,NGRID)
            print(xcoord.shape)
    
            Eps_nux1 = np.array(Eps_nux1Dx.data_xyz)
    
            # ------------------------------
            Eps_nux1x_df = pd.DataFrame(Eps_nux1)
            Filename1 = "./Data_extract/ZSaved_Eps_nux_x" + EoS_name + tms_name + ".txt"
            #Press1x_df.to_csv(Filename1,header = False, index= False)
    
            ax1.plot(xcoord, Eps_nux1, color='blue')
    
            ax1.set_xlim(0.,upboundplot)#upbound)
    
            ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
            ax1.set_ylabel(r'$ \epsilon_{\nu,x} $',fontsize=22, style = 'italic')
    
            plt.plot()
            Filename = "./Plots/Eps_nux_x_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()

        ################################ Bvec0 1D in x #################################
    if(magfields ==1):
        
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        yslice = 0.0
        Bvec01Dx = Bvec0_2D_xy.sliced([None,yslice])

        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        Bvec01 = np.array(Bvec01Dx.data_xyz)

        # ------------------------------
        Bvec01x_df = pd.DataFrame(Bvec01)
        Filename1 = "./Data_extract/ZSaved_Bvec0_x" + EoS_name + tms_name + ".txt"
        #Press1x_df.to_csv(Filename1,header = False, index= False)

        ax1.plot(xcoord, Bvec01, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$ B_{0} $',fontsize=22, style = 'italic')

        plt.plot()
        Filename = "./Plots/Bvec0_x_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()

        ################################ Bvec1 1D in x #################################
    if(magfields ==1):
        
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        yslice = 0.0
        Bvec11Dx = Bvec1_2D_xy.sliced([None,yslice])

        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        Bvec11 = np.array(Bvec11Dx.data_xyz)

        # ------------------------------
        Bvec11x_df = pd.DataFrame(Bvec11)
        Filename1 = "./Data_extract/ZSaved_Bvec1_x" + EoS_name + tms_name + ".txt"
        #Press1x_df.to_csv(Filename1,header = False, index= False)

        ax1.plot(xcoord, Bvec11, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$ B_{1} $',fontsize=22, style = 'italic')

        plt.plot()
        Filename = "./Plots/Bvec1_x_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()

        ################################ Bvec2 1D in x #################################
    if(magfields ==1):
        
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax1 = fig.add_subplot()

        yslice = 0.0
        Bvec21Dx = Bvec2_2D_xy.sliced([None,yslice])

        # ---> make 1D grid data to a array which can be plotted ----
        xcoord = np.linspace(lowbound,upbound,NGRID)
        print(xcoord.shape)

        Bvec21 = np.array(Bvec21Dx.data_xyz)

        # ------------------------------
        Bvec21x_df = pd.DataFrame(Bvec21)
        Filename1 = "./Data_extract/ZSaved_Bvec2_x" + EoS_name + tms_name + ".txt"
        #Press1x_df.to_csv(Filename1,header = False, index= False)

        ax1.plot(xcoord, Bvec21, color='blue')

        ax1.set_xlim(0.,upboundplot)#upbound)

        ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
        ax1.set_ylabel(r'$ B_{2} $',fontsize=22, style = 'italic')

        plt.plot()
        Filename = "./Plots/Bvec2_x_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()
      
                       
    ################################ W lorentz 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0
    Wlor1Dx = Wlor_2D_xy.sliced([None,yslice]) 

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)

    Wlor1 = np.array(Wlor1Dx.data_xyz)

    ax1.plot(xcoord, Wlor1, color='blue')

    ax1.set_xlim(0.,upboundplot)#upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$W$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Wlor_x_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ u^t 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0
    ut1Dx = ut_2D_xy.sliced([None,yslice]) 

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)

    ut1 = np.array(ut1Dx.data_xyz)

    ax1.plot(xcoord, ut1, color='blue')

    ax1.set_xlim(0.,upboundplot)#upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$u^t$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Ut_x_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ Vphi 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0

    vphi_1Dx = vphi_2D_xy_griddat.sliced([None,yslice])

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)

    Data_type = object

    #Omega1 = np.array(Omega_1Dx.data, dtype=Data_type)#.data_xyz)
    vphi1 = np.array(vphi_1Dx.data_xyz)

    # ------------------------------
    Vphi1x_df = pd.DataFrame(vphi1)
    Filename1 = "./Data_extract/ZSaved_Vphi_x" + EoS_name + tms_name + ".txt"
    #Vphi1x_df.to_csv(Filename1,header = False, index= False)


    ax1.plot(xcoord, vphi1, color='blue')

    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$v_{\phi}$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Vphi_x_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ Vphi 1D in y #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    xslice = 0.0

    vphi_1Dy = vphi_2D_xy_griddat.sliced([xslice,None])

    print("here")
    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)

    Data_type = object

    #Omega1 = np.array(Omega_1Dx.data, dtype=Data_type)#.data_xyz)
    vphi1 = np.array(vphi_1Dy.data_xyz)

    ax1.plot(xcoord, vphi1, color='blue')
    #ax1.plot(xcoord, Omega1, color='blue')

    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$v_{\phi}$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Vphi_y_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()


    ################################ Vx 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0

    vx_1Dx = vel0_2D_xy.sliced([None,yslice])
    vz_1Dx = vel2_2D_xy.sliced([None,yslice])


    print("here")
    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(-15,15,NGRID)
    print(xcoord.shape)


    vx1 = np.array(vx_1Dx.data_xyz)
    vz1 = np.array(vz_1Dx.data_xyz)


    ax1.plot(xcoord, vx1, color='blue')
    #ax1.plot(xcoord, Omega1, color='blue')

    ax1.set_xlim(0.,upboundplot)#upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$v_x$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Vx_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ Vy 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0

    vx_1Dx = vel0_2D_xy.sliced([None,yslice])
    vy_1Dx = vel1_2D_xy.sliced([None,yslice])
    vz_1Dx = vel2_2D_xy.sliced([None,yslice])


    print("here")
    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)


    vx1 = np.array(vx_1Dx.data_xyz)
    vy1 = np.array(vy_1Dx.data_xyz)
    vz1 = np.array(vz_1Dx.data_xyz)


    ax1.plot(xcoord, vx1, color='blue')
    #ax1.plot(xcoord, Omega1, color='blue')

    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$v_y$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Vy_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()


    ################################ Vz 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0

    vx_1Dx = vel0_2D_xy.sliced([None,yslice])
    vz_1Dx = vel2_2D_xy.sliced([None,yslice])

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)


    vx1 = np.array(vx_1Dx.data_xyz)
    vz1 = np.array(vz_1Dx.data_xyz)

    ax1.plot(xcoord, vz1, color='blue')
    #ax1.plot(xcoord, Omega1, color='blue')

    ax1.set_xlim(0.,upboundplot)#upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$v_z$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Vz_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()


    ################################ Gphiphi 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0
    #just slicing 2D Omega
    Gphiphi1Dx = gphiphi_2D_xy_griddat.sliced([None,yslice])

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)

    Gphiphi1 = np.array(Gphiphi1Dx.data_xyz)
    ax1.plot(xcoord, Gphiphi1, color='blue')

    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$g_{\phi \phi}$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Gphiphi_x_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ Gphiphi 1D in y #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    xslice = 0.0
    #just slicing 2D Omega
    Gphiphi1Dx = gphiphi_2D_xy_griddat.sliced([xslice,None])

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)

    Gphiphi1 = np.array(Gphiphi1Dx.data_xyz)
    ax1.plot(xcoord, Gphiphi1, color='blue')

    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$g_{\phi \phi}$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Gphiphi_y_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()


    ################################ OMEGA 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0

    #works and gives same result
    alpha_1Dx = alpha_2D_xy.sliced([None,yslice])
    vphi_1Dx = vphi_2D_xy_griddat.sliced([None,yslice])
    betaphi_1Dx = betaphi_2D_xy_griddat.sliced([None,yslice])
    #Omega1Dx = alpha_1Dx*vphi_1Dx - betaphi_1Dx

    #just slicing 2D Omega
    Omega1Dx = Omega_2D_xy_griddat.sliced([None,yslice])

    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)

    Data_type = object

    Omega1 = np.array(Omega1Dx.data_xyz)
   
    ax1.plot(xcoord, Omega1, color='blue')

    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$\Omega$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Omega_x_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
         plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ OMEGA 1D in -x,x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0

    #just slicing 2D Omega
    Omega1Dx = Omega_2D_xy_griddat.sliced([None,yslice])
    
    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)


    Omega1 = np.array(Omega1Dx.data_xyz)
   
    ax1.plot(xcoord, Omega1, color='blue')


    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$\Omega$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Omega_xpm_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ OMEGA 1D in y #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    xslice = 0.0

    #just slicing 2D Omega
    Omega1Dx = Omega_2D_xy_griddat.sliced([xslice,None])


    print("here")
    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)


    Omega1 = np.array(Omega1Dx.data_xyz)

    # ------------------------------
    Om1y_df = pd.DataFrame(Omega1)
    Filename1 = "./Data_extract/ZSaved_Omega_y" + EoS_name + tms_name + ".txt"
    #Om1y_df.to_csv(Filename1,header = False, index= False)
    #------------------------------------

    ax1.plot(xcoord, Omega1, color='blue')

    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$\Omega$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Omega_y_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ Jmom 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    yslice = 0.0
    #just slicing 2D Omega
    Jmom1Dx = Jmom_2D_xy_griddat.sliced([None,yslice])

    print("here")
    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)

    Data_type = object

    Jmom1 = np.array(Jmom1Dx.data_xyz)

    ax1.plot(xcoord, Jmom1, color='blue')

    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$x \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$j$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Jmom_x_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################ Jmom 1D in x #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax1 = fig.add_subplot()

    xslice = 0.0
    #just slicing 2D Omega
    Jmom1Dy = Jmom_2D_xy_griddat.sliced([xslice,None])


    #print(rho1Dx.data_xyz)
    #print(rho1Dx.shape)
    print("here")
    # ---> make 1D grid data to a array which can be plotted ----
    xcoord = np.linspace(lowbound,upbound,NGRID)
    print(xcoord.shape)


    Jmom1 = np.array(Jmom1Dy.data_xyz)

    #--------------------------------
    Jmom1y_df = pd.DataFrame(Jmom1)
    Filename1 = "./Data_extract/ZSaved_Jmom_y" + EoS_name + tms_name + ".txt"
    #Jmom1y_df.to_csv(Filename1,header = False, index= False)
    #-----------------------------------


    ax1.plot(xcoord, Jmom1, color='blue')

    #if(t_curr_int1==50):
    #    ax1.plot(data[0][:,0], data[0][:,2], color='red')
    #if(t_curr_int1==20):
    #    ax1.plot(data[1][:,0], data[1][:,2], color='red')


    ax1.set_xlim(0.,upbound)

    ax1.set_xlabel(r'$y \,[\rm km]$',fontsize=22, style = 'italic')
    ax1.set_ylabel(r'$j$',fontsize=22, style = 'italic')


    plt.plot()
    Filename = "./Plots/Jmom_y_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()


    ################################ OMEGA 1D in Jmom #################################

    fig = plt.figure(figsize=(6,6), dpi=128)
    ax2 = fig.add_subplot()

    yslice = 0.0
    #just slicing 2D Omega
    Omega1Dx_new = Omega_2D_xy_griddat.sliced([None,yslice])
    Jmom1Dx = Jmom_2D_xy_griddat.sliced([None,yslice])

    Omega1_new = np.array(Omega1Dx_new.data_xyz)
    Jmom1 = np.array(Jmom1Dx.data_xyz)


    ## --- trying to cut data in positiive x and negative x
    Nhalf = int(NGRID/2)
    Omega1_positive = Omega1_new[Nhalf:NGRID]
    #Omega1_negative = Omega1_new[0:Nhalf]

    #Omega1_positive =[]
    #for i in range(0,NGRID):
    #    if(i%2==0):
    #        Omega1_positive.append(Omega1_new[i])

    Omega1_pos_arr = []
    for ii in range(0,Nhalf-1):
        Omega1_pos_arr.append(Omega1_positive[ii])
        OmInterpol = (Omega1_positive[ii] + Omega1_positive[ii+1])/2
        Omega1_pos_arr.append(OmInterpol)
    #    #if(i<=Nhalf-1):
    #    #    OmInterpol = (Omega1_positive[i] + Omega1_positive[i+1])/2
    #    #    Omega1_pos_arr.append(OmInterpol)
    #    #if(i==Nhalf):
    #    #    Omega1_pos_arr.append(Omega1_positive[i])
    Omega1_pos_arr.append(Omega1_positive[Nhalf-1])
    if(NGRID%2==0):
        Omega1_pos_arr.append(Omega1_positive[Nhalf-1])

    #----------------------------

    print(len(Omega1_pos_arr))

    ax2.plot(Jmom1, Omega1_new, color='blue')
    #ax2.plot(Jmom1, Omega1_pos_arr, color='blue')
    #ax2.plot(data[0][:,2], data[0][:,3], color='red')

    #if(t_curr_int1==50):
    #    ax1.plot(data[0][:,2], data[0][:,3], color='red')
    #if(t_curr_int1==20):
    #    ax1.plot(data[1][:,2], data[1][:,3], color='red')

    #ax2.set_xlim(0.,upbound)

    ax2.set_xlabel(r'$j$',fontsize=22, style = 'italic')
    ax2.set_ylabel(r'$\Omega$',fontsize=22, style = 'italic')

    # ------------------------------
    Om1x_df = pd.DataFrame(Omega1_new)
    Filename1 = "./Data_extract/ZSaved_Omega_x" + EoS_name + tms_name + ".txt"
   # Om1x_df.to_csv(Filename1,header = False, index= False)

    Jmom1x_df = pd.DataFrame(Jmom1)
    Filename2 = "./Data_extract/ZSaved_Jmom_x" + EoS_name + tms_name + ".txt"
    #Jmom1x_df.to_csv(Filename2,header = False, index= False)
    #-----------------------------------

    plt.plot()
    Filename = "./Plots/Omega_Jmom_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    if(Flag_no_2D_Plots == 1):
        exit(1)
    ###############################################################################
    ######################################## RHO  2D ###################################
    fig = plt.figure(figsize=(6,6), dpi=128)
    ax = fig.add_subplot()
     #ax = fig.add_subplot(projection='3d', computed_zorder=True)
     
    #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
    pic1 = ax.imshow(rho_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
    plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
    cbar.set_label(r'$\rho $',fontsize=18, style = 'italic',fontweight='bold',color="black")

    plt.plot()
    Filename = "./Plots/Rho_xy_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

        ###############################################################################
        ######################################## Press  2D ###################################
    if(Flag_press==1):    
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax = fig.add_subplot()
         #ax = fig.add_subplot(projection='3d', computed_zorder=True)
         
        #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
        pic1 = ax.imshow(press_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
        plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
        cbar.set_label(r'$P $',fontsize=18, style = 'italic',fontweight='bold',color="black")

        plt.plot()
        Filename = "./Plots/Press_xy_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()

        ###############################################################################
        ######################################## Temp 2D ###################################
    if(Flag_temp==1):    
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax = fig.add_subplot()
         #ax = fig.add_subplot(projection='3d', computed_zorder=True)
         
        #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
        pic1 = ax.imshow(temp_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
        plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
        cbar.set_label(r'$T [MeV] $',fontsize=18, style = 'italic',fontweight='bold',color="black")

        plt.plot()
        Filename = "./Plots/Temp_xy_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()
        
        ###############################################################################
        ######################################## Entropy 2D ###################################
    if(Flag_entropy==1):    
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax = fig.add_subplot()
         #ax = fig.add_subplot(projection='3d', computed_zorder=True)
         
        #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
        pic1 = ax.imshow(entropy_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
        plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
        cbar.set_label(r'$S $',fontsize=18, style = 'italic',fontweight='bold',color="black")

        plt.plot()
        Filename = "./Plots/Entropy_xy_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()
        
        ###############################################################################
        ######################################## Ye 2D ###################################
    if(Flag_Ye==1):    
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax = fig.add_subplot()
         #ax = fig.add_subplot(projection='3d', computed_zorder=True)
         
        #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
        pic1 = ax.imshow(Ye_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
        plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
        cbar.set_label(r'$Y_e $',fontsize=18, style = 'italic',fontweight='bold',color="black")

        plt.plot()
        Filename = "./Plots/Y_e_xy_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()
        
        ###############################################################################
        ######################################## Enue 2D ###################################
    if(M1_scheme==1):  
        if(Flag_M1_Enu==1):  
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax = fig.add_subplot()
             #ax = fig.add_subplot(projection='3d', computed_zorder=True)
             
            #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
            pic1 = ax.imshow(Enue_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))
    
            ax.get_yaxis().set_tick_params(which='both', direction='in')
            ax.get_xaxis().set_tick_params(which='both', direction='in')
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
            plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
            cbar.set_label(r'$E_{\nu,e} $',fontsize=18, style = 'italic',fontweight='bold',color="black")
    
            plt.plot()
            Filename = "./Plots/Enue_xy_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()
        
        ###############################################################################
        ######################################## Enua 2D ###################################
    if(M1_scheme==1): 
        if(Flag_M1_Enu==1):     
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax = fig.add_subplot()
             #ax = fig.add_subplot(projection='3d', computed_zorder=True)
             
            #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
            pic1 = ax.imshow(Enua_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))
    
            ax.get_yaxis().set_tick_params(which='both', direction='in')
            ax.get_xaxis().set_tick_params(which='both', direction='in')
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
            plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
            cbar.set_label(r'$E_{\nu,a} $',fontsize=18, style = 'italic',fontweight='bold',color="black")
    
            plt.plot()
            Filename = "./Plots/Enua_xy_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()

        ###############################################################################
        ######################################## Enux 2D ###################################
    if(M1_scheme==1): 
        if(Flag_M1_Enu==1):      
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax = fig.add_subplot()
             #ax = fig.add_subplot(projection='3d', computed_zorder=True)
             
            #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
            pic1 = ax.imshow(Enux_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))
    
            ax.get_yaxis().set_tick_params(which='both', direction='in')
            ax.get_xaxis().set_tick_params(which='both', direction='in')
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
            plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
            cbar.set_label(r'$E_{\nu,x} $',fontsize=18, style = 'italic',fontweight='bold',color="black")
    
            plt.plot()
            Filename = "./Plots/Enux_xy_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()

        ###############################################################################
        ######################################## Eps_nue 2D ###################################
    if(M1_scheme==1):  
        if(Flag_M1_eps_nu==1):    
          
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax = fig.add_subplot()
             #ax = fig.add_subplot(projection='3d', computed_zorder=True)
             
            #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
            pic1 = ax.imshow(eps_nue_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))
    
            ax.get_yaxis().set_tick_params(which='both', direction='in')
            ax.get_xaxis().set_tick_params(which='both', direction='in')
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
            plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
            cbar.set_label(r'$\epsilon_{\nu,e} $',fontsize=18, style = 'italic',fontweight='bold',color="black")
    
            plt.plot()
            Filename = "./Plots/Eps_nue_xy_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()

        ###############################################################################
        ######################################## Eps_nua 2D ###################################
    if(M1_scheme==1):   
        if(Flag_M1_eps_nu==1):  
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax = fig.add_subplot()
             #ax = fig.add_subplot(projection='3d', computed_zorder=True)
             
            #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
            pic1 = ax.imshow(eps_nua_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))
    
            ax.get_yaxis().set_tick_params(which='both', direction='in')
            ax.get_xaxis().set_tick_params(which='both', direction='in')
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
            plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
            cbar.set_label(r'$\epsilon_{\nu,a} $',fontsize=18, style = 'italic',fontweight='bold',color="black")
    
            plt.plot()
            Filename = "./Plots/Eps_nua_xy_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()
      
        ###############################################################################
        ######################################## Eps_nux 2D ###################################
    if(M1_scheme==1):  
        if(Flag_M1_eps_nu==1):   
            fig = plt.figure(figsize=(6,6), dpi=128)
            ax = fig.add_subplot()
             #ax = fig.add_subplot(projection='3d', computed_zorder=True)
             
            #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
            pic1 = ax.imshow(eps_nux_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))
    
            ax.get_yaxis().set_tick_params(which='both', direction='in')
            ax.get_xaxis().set_tick_params(which='both', direction='in')
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
            plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
            cbar.set_label(r'$\epsilon_{\nu,x} $',fontsize=18, style = 'italic',fontweight='bold',color="black")
    
            plt.plot()
            Filename = "./Plots/Eps_nux_xy_it_" + str(it) + ".jpg"
            if(Flag_Save_Plots ==1):
                plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
            plt.show()
 
        ###############################################################################
        ######################################## Bvec0 2D ###################################
    if(magfields==1):    
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax = fig.add_subplot()
         #ax = fig.add_subplot(projection='3d', computed_zorder=True)
         
        #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
        pic1 = ax.imshow(Bvec0_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
        plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
        cbar.set_label(r'$B_{0} $',fontsize=18, style = 'italic',fontweight='bold',color="black")

        plt.plot()
        Filename = "./Plots/Bvec0_xy_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()
 
        ###############################################################################
        ######################################## Bvec1 2D ###################################
    if(magfields==1):    
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax = fig.add_subplot()
         #ax = fig.add_subplot(projection='3d', computed_zorder=True)
         
        #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
        pic1 = ax.imshow(Bvec1_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
        plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
        cbar.set_label(r'$B_{1} $',fontsize=18, style = 'italic',fontweight='bold',color="black")

        plt.plot()
        Filename = "./Plots/Bvec1_xy_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()
 
        ###############################################################################
        ######################################## Bvec2 2D ###################################
    if(magfields==1):    
        fig = plt.figure(figsize=(6,6), dpi=128)
        ax = fig.add_subplot()
         #ax = fig.add_subplot(projection='3d', computed_zorder=True)
         
        #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
        pic1 = ax.imshow(Bvec2_2D_xy.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))# extent=(-15,15,-15,15))

        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
        plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
        cbar.set_label(r'$B_{2} $',fontsize=18, style = 'italic',fontweight='bold',color="black")

        plt.plot()
        Filename = "./Plots/Bvec2_xy_it_" + str(it) + ".jpg"
        if(Flag_Save_Plots ==1):
            plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
        plt.show()
                            
                                    
    ################################################################
    ###############################################################################
    ######################################## W lorentz 2D ###################################
    fig = plt.figure(figsize=(6,6), dpi=128)
    ax = fig.add_subplot()
     #ax = fig.add_subplot(projection='3d', computed_zorder=True)
     
    pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower",  extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))#extent=(-15,15,-15,15))
    #pic1 = ax.imshow(Omega_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
    plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
    cbar.set_label(r'$W $',fontsize=18, style = 'italic',fontweight='bold',color="black")

    plt.plot()
    Filename = "./Plots/Wlor_xy_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################################################
    ###############################################################################
    ######################################## Vphi 2D ###################################
    fig = plt.figure(figsize=(6,6), dpi=128)
    ax = fig.add_subplot()
     #ax = fig.add_subplot(projection='3d', computed_zorder=True)
     
    #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
    pic1 = ax.imshow(vphi_2D_xy_griddat.data_xyz, origin="lower", vmin = 0.0, vmax = 0.06,  extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))#extent=(-15,15,-15,15))

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
    plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
    cbar.set_label(r'$v_{\phi} $',fontsize=18, style = 'italic',fontweight='bold',color="black")
 
    plt.plot()
    Filename = "./Plots/Vphi_xy_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################################################

    ###############################################################################
    ######################################## JMOM 2D ###################################
    fig = plt.figure(figsize=(6,6), dpi=128)
    ax = fig.add_subplot()
     #ax = fig.add_subplot(projection='3d', computed_zorder=True)
     
    #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
    pic1 = ax.imshow(Jmom_2D_xy_griddat.data_xyz, origin="lower", extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))#extent=(-15,15,-15,15))

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
    plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
    cbar.set_label(r'$j $',fontsize=18, style = 'italic',fontweight='bold',color="black")
 
    plt.plot()
    Filename = "./Plots/Jmom_xy_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################################################

    ###############################################################################
    ######################################## OMEGA 2D ###################################
    fig = plt.figure(figsize=(6,6), dpi=128)
    ax = fig.add_subplot()
     #ax = fig.add_subplot(projection='3d', computed_zorder=True)
     
    #pic1 = ax.imshow(Wlor_2D_xy.data_xyz, origin="lower", extent=(-15,15,-15,15))
    pic1 = ax.imshow(Omega_2D_xy_griddat.data_xyz, origin="lower", vmin = 0.0, vmax = 0.06,  extent=(lowbound2D,upbound2D,lowbound2D,upbound2D))#extent=(-15,15,-15,15))

    # ---- > plot rho_2D at given slice : -------------------
    #zslice=7.7 #0.8 #0.0#5.0 #10.0 #2.0
    #rho_2D=rho_3D.sliced([None,None,zslice])
    #ax.imshow(rho_2D.data_xyz, origin="lower")


    # ----> plot rho_2D as 2D grid data : -------------------
    #ax.imshow(rho_2D_xz.data_xyz, origin="lower")
    #ax.imshow(rho_2D_xy.data_xyz, origin="lower")



    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_xticklabels(),rotation=0,fontsize=nFontSize)
    plt.setp(ax.get_yticklabels(),rotation=0,fontsize=nFontSize)
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
    cbar.set_label(r'$\Omega $',fontsize=18, style = 'italic',fontweight='bold',color="black")


    # ---> remove axes : -------------------------
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)

    # X AXIS -BORDER
    #ax.spines['bottom'].set_visible(False)
    # pad of xlabels
    #ax.set_xticklabels([])
    # xticks
    #ax.set_xticks([])
    # xlabels and ticks
     #  ax1.axes.get_xaxis().set_visible(False)

    # Y AXIS -BORDER
    #ax.spines['left'].set_visible(False)
    # pad of y label
    #ax.set_yticklabels([])
    # yticks
    #ax.set_yticks([])
    #ax.set_frame_on(False)
       
    plt.plot()
    Filename = "./Plots/Omega_xy_it_" + str(it) + ".jpg"
    if(Flag_Save_Plots ==1):
        plt.savefig(Filename, dpi=400, pad_inches=0,bbox_inches="tight",format="jpg")
    plt.show()

    ################################################################
    i +=1
    continue
    ################################################################  
print("I guess I ended again")
