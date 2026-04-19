!================================================================================
!
!    BHAC (The Black Hole Accretion Code) solves the equations of
!    general relativistic magnetohydrodynamics and other hyperbolic systems
!    in curved spacetimes.
!
!    Copyright (C) 2019 Oliver Porth, Hector Olivares, Yosuke Mizuno, Ziri Younsi,
!    Luciano Rezzolla, Elias Most, Bart Ripperda and Fabio Bacchini
!
!    This file is part of BHAC.
!
!    BHAC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BHAC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BHAC.  If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

!=============================================================================
! include file amrvacdef.t
!use mpi ! Well, this seemed like a good idea but does not compile on iboga.
{#IFDEF UNIT_TESTS
module amrvacdef
}
use mod_indices
use mod_physicaldata
use mod_connectivity
use amrvacpar
{#IFDEF UNIT_TESTS
 use mod_m1_tests
}
IMPLICIT NONE

! DEFINITIONS OF GLOBAL PARAMETERS AND VARIABLES
! Parameters:

! commented r_,z_,pphi_zz_ here as it is defined in mod_physicaldata

! Indices for cylindrical coordinates FOR TESTS, negative value when not used:
INTEGER,PARAMETER:: r_=1, phi_=^PHI, z_=^Z

! Indices for cylindrical coordinates FOR INDEXING, always positive
{^IFTHREED
!INTEGER,PARAMETER:: pphi_=^PPHI, zz_=^ZZ
}

!include 'amrvacpar.f'
{#IFNDEF UNIT_TESTS
INTEGER,PARAMETER:: ndim=^ND, ndir=^NC
}

include 'amrvacsettings.f'

INTEGER, PARAMETER :: ixGlo^D=1
COMMON, INTEGER :: ixGhi^D
COMMON, INTEGER :: block_nx^D
INTEGER,PARAMETER:: filelog_=1,fileout_=2,fileslice_=3,filecollapse_=4,fileanalysis_=5,fileshell_=6,fileofdet_=7,nfile=7 ! outputfiles
INTEGER,PARAMETER:: nslicemax=1000, nshellmax=1000

INTEGER,PARAMETER:: unitstdin=5,unitterm=6,uniterr=6 ! Unit names.

! Units reserved for files:
INTEGER,PARAMETER:: unitpar=9
INTEGER,PARAMETER:: unitconvert=10
INTEGER,PARAMETER:: unitslice=11
INTEGER,PARAMETER:: unitsnapshot=12
INTEGER,PARAMETER:: unitcollapse=13
INTEGER,PARAMETER:: unitanalysis=14
INTEGER,PARAMETER:: unitshell=15

INTEGER,PARAMETER:: biginteger=10000000

! Note: smalldouble must be above machine precision 
DOUBLE PRECISION,PARAMETER:: smalldouble=1.D-16, bigdouble=1.D+99
DOUBLE PRECISION,PARAMETER:: zero=0D0,one=1D0,two=2D0,half=0.5D0,quarter=0.25D0,third=0.33333333333333333333d0
DOUBLE PRECISION,PARAMETER:: dpi=3.141592653589793238462643383279502884197169399375105d0
DOUBLE PRECISION,PARAMETER:: de=2.71828182845904523536028747135

! Physical scaling parameters:
COMMON, DOUBLE PRECISION:: UNIT_LENGTH, UNIT_DENSITY, UNIT_VELOCITY

{#IFDEF UNIT_TESTS
INTEGER,PARAMETER:: swap_=4+1, nspecialpar=1
}
{#IFNDEF UNIT_TESTS
include 'amrvacusrpar.f'
}

! For transform variables and save selected data
COMMON, INTEGER :: nwtf
COMMON, INTEGER :: neqpartf

!Kronecker delta and Levi-Civita tensors
COMMON, INTEGER:: kr(0:3,0:3),lvc(1:3,1:3,1:3)

!Equation and method parameters
COMMON, DOUBLE PRECISION:: eqpar(neqpar+nspecialpar)
COMMON, DOUBLE PRECISION:: time_between_print
! Time step control parameters
COMMON, DOUBLE PRECISION :: courantpar, dtpar, dtdiffpar, dtTCpar
COMMON, CHARACTER*131 :: typecourant,typeresid
COMMON, LOGICAL :: time_accurate, addmpibarrier, use_t_dependent_output, extra_2d_output
COMMON, logical :: PPM_extrema, PPM_flatcd, PPM_flatsh, PPM_flatten
COMMON, logical :: use_4th_order_fd, use_6th_order_fd, use_del_zanna

! restart
COMMON, CHARACTER*131 :: restart_from_file 
COMMON, LOGICAL :: need_restart_output

!Time parameters
INTEGER,PARAMETER:: nsavehi=100       ! maximum No. saves into outputfiles
                                      ! defined by arrays of tsave or itsave
COMMON, DOUBLE PRECISION:: t,tmax,continue_t,dtmin,tin,residmin,residmax,residual
COMMON, DOUBLE PRECISION:: tfixgrid
COMMON, DOUBLE PRECISION:: tsave(nsavehi,nfile),tsavelast(nfile),dtsave(nfile),slicecoord(nslicemax),shellcoord(nshellmax)
COMMON, DOUBLE PRECISION:: tsavelast_rt(nfile),dtsave_rt(nfile)
COMMON, LOGICAL:: tmaxexact,treset,itreset,firstprocess,resetgrid,fixprocess,changeglobals,collapse(ndim)
COMMON, INTEGER:: it_start,it,itmax,itmin,continue_it,slowsteps{#IFDEF MAGNETOFRICTION , mfitmax}
COMMON, INTEGER:: itsave(nsavehi,nfile),itsavelast(nfile),ditsave(nfile)
COMMON, INTEGER:: isavet(nfile),isaveit(nfile), nslices, nshells, slicedir(nslicemax), collapseLevel{^NOONED, nxShell^DE}
COMMON, INTEGER:: n_saves(1:nfile)
COMMON, INTEGER:: n_saves_rt(1:nfile)
COMMON, INTEGER:: typeparIO
COMMON, INTEGER:: itfixgrid,ditregrid
COMMON, INTEGER:: nwauxio
COMMON, INTEGER:: istep, nstep

!Method switches
COMMON, CHARACTER*131 :: typeadvance
COMMON, CHARACTER*131 :: typelow1(nlevelshi),typelimited,typesourcesplit
COMMON, CHARACTER*131 :: typefull1(nlevelshi), typepred1(nlevelshi)
COMMON, CHARACTER*131 :: typelimiter1(nlevelshi),typegradlimiter1(nlevelshi)
COMMON, INTEGER       :: type_limiter(nlevelshi), type_gradient_limiter(nlevelshi)
COMMON, CHARACTER*131 :: typetvd,typeaverage
COMMON, LOGICAL       :: PP_limiter
COMMON, INTEGER       :: fv_n_interp
COMMON, LOGICAL       :: ssplitdivb,ssplitresis,ssplituser,useprimitive,dimsplit,fixp_after_reconstruct
COMMON, LOGICAL       :: old_bhac_safety, use_backup_c2p
COMMON, LOGICAL       :: reconstruct_eps_for_T
COMMON, INTEGER       :: face_to_center_order
!character*131, allocatable :: typeentropy(:)

COMMON, CHARACTER*131 :: typedimsplit,typeaxial,typecoord,typepoly,typeinversion,typeinversionresis
COMMON, INTEGER:: errorestimate,nxdiffusehllc,typespherical,ncyclemax
COMMON, CHARACTER*131 :: typeprolonglimit, clean_init_divB
INTEGER       :: typelimiter, typegradlimiter
!double precision, allocatable :: entropycoef(:)

COMMON, DOUBLE PRECISION:: tvdlfeps, mcbeta, parastsnu, TCphi
COMMON, LOGICAL:: sourceparasts,sourceimpl,sourceimplcycle,conduction,TCsaturate,energyonly
COMMON, LOGICAL:: flathllc,flatcd,flatsh

!logical, allocatable :: loglimit(:), logflag(:)

COMMON, LOGICAL:: restrictprimitive,prolongprimitive, &
                  coarsenprimitive,useprimitiveRel, &
                  amrentropy
COMMON, LOGICAL:: extra_buffer_refine

COMMON, LOGICAL:: BnormLF
COMMON, DOUBLE PRECISION:: smallT,smallp,smallrho,amr_wavefilter(nlevelshi)
COMMON, CHARACTER*131 :: typediv,typegrad,typeemf,DYSP_solver_type

COMMON, LOGICAL:: strictnr,strictgetaux, evolve_hydro, metric_vars_interpolation
COMMON, DOUBLE PRECISION::dmaxvel,tolernr,absaccnr,tlow,lfac_max
COMMON, INTEGER:: maxitnr,nflatgetaux

COMMON, LOGICAL:: convert_nocartesian, slice_nocartesian
COMMON, LOGICAL:: writelevel(nlevelshi)

!logical, allocatable :: writew(:)

COMMON, DOUBLE PRECISION:: writespshift(ndim,2)
COMMON, INTEGER:: level_io, level_io_min, level_io_max

! local and global fastest wave speed (computed in setdt):
COMMON, DOUBLE PRECISION :: cmax_mype, cmax_global

!Boundary region parameters
INTEGER,PARAMETER:: nhiB=2*ndim         ! maximum No. boundary sections
COMMON, LOGICAL:: periodB(ndim), poleB(2,ndim), aperiodB(ndim), primitiveB(2,ndim)

!character*131, allocatable :: typeB(:,:)

COMMON, CHARACTER*131 :: typeghostfill,typegridfill
COMMON, DOUBLE PRECISION::ratebdflux
COMMON, LOGICAL:: internalboundary

!File parameters
COMMON, CHARACTER(len=131) :: inifile,filenameout,filenameini,filenamelog
COMMON, CHARACTER*131 :: fileheadout
COMMON, CHARACTER*1024 :: wnames,primnames,wnameslog
COMMON, CHARACTER*131 :: typefilelog{#IFDEF HDF5 ,xdmf_type \}

COMMON, character*131 :: small_values_method

COMMON, INTEGER :: snapshotini
COMMON, LOGICAL :: sliceascii{#IFDEF HDF5 ,hdf5_ini,fastIO,write_xdmf,save_GZ \}

! healpix spherical outflow detector params
INTEGER, PARAMETER :: max_healpix_det_num = 20
COMMON, INTEGER    :: healpix_var_num, healpix_nside(1:max_healpix_det_num), healpix_det_num
COMMON, LOGICAL    :: healpix_det_open, healpix_det_zsymm
COMMON, CHARACTER*2048 :: healpix_vname_str
COMMON, DOUBLE PRECISION :: healpix_det_radii(1:max_healpix_det_num)

!M1 neutrino radiation params
COMMON, character*512 :: fileWeakhub
COMMON, character*131 :: m1_closure_type
COMMON, DOUBLE PRECISION :: m1_E_atmo, m1_N_atmo, TESTqtC, m1_tset, m1_tset_backreact, m1_rho_floor
COMMON, LOGICAL    :: m1_radice_speeds,m1_actual_speeds,m1_parabolic,&
                      m1_frad_A2_A3,m1_frad_A2_A,m1_frad_LLF,m1_erad_LLF, &
                      M1_FLUID_BACKREACT,m1_use_neutrinos,m1_use_muons,m1_use_photons,m1_2_eas_updates
!!m1_i_nue, m1_i_nuebar, m1_i_nux, m1_i_mu, m1_i_mubar, m1_i_photon

!Convert parameters
COMMON, LOGICAL :: convert,autoconvert,saveprim,uselimiter,endian_swap
COMMON, CHARACTER*131 :: convert_type, dxfiletype, collapse_type, slice_type, shell_type
COMMON, DOUBLE PRECISION :: normt
!double precision, allocatable :: normvar(:)

! --------------------------------------------
INTEGER:: saveigrid
! Stores the memory and load imbalance, to be used in printlog:
COMMON, DOUBLE PRECISION :: Xload, Xmemory

COMMON, DOUBLE PRECISION:: time_bc

integer,parameter:: nodehi=^ND+1
integer,parameter:: plevel_=1
integer,parameter:: pig^D_=plevel_+^D

integer,parameter:: rnodehi=3*^ND
integer,parameter:: rpxmin0_=0
integer,parameter:: rpxmin^D_=rpxmin0_+^D 
integer,parameter:: rpxmax0_=^ND
integer,parameter:: rpxmax^D_=rpxmax0_+^D 
integer,parameter:: rpdx^D_=2*^ND+^D

! parameters for bc_phys
integer,parameter:: ismin^D=-1+2*^D
integer,parameter:: ismax^D=2*^D



{#IFNDEF UNIT_TESTS
include 'mpif.h'
}
!-----------------------------------------------------------------------------
! common block variables
!-----------------------------------------------------------------------------

! nodal info per grid: 
! rnode:            corner coordinates,dx(idim)
! node:             dimensions and pointers 

! mxnest:           maximal number of levels

! tol:              error tolerance used for refinement decision
! nbufferx D:       number of cells as buffer zone

! dixB:             number of ghost cells surrounding a grid for boundary cndt.

double precision :: tol, tolratio
double precision :: rnode, rnode_sub, dx, dt, dtimpl, dt_grid, dxlevel, dt_old, dt_old2, dt_old3
integer          :: mxnest, dixB, nbufferx^D, nxlone^D
integer          :: node, node_sub, levmin, levmax, levmax_sub
logical          :: skipfinestep, time_advance
common /nodalr/     rnode(rnodehi,ngridshi), rnode_sub(rnodehi,ngridshi), tol(nlevelshi), tolratio(nlevelshi), &
                    dx(ndim,nlevelshi), dt, dtimpl, dt_grid(ngridshi), dt_old, dt_old2, dt_old3
common /nodali/     node(nodehi,ngridshi), node_sub(nodehi,ngridshi), &
                    nbufferx^D, mxnest, dixB, levmin, levmax, levmax_sub, nxlone^D
common /nodall/     skipfinestep, time_advance

common /ompdoub/ dxlevel(ndim)
common /ompinte/ saveigrid, typelimiter, typegradlimiter

! from amrvacpar originally
COMMON, DOUBLE PRECISION::minp,minrho,smallxi,smalltau{#IFNDEF SYNGE , govergminone}
COMMON, DOUBLE PRECISION::limitvalue

! xprob: problem box; iprob: problem
COMMON, INTEGER:: iprob
COMMON, DOUBLE PRECISION:: xprob^L

! Iij(1:6)
COMMON, double precision:: Iij_dot_old(1:6), Iij_dot_old2(1:6)
COMMON, double precision:: Iij_old(1:6), Iij_old2(1:6), Iij_old3(1:6)
COMMON, double precision:: I3ij_new_global(1:6) 
!COMMON, DOUBLE PRECISION:: Iij_dot_old(1:6), Iij_dot_old2(1:6)
!COMMON, DOUBLE PRECISION:: Iij_old(1:6), Iij_old2(1:6), Iij_old3(1:6)
!COMMON, DOUBLE PRECISION:: I3ij_new_global(1:6) 

!$ OMP THREADPRIVATE(/ompdoub/,/ompinte/)
!  end include file amrvacdef.t
{#IFDEF UNIT_TESTS
end module amrvacdef
}
!============================================================================
