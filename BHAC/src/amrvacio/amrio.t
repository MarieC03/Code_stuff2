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
subroutine readcommandline

!use M_kracken
include 'amrvacdef.f'
integer           :: len, ier
logical           :: help

INTEGER :: i,stat
CHARACTER(len=32) :: arg
LOGICAL :: unknown_arg

!----------------------------------------------------------------------------

! =================== Kracken command line reading =================
!
!!defaults and usage:
!call kracken('cmd','-i amrvac.par -if data -restart -1 -slice 0 -collapse 0 '//&
!  '-shell 0 --help .false. -convert .false.')

!! Getting the filename
!      call retrev('cmd_i',inifile,len,ier)
!      call retrev('cmd_if',filenameini,len,ier)
!      snapshotini = iget('cmd_restart')
!      snapshotnext= snapshotini+1
!      slicenext   = iget('cmd_slice')
!      collapseNext   = iget('cmd_collapse')
!      shellNext      = iget('cmd_shell')
!      help = lget('cmd_-help')                    ! --help present?
!      convert = lget('cmd_convert')               ! -convert present?
!
! ==================================================================

! =============== Fortran 2003 command line reading ================

  ! Default command line arguments

  inifile="amrvac.par"

  filenameini="data"

  snapshotini=-1

  snapshotnext=0

  slicenext=0

  collapsenext=0

  shellnext=0

  help=.false.

  convert=.false.

  ! Argument 0 is program name, so we start from one
  i = 1

  unknown_arg=.false.

  DO
    CALL get_command_argument(i, arg)

    IF (LEN_TRIM(arg) == 0) EXIT

    select case(arg)
    case("-i")

        i = i+1
    
        CALL get_command_argument(i, arg)
    
        inifile=TRIM(arg)

    case("-if")

        i = i+1
    
        CALL get_command_argument(i, arg)
    
        filenameini=TRIM(arg)    

    case("-restart")

        i = i+1
    
        CALL get_command_argument(i, arg)
   
        read(arg,*,iostat=stat) snapshotini

        snapshotnext = snapshotini+1
 
    case("-slice")

        i = i+1
    
        CALL get_command_argument(i, arg)
   
        read(arg,*,iostat=stat) slicenext


    case("-collapse")

        i = i+1
    
        CALL get_command_argument(i, arg)
   
        read(arg,*,iostat=stat) collapsenext

    case("-shell")

        i = i+1
    
        CALL get_command_argument(i, arg)
   
        read(arg,*,iostat=stat) shellnext

    case("-convert")

        convert=.true.

    case("--help")

        help=.true.

        EXIT

    case default

        unknown_arg=.true.

        help=.true.

        EXIT

    end select

    i = i+1

  END DO

  if (unknown_arg) then
       print*,"======================================="
       print*,"Error: Command line argument ' ",TRIM(arg)," ' not recognized"
       print*,"======================================="

        help=.true.

  end if

! ==================================================================


if (mype==0) then
   print*,'-----------------------------------------------------------------------------'
   print*,'-----------------------------------------------------------------------------'
   print*,'                         ____  _    _          _____                         '
   print*,'                        |  _ \| |  | |   /\   / ____|                        '
   print*,'                        | |_) | |__| |  /  \ | |                             '
   print*,'                        |  _ <|  __  | / /\ \| |                             '
   print*,'                        | |_) | |  | |/ ____ \ |____                         '
   print*,'                        |____/|_|  |_/_/    \_\_____|                        '
   print*,'                        -----------------------------                        '
   print*,'                        The Black Hole Accretion Code                        '
   print*,'-----------------------------------------------------------------------------'
   print*,'-----------------------------------------------------------------------------'
end if

if (help) then
   if (mype==0) then 
      print*,'calling example:                                         '
      print*,'./bhac -i parameterfile -restart 100 -convert -if datamr/data'
      print*,'default parameterfile is amrvac.par                            '
      print*,'Note that parameterfile parameters overwrite the commandline   '
      print*,'-----------------------------------------------------------------------------'
      print*,'Available options are:'
      print*,'-i parfile'
      print*,'-if filenameout'
      print*,'-restart #snapshot'
      print*,'-convert'
      print*,'-slice #slicenext'
      print*,'-collapse #collapsenext'
      print*,'-shell #shellnext'
      print*,'--help 	Prints the help message'
      print*,'-----------------------------------------------------------------------------' 
      print*,'See COMMAND LINE PARAMETERS entry in the manual for more information'
      print*,'(https://bhac.science/documentation_html/commandline.html).'
      print*,'                                                         '
   endif
   call comm_finalize
   STOP
endif

end subroutine readcommandline
!=============================================================================
subroutine readparameters
use mod_imhd_con2prim
use mod_limiter
use mod_cfc_parameters
use mod_cfc, only: cfc_solver_activate
use mod_gw_br
use mod_eos_tabulated_parameters
use mod_eos_tabulated
use mod_eos_idealgas
use mod_eos_hybrid
use mod_eos_polytrope
use mod_eos
use mod_variables
include 'amrvacdef.f'

logical :: fileopen
integer :: ifile, iB, isave, iw, level, idim, islice
double precision :: dxlone^D
double precision :: cfc_tol1, cfc_tol2, cfc_tol3
double precision :: cfc_tol_evolve1, cfc_tol_evolve2, cfc_tol_evolve3
double precision :: gw_br_tol1, gw_br_tol2, gw_br_tol3
double precision :: gw_br_tol_evolve1, gw_br_tol_evolve2, gw_br_tol_evolve3
character*30     :: cfc_coordinate

namelist /filelist/   filenameout,filenameini,filenamelog, &
                      snapshotini,typefilelog,firstprocess,resetgrid,changeglobals,snapshotnext, &
                      convert,convert_type,slice_type,dxfiletype,saveprim,&
                      primnames,typeparIO,uselimiter,nwauxio,convert_nocartesian,slice_nocartesian,addmpibarrier, &
                      writew,writelevel,writespshift,endian_swap, &
                      normvar,normt,level_io,level_io_min,level_io_max, &
                      autoconvert,slicenext,collapseNext,collapse_type, &
                      shellNext,shell_type{#IFDEF HDF5 ,hdf5_ini,fastIO,write_xdmf,save_GZ,xdmf_type\}
namelist /savelist/   tsave,itsave,dtsave,ditsave,dtsave_rt,nslices,slicedir,slicecoord,collapse,collapseLevel,&
                      nshells,shellcoord,nxShell^DE,time_between_print, use_t_dependent_output, extra_2d_output, &
                      need_restart_output
namelist /stoplist/   itmax,tmax,continue_it,continue_t,tmaxexact,dtmin,t,it,treset,itreset,residmin,residmax,typeresid
namelist /methodlist/ wnames,fileheadout,typeadvance, &
                      ssplitdivb,ssplitresis,ssplituser,typesourcesplit,&
                      sourceimpl,sourceimplcycle,conduction,TCsaturate,TCphi,ncyclemax,sourceparasts,parastsnu,&
                      dimsplit,typedimsplit,typeaxial,typecoord,&
                      typefull1,typepred1,typelow1,&
                      typelimiter1,mcbeta,typegradlimiter1,&
                      flatcd,flatsh,&
                      typelimited,useprimitive, &
                      typetvd,typeentropy,entropycoef,typeaverage, &
                      B0field,Bdip,Bquad,Boct,Busr,&
                      useprimitiveRel,maxitnr,dmaxvel,typepoly,&
                      tvdlfeps,BnormLF,&
                      smallT,smallp,smallrho,typegrad,typediv,tolernr,absaccnr,&
                      strictnr,strictgetaux,nflatgetaux,&
                      nxdiffusehllc,typespherical,&
                      fixprocess,flathllc, &
                      tlow,nwtf,neqpartf,typeinversion,&
                      typeemf,typeinversionresis, clean_init_divB, fv_n_interp, &
                      lfac_max, evolve_hydro, small_values_method, PP_limiter, Kc2p_tolerance, &
                      PPM_extrema,PPM_flatcd,PPM_flatsh,PPM_flatten,usr_W_limit_scheme, old_bhac_safety, &
                      fixp_after_reconstruct, reconstruct_eps_for_T, &
                      use_4th_order_fd, use_6th_order_fd, use_del_zanna, face_to_center_order
namelist /boundlist/  dixB,typeB,typeghostfill,typegridfill,ratebdflux,&
                      internalboundary,primitiveB
namelist /amrlist/    mxnest,nbufferx^D,tol,tolratio,errorestimate, &
                      amr_wavefilter,nxlone^D,dxlone^D,iprob,xprob^L, &
                      skipfinestep,wflags,flags,&
                      restrictprimitive,prolongprimitive,coarsenprimitive, &
                      typeprolonglimit, &
                      amrentropy,logflag,tfixgrid,itfixgrid,ditregrid, &
                      block_nx^D, extra_buffer_refine
namelist /paramlist/  time_accurate, courantpar, dtpar, dtdiffpar, dtTCpar,&
                      typecourant, slowsteps

namelist /cfclist/    cfc_evolve, cfc_smallest_dt, cfc_dit_update, cfc_dt_update, &
                      cfc_tol1, cfc_tol2, cfc_tol3, cfc_tol_evolve1, &
                      cfc_tol_evolve2, cfc_tol_evolve3, initialize_metric, restart_init_metric, &        
                      metric_vars_interpolation, cfc_coordinate, use_lfac_profile_to_init, &
                      cfc_t_end_of_stage1, cfc_t_end_of_stage2, &
                      cfc_dit_update_stage2, cfc_dit_update_stage3, & 
                      cfc_dt_update_stage2, cfc_dt_update_stage3, cfc_beta_robin_bc, cfc_n_cycle, &
                      fix_prim_init, use_check_psi, use_gw_br, &
                      gw_br_tol1, gw_br_tol2, gw_br_tol3, &
                      gw_br_tol_evolve1, gw_br_tol_evolve2, gw_br_tol_evolve3, gw_br_include_dU_i, &
                      gw_br_use_I3_ij, gw_br_use_wi, use_h00_coupling, &
                      gw_br_dt_update, use_index_contract, use_3rd_d_Iij, use_h00_cfc_betai, use_h00_cfc_alp, &
                      gw_br_couple_weakly
namelist /eoslist/    eos_type_input, atmo_type, table_type_input, eos_table_name, &
                      eos_gamma, eos_adiab, eos_gamma_th, massn_cgs, &
                      small_rho, small_rho_thr, small_rho_fac, &
                      small_temp, atmo_gamma, atmo_adiab, use_realistic_mp_table

namelist /outflowlist/  healpix_det_open, healpix_det_num, healpix_det_radii, healpix_det_zsymm, &
                      healpix_var_num, healpix_vname_str, healpix_nside

namelist /m1list/  fileWeakhub, m1_closure_type, m1_actual_speeds, m1_radice_speeds, m1_parabolic, &
                    m1_frad_A2_A, m1_frad_A2_A3, m1_erad_LLF, m1_frad_LLF, M1_FLUID_BACKREACT, &
                    m1_use_neutrinos, m1_use_muons, m1_use_photons, m1_2_eas_updates, TESTqtC, m1_E_atmo, &
                    m1_N_atmo, m1_tset, m1_tset_backreact, m1_rho_floor

! m1_i_nue, m1_i_nuebar, m1_i_nux, m1_i_mu, m1_i_mubar, m1_i_photon
!----------------------------------------------------------------------------

! defaults for boundary treatments
ratebdflux=one
typeghostfill='linear' 
dixB=2
typeB(1:nw,1:nhiB)='cont'
internalboundary=.false.
primitiveB(1:2,1:ndim) = .false.
fv_n_interp = 4
fv_n_interp = min(fv_n_interp, dixB * 2) 

! code behavior and testing defaults
addmpibarrier=.false.

! defaults for specific options
fixprocess=.false.
typegrad='central'
typediv='central'
smallT=-one
smallp=-one
smallrho=-one

! relativistic module defaults
useprimitiveRel=.true.
strictnr=.true.
strictgetaux=.false.
nflatgetaux=1
typepoly='gammie'
tolernr=1.0d-13
absaccnr=1.0d-13
maxitnr=100
dmaxvel=1.0d-8
tlow=zero
fixp_after_reconstruct = .false.
reconstruct_eps_for_T = .false.

! defaults for convert behavior
nwauxio=0
convert_nocartesian=.false.
slice_nocartesian=.false.
saveprim=.false.
!saveprim=.true.
autoconvert=.false.
endian_swap=.false.
convert_type='vtuBCCmpi'
collapse_type='vti'
dxfiletype='lsb'
writew(1:nw)=.true.
writelevel(1:nlevelshi)=.true.
writespshift(1:ndim,1:2)=zero
level_io=-1
level_io_min=1
level_io_max=nlevelshi
use_t_dependent_output = .false.
extra_2d_output = .false.
need_restart_output = .false.

! normalization of primitive variables: only for output
! note that normvar(0) is for length
! this scaling is optional, and must be set consistently if used
normvar(0:nw)=one
normt=one

! residual defaults
residmin=-1.0d0
residmax=bigdouble
typeresid='relative'

! AMR related defaults
{^D& block_nx^D = 4 \}
mxnest=1
nbufferx^D=0;
tol(1:nlevelshi)=0.1d0
tolratio(1:nlevelshi)=1.0d0/8.0d0
typegridfill='linear'
amrentropy=.false.
restrictprimitive=.true.
coarsenprimitive=.true.
prolongprimitive=.true.
typeprolonglimit='default'
extra_buffer_refine = .false. 
errorestimate=3
flags(1:nflag_)=0
wflags(1:nflag_)=zero
flags(nflag_)=1
flags(1)=1
!wflags(1)=one
!code test
wflags(rho_)=one  !using prim instead of D_ to refine as default


logflag(1:nw)=.false.
amr_wavefilter(1:nlevelshi)=1.0d-2
skipfinestep=.false.
tfixgrid=bigdouble
itfixgrid=biginteger
ditregrid=1

! MHD specific defaults
B0field=.false.
Bdip=zero
Bquad=zero
Boct=zero
Busr=zero


! IO defaults
itmax=biginteger
{#IFDEF MAGNETOFRICTION
mfitmax=0
}
continue_it = biginteger
tmax=bigdouble
continue_t = bigdouble 
tmaxexact=.true.
dtmin=1.0d-10
typeparIO=0
nslices=0
collapse=.false.
collapseLevel=1
sliceascii=.true.
slice_type='vtuCC'
{#IFDEF HDF5
hdf5_ini=.true.
fastIO=.true.
write_xdmf=.true.
xdmf_type='Node'
save_GZ=.False.
\}

nshells=0
shell_type='vtu'

do ifile=1,nfile
   do isave=1,nsavehi
      tsave(isave,ifile)=bigdouble   ! t  of saves into the output files
      itsave(isave,ifile)=biginteger ! it of saves into the output files
   end do
   dtsave(ifile)=bigdouble           ! time between saves
   dtsave_rt(ifile)=bigdouble        ! reality time between saves
   ditsave(ifile)=biginteger         ! timesteps between saves
   isavet(ifile)=1                   ! index for saves by t
   isaveit(ifile)=1                  ! index for saves by it
end do
typefilelog='default'
fileheadout = 'AMRVAC'
nwtf=0
neqpartf=0

! defaults for input 
firstprocess=.false.
resetgrid=.false.
changeglobals=.false.
treset=.false.
itreset=.false.
filenameout='data'
filenamelog='amrvac'

! Defaults for discretization methods
! extreme needs an extra cell in FV and FD!!
PPM_extrema=.false.
PPM_flatcd =.false.
PPM_flatsh =.false.
PPM_flatten = .false.
typeaverage='default'
tvdlfeps=one
BnormLF=.true.
nxdiffusehllc=0
flathllc=.false.
typeaxial='slab'
typecoord='covariant'
typespherical=1
slowsteps=-1
courantpar=0.8d0
typecourant='minimum'
dimsplit=.false.
typedimsplit='default'
typelimited='predictor'
typeinversion='KastaunC2P'
typeinversionresis='4DXIUNR'
use_backup_c2p = .false.
mcbeta=1.4d0
useprimitive=.true.
typetvd='roe'
typeemf='none'
old_bhac_safety = .false.
usr_W_limit_scheme = .false.
fixp_after_reconstruct = .false.
sourceimpl=.false.
sourceimplcycle=.false.
conduction=.false.
TCsaturate=.false.
TCphi=1.d0
energyonly=.false.
ncyclemax=1000
sourceparasts=.false.
parastsnu=0.001d0
ssplitdivb=.false.
{^IFMHDPHYS
ssplitdivb=.true.
}
ssplitresis=.false.
ssplituser=.false.
clean_init_divB='vecpot'
typeadvance='twostep'
evolve_hydro = .true.
do level=1,nlevelshi
   typefull1(level)='tvdlf'
   typepred1(level)='default'
   typelow1(level)='default'
   typelimiter1(level)='minmod'
   typegradlimiter1(level)='minmod'
enddo
use_del_zanna = .false.
use_4th_order_fd = .false.
use_6th_order_fd = .false.
face_to_center_order = 2
flatcd=.false.
flatsh=.false.
typesourcesplit='sfs'
do iw=1,nw
   typeentropy(iw)='nul'      ! Entropy fix type
end do
dtdiffpar=0.5d0
dtTCpar=0.5d0
dtpar=-1.d0
time_accurate=.true.
{#IFDEF MAGNETOFRICTION
cmf_c=0.3d0
cmf_y=0.2d0
cmf_divb=0.1d0
}

healpix_det_open = .false.
healpix_det_zsymm = .true.
healpix_det_num = 1
healpix_var_num = 0
healpix_nside = 16

!> --- m1 defaults
m1_closure_type = "Minerbo"
fileWeakhub = "/mnt/rafast/mcassing/Weakhub/dURCA_convenforbetaproc_DD2_grey_npe_20240723_Weakhub_grey.h5"
!for wavespeeds:
m1_actual_speeds = .true.
m1_radice_speeds = .false.
!for asymptotic flux correction
m1_parabolic = .false.
m1_frad_A2_A3 = .false.
m1_frad_A2_A = .true.
m1_erad_LLF = .false.
m1_frad_LLF = .false.
M1_FLUID_BACKREACT = .false.
m1_use_neutrinos = .true.
m1_use_muons = .false.
m1_use_photons = .false.
TESTqtC = 1.0d+10
m1_2_eas_updates = .true.
m1_E_atmo = 1.0d-15
m1_N_atmo = 1.0d-15
 !m1_i_nue = 1
 !m1_i_nuebar = 2
 !m1_i_nux = 3
 !m1_i_mu = 4
 !m1_i_mubar = 5
 !m1_i_photon = 6
m1_tset = 1.0d0 !0.4d0
m1_tset_backreact = 2.0d0
m1_rho_floor = 1.0d-12

! problem setup defaults
dxlone^D=zero;
nxlone^D=0;
iprob=1

! end defaults

! MPI reads from a file
open(unitpar,file=inifile,status='old')

! Start reading from standard input
primnames='default'

read(unitpar,filelist)
if (saveprim) call mpistop('saveprim is true, the code now separated into prim and cons, &
                                      why you need to saveprim to call one more c2p?')

if(TRIM(primnames)=='default'.and.mype==0) write(uniterr,*) &
   'Warning in ReadParameters: primnames not given!'

if(firstprocess .and. snapshotini<0) &
  call mpistop("Please restart from a snapshot when firstprocess=T")

read(unitpar,savelist)

if (mype.eq.0) then
   do ifile=1,nfile
      if(dtsave(ifile)<bigdouble/2) &
           write(unitterm,'(" DTSAVE  for file",i2," =",g12.5)') &
           ifile,dtsave(ifile)
      if(ditsave(ifile)<biginteger) &
           write(unitterm,'(" DITSAVE for file",i2," =",i10)') &
           ifile,ditsave(ifile)
      if(dtsave_rt(ifile)<bigdouble/2) &
           write(unitterm,'(" DTSAVE_RT  for file",i2," =",g12.5)') &
           ifile,dtsave_rt(ifile)
      if(tsave(1,ifile)==bigdouble.and.itsave(1,ifile)==biginteger.and. &
           dtsave(ifile)==bigdouble.and.ditsave(ifile)==biginteger.and.mype==0) &
           write(uniterr,*)'Warning in ReadParameters: ', &
           'No save condition for file ',ifile
   enddo
end if

! Consistency checks for the slices:
if (index(slice_type,'dat')>=1) sliceascii = .false.

do islice=1,nslices
   if(slicedir(islice) > ndim) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice,' direction',slicedir(islice),'larger than ndim=',ndim
   if(slicedir(islice) < 1) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice,' direction',slicedir(islice),'too small, should be [',1,ndim,']'
end do



read(unitpar,stoplist)
if (mype.eq.0) then
   if(itmax<biginteger)       write(unitterm,*) 'ITMAX=',itmax
   if(tmax<bigdouble)         write(unitterm,*) 'TMAX=',tmax
   if(dtmin>smalldouble)      write(unitterm,*) 'DTMIN=',dtmin
end if

if(itmax==biginteger .and. tmax==bigdouble.and.mype==0) write(uniterr,*) &
   'Warning in ReadParameters: itmax or tmax not given!'

if(residmin>=zero) then
   if (mype==0) write(unitterm,*)"SS computation with input value residmin"
   if(residmin<=smalldouble) call mpistop("Provide value for residual above smalldouble")
end if

wnames='default'

read(unitpar,methodlist)
if (B0field) call mpistop('Why you activated B0field?')
if (reconstruct_eps_for_T) call mpistop('reconstruct_eps_for_T is not ready yet')
if (typeinversion .ne. 'KastaunC2P') then
   call mpistop('Typeinversion: BHAC only has KastaunC2P primitive recovery, &
                 please select KastaunC2P')
endif
!if (typeinversionresis .ne. 'KastaunC2P') then
!   call mpistop('Typeinversionresis: BHAC only has KastaunC2P primitive recovery, &
!                 please select KastaunC2P')
!endif

if(TRIM(wnames)=='default') call mpistop("Provide wnames and restart code")
wnameslog=wnames

do level=1,nlevelshi
   !if(typefull1(level)=='tvdlf1'.and.typeadvance=='twostep') &
   !   call mpistop(" tvdlf1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hll1'.and.typeadvance=='twostep') &
   !   call mpistop(" hll1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hllc1'.and.typeadvance=='twostep') &
   !   call mpistop(" hllc1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hllcd1'.and.typeadvance=='twostep') &
   !   call mpistop(" hllcd1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='tvdmu1'.and.typeadvance=='twostep') &
   !   call mpistop(" tvdmu1 is onestep method, reset typeadvance=onestep!")
   if(typefull1(level)=='tvd'.and.typeadvance=='twostep') &
      call mpistop(" tvd is onestep method, reset typeadvance=onestep!")
   if(typefull1(level)=='tvd1'.and.typeadvance=='twostep') &
      call mpistop(" tvd1 is onestep method, reset typeadvance=onestep!")
   if(typefull1(level)=='tvd'.or.typefull1(level)=='tvd1')then 
      if(mype==0.and.(.not.dimsplit)) write(unitterm,*) &
         'Warning: setting dimsplit=T for tvd, as used for level=',level
      dimsplit=.true.
   endif

   if (typepred1(level)=='default') then
      select case (typefull1(level))
      case ('fd_tvdlf')
         typepred1(level)='fd_tvdlf'
      case ('fd_hll')
         typepred1(level)='fd_hll'
      case ('tvdlf')
         typepred1(level)='tvdlf'
      case ('tvdlfpos')
         typepred1(level)='tvdlfpos'
      case ('hll')
         typepred1(level)='hll'
      case ('tvdlf1','tvdmu1','tvd1','tvd','hll1','hllc1', &
            'hlld1','hllcd1','hlldd1','nul','source')
         typepred1(level)='nul'
      case default
         call mpistop("No default predictor for this full step")
      end select
   end if
end do

if (use_6th_order_fd) call mpistop('6th order fd is not working yet')

if (use_del_zanna .and. use_4th_order_fd) call mpistop('del zanna 4th order FD is not working yet &
                                                please use 4th order without del-zanna flux')
if (.not. use_4th_order_fd .and. face_to_center_order == 4) &
        call mpistop('Only if using 4th_order_fd --> 4th order face_to_center for CT')

select case (typeadvance)
case ("onestep")
   nstep=1
case ("twostep")
   nstep=2
case ("ImEx12")
   nstep=2
case ("ImEx112")
   nstep=1
case ("ImEx122")
   nstep=2
case ("ImEx222")
   nstep=2
case ("threestep")
   nstep=3
case ("fourstep","rk4","jameson","ssprk43")
   nstep=4
case ("ssprk54")
   nstep=5
case default
   call mpistop("Unknown typeadvance")
end select


do level=1,nlevelshi
  if (typelow1(level)=='default') then
   select case (typefull1(level))
   case ('fd_tvdlf')
      typelow1(level)='tvdlf1'
      !if (typeaxial .ne. 'cartesian') call mpistop('Strongly suggest you to not use &
      !                                   finite differencing in Curvilinear coordinates')
   case ('fd_hll')
      typelow1(level)='hll1'
      !if (typeaxial .ne. 'cartesian') call mpistop('Strongly suggest you to not use &
      !                                   finite differencing in Curvilinear coordinates')
   case ('tvdlf','tvdlf1','tvdmu','tvdmu1','tvd1','tvd','tvdlfpos')
      typelow1(level)='tvdlf1'
   case ('hll','hll1')
      typelow1(level)='hll1'
   case ('nul')
      typelow1(level)='nul'
   case ('source')
      typelow1(level)='source'
   case default
      call mpistop("No default typelow for this full step")
   end select
  end if
enddo


!code test
!write(*,*) PPM_extrema,PPM_flatcd,PPM_flatsh
!stop

! Harmonize the parameters for dimensional splitting and source splitting
if(typedimsplit   =='default'.and.     dimsplit)   typedimsplit='xyyx'
if(typedimsplit   =='default'.and..not.dimsplit)   typedimsplit='unsplit'
dimsplit   = typedimsplit   /='unsplit'

if (index(typecoord,'covariant')>=1) then
   covariant=.true.
else
   covariant=.false.
end if

if (typeaxial=="slab".and..not.covariant) then
   slab=.true.
else
   slab=.false.
end if

if (typeaxial=='spherical') then
   if (dimsplit) then
      if(mype==0)print *,'Warning: spherical symmetry needs dimsplit=F, resetting'
      dimsplit=.false.
   end if
end if

if (typecoord=='default') then
   typecoord = typeaxial
end if

if (ndim==1) dimsplit=.false.
if (.not.dimsplit.and.ndim>1) then
   select case (typeadvance)
   case ("ssprk54","ssprk43","fourstep", "rk4", "threestep", "twostep", "ImEx12", "ImEx122", "ImEx222")
      ! Runge-Kutta needs predictor
      typelimited="predictor"
      if(mype==0)print *,'typelimited to predictor for RK'
   end select
end if

{#IFDEF GLM
if (ssplitdivb .eqv. .false.) call mpistop('GLM needs ssplitdivb = .true.')
if (typeemf .ne. 'none') call mpistop('GLM and EMF averaging incompatible. Set typeemf= "none"')
}
{#IFDEF STAGGERED
if (typeemf .eq. 'none') typeemf = 'average' ! Defaulting to average Balsara Spicer
}
{#IFNDEF GLM
{#IFNDEF STAGGERED
if (typeemf .ne. 'average') call mpistop('FCT needs to be run with typeemf="average" ')
}
}



if (B0field) then
   if(mype==0)print *,'B0+B1 split for MHD'
   if (.not.typephys=='mhd') call mpistop("B0+B1 split for MHD only")
end if

do level=1,nlevelshi
   type_limiter(level) = limiter_type(typelimiter1(level))
   type_gradient_limiter(level) = limiter_type(typegradlimiter1(level))
end do

if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatsh.and.typephys=='rho')) then
    call mpistop(" PPM with flatsh=.true. can not be used with typephys='rho'!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatsh.and.typephys=='hdadiab')) then
     call mpistop(" PPM with flatsh=.true. can not be used with typephys='hdadiab'!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatcd.and.typephys=='hdadiab')) then
     call mpistop(" PPM with flatcd=.true. can not be used with typephys='hdadiab'!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatsh.and..not.useprimitive)) then
     call mpistop(" PPM with flatsh=.true. needs useprimitive=T!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatcd.and..not.useprimitive)) then
     call mpistop(" PPM with flatcd=.true. needs useprimitive=T!")
end if


{#IFDEF STAGGERED
  ! only ppm needs reconstruct for cell center B-field
  !if (any(typelimiter1(1:nlevelshi)== 'ppm')) then
  !  var_reconstruct(b1_) = .true.
  !  var_reconstruct(b2_) = .true.
  !  var_reconstruct(b3_) = .true.
  !endif
}

if (mype == 0) then
    write(*,*) '--------------The indices of Variables------------------'

    print*,'d_            :', d_           , 'reconstruct?', var_reconstruct(d_)
    print*,'dye_          :', dye_         , 'reconstruct?', var_reconstruct(dye_)
    print*,'s1_           :', s1_          , 'reconstruct?', var_reconstruct(s1_)
    print*,'s2_           :', s2_          , 'reconstruct?', var_reconstruct(s2_)
    print*,'s3_           :', s3_          , 'reconstruct?', var_reconstruct(s3_)
    print*,'tau_          :', tau_         , 'reconstruct?', var_reconstruct(tau_)
    print*,'Ds_           :', Ds_          , 'reconstruct?', var_reconstruct(Ds_)
    print*,'Dtr1_         :', Dtr1_        , 'reconstruct?', var_reconstruct(Dtr1_)
    print*,'pp_           :', pp_          , 'reconstruct?', var_reconstruct(pp_)
    print*,'ye_           :', ye_          , 'reconstruct?', var_reconstruct(ye_)
    print*,'rho_          :', rho_         , 'reconstruct?', var_reconstruct(rho_)
    print*,'u1_           :', u1_          , 'reconstruct?', var_reconstruct(u1_)
    print*,'u2_           :', u2_          , 'reconstruct?', var_reconstruct(u2_)
    print*,'u3_           :', u3_          , 'reconstruct?', var_reconstruct(u3_)
    print*,'b1_           :', b1_          , 'reconstruct?', var_reconstruct(b1_)
    print*,'b2_           :', b2_          , 'reconstruct?', var_reconstruct(b2_)
    print*,'b3_           :', b3_          , 'reconstruct?', var_reconstruct(b3_)
    print*,'T_eps_        :', T_eps_       , 'reconstruct?', var_reconstruct(T_eps_)
    print*,'cs2_          :', cs2_         , 'reconstruct?', var_reconstruct(cs2_)
    print*,'tr1_          :', tr1_         , 'reconstruct?', var_reconstruct(tr1_)
    print*,'lfac_         :', lfac_        , 'reconstruct?', var_reconstruct(lfac_)
    print*,'xi_           :', xi_          , 'reconstruct?', var_reconstruct(xi_)
    print*,'alp_metric_   :', alp_metric_  , 'reconstruct?', var_reconstruct(alp_metric_)
    print*,'beta_metric1_ :', beta_metric1_, 'reconstruct?', var_reconstruct(beta_metric1_)
    print*,'beta_metric2_ :', beta_metric2_, 'reconstruct?', var_reconstruct(beta_metric2_)
    print*,'beta_metric3_ :', beta_metric3_, 'reconstruct?', var_reconstruct(beta_metric3_)
    print*,'Xvec_metric1_ :', Xvec_metric1_, 'reconstruct?', var_reconstruct(Xvec_metric1_)
    print*,'Xvec_metric2_ :', Xvec_metric2_, 'reconstruct?', var_reconstruct(Xvec_metric2_)
    print*,'Xvec_metric3_ :', Xvec_metric3_, 'reconstruct?', var_reconstruct(Xvec_metric3_)
    print*,'psi_metric_   :', psi_metric_  , 'reconstruct?', var_reconstruct(psi_metric_)
    print*,'U_br_         :', U_br_        , 'reconstruct?', var_reconstruct(U_br_)
    print*,'U_br1_        :', U_br1_       , 'reconstruct?', var_reconstruct(U_br1_)
    print*,'U_br2_        :', U_br2_       , 'reconstruct?', var_reconstruct(U_br2_)
    print*,'U_br3_        :', U_br3_       , 'reconstruct?', var_reconstruct(U_br3_)
    print*,'h_00_         :', h_00_        , 'reconstruct?', var_reconstruct(h_00_)
    print*,'delta_br_     :', delta_br_    , 'reconstruct?', var_reconstruct(delta_br_)
    {#IFDEF M1
    {^KSP&
    print*,'nrad',^KSP,'_   :',nrad^KSP_    , 'reconstruct?', var_reconstruct(nrad^KSP_)
    print*,'erad',^KSP,'_   :',erad^KSP_    , 'reconstruct?', var_reconstruct(erad^KSP_)
    print*,'frad',^KSP,'1_  :',frad^KSP1_    , 'reconstruct?', var_reconstruct(frad^KSP1_)
    print*,'frad',^KSP,'2_  :',frad^KSP2_    , 'reconstruct?', var_reconstruct(frad^KSP2_)
    print*,'frad',^KSP,'3_  :',frad^KSP3_    , 'reconstruct?', var_reconstruct(frad^KSP3_)
    \}
    }
endif



read(unitpar,boundlist)
{#IFDEF DY_SP
If (dixB /= 2 .and. dixB /= 4) call mpistop('Dynamical spacetime only supports &
                                    ghostcells dixB = 2 or 4')
}
do idim=1,ndim
   periodB(idim)=(any(typeB(:,2*idim-1:2*idim)=='periodic'))
   aperiodB(idim)=(any(typeB(:,2*idim-1:2*idim)=='aperiodic'))
   if (periodB(idim).or.aperiodB(idim)) then
      do iw=1,nw
         if (typeB(iw,2*idim-1) .ne. typeB(iw,2*idim)) &
              call mpistop("Wrong counterpart in periodic boundary")
         if (typeB(iw,2*idim-1) /= 'periodic' .and. typeB(iw,2*idim-1) /= 'aperiodic') &
              call mpistop("Each dimension should either have all &
              or no variables periodic, some can be aperiodic")
      end do
   end if
end do

if (any(typelimiter1(1:nlevelshi)=='ppm').and.(dixB<4)) then
    call mpistop(" PPM works only with dixB>=4 !")
end if

if (any(typelimiter1(1:nlevelshi)=='mp5') .and. (dixB<3)) then
 call mpistop("mp5 needs at at least 3 ghost cells! Set dixB=3 in boundlist.")
end if
fv_n_interp = min(fv_n_interp, dixB * 2)

read(unitpar,amrlist)

select case (typeaxial)
{^NOONED
case ("spherical")
   xprob^LIM^DE=xprob^LIM^DE*two*dpi;
\}
case ("cylindrical")
   {if (^D==^PHI) then
      xprob^LIM^D=xprob^LIM^D*two*dpi;
   end if\}
end select
{
if(nxlone^D>1 .and. mod(nxlone^D,2)==0)then
   dxlone^D=(xprobmax^D-xprobmin^D)/dble(nxlone^D)
   if (mype==0) then
      write(unitterm,*)'Using ',nxlone^D,' cells in dimension ',^D
      write(unitterm,*)'level one dx(',^D,')=',dxlone^D
   end if
end if
\}

if({dxlone^D*}<smalldouble)then
   write(unitterm,*)'Wrong value(s) for level one dx:',dxlone^D
   call mpistop("Reset nxlone or dxlone!")
endif
^D&dx(^D,1)=dxlone^D;
if(mxnest>nlevelshi.or.mxnest<1)then
   write(unitterm,*)'Error: mxnest',mxnest,'>nlevelshi ',nlevelshi
   call mpistop("Reset nlevelshi and recompile!")
endif

if (flags(nflag_)>nw) then
   write(unitterm,*)'Error: flags(nw+1)=',flags(nw+1),'>nw ',nw
   call mpistop("Reset flags(nw+1)!")
end if
if (flags(nflag_)==0) errorestimate=0
if (flags(nflag_)<0) then
   if (mype==0) then
      write(unitterm,*) "flags(",nflag_,") can not be negative"
      call mpistop("")
   end if
end if
select case (errorestimate)
case (0)
   if (mype==0) write(unitterm,*)"Error estimation is user defined"
case (1)
   if (mype==0) write(unitterm,*)"Error estimation is richardson procedure"
case (2)
   if (mype==0) write(unitterm,*)"Error estimation is relative error"
case (3)
   if (mype==0) write(unitterm,*)"Error estimation is Lohner's scheme"
case (4)
   if (mype==0) write(unitterm,*)"Error estimation is Lohner's original scheme"
case default
   call mpistop("Unknown error estimator, change errorestimate")
end select
if (B0field.and.errorestimate==1) then
   call mpistop("No Richardson procedure in combination with B0")
end if

if (tfixgrid<bigdouble/2.0d0) then
   if(mype==0)print*,'Warning, at time=',tfixgrid,'the grid will be fixed'
end if
if (itfixgrid<biginteger/2) then
   if(mype==0)print*,'Warning, at iteration=',itfixgrid,'the grid will be fixed'
end if
if (ditregrid>1) then
   if(mype==0)print*,'Note, Grid is reconstructed once every',ditregrid,'iterations'
end if

do islice=1,nslices
select case(slicedir(islice))
{case(^D)
   if(slicecoord(islice)<xprobmin^D.or.slicecoord(islice)>xprobmax^D) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice, ' coordinate',slicecoord(islice),'out of bounds for dimension ',slicedir(islice)
\}
end select
end do

! now, -g depends on block_nx and dixB which can set in parafile
{^D& ixGhi^D = block_nx^D + 2*dixB \}

!write(*,*) 'ixGhi1, ixGhi2'
!write(*,*) ixGhi1, ixGhi2
!write(*,*) 'ixGlo1, ixGlo2'
!write(*,*) ixGlo1, ixGlo2
!stop
!{^D& ixGshi^D = ixGhi^D \}   !this is default, will be different in mod_physicaldata
if ({^D& ixGhi^D-ixGlo^D+1-2*dixB .lt. 2*dixB|.or.}) then
   call mpistop("For AMR you will need a physical cell to ghost-cell ratio of =>2, &
                 as too small physical cell, extremely slower for the sim")
end if



if (mype == 0) then
   write(*,*) '================Printing out the Grid parameters================'
   write(unitterm,*) 'mxnest :', mxnest 
   {^D&
   write(unitterm,*) 'nxlone^D :', nxlone^D
   write(unitterm,*) 'xprobmin^D :', xprobmin^D 
   write(unitterm,*) 'xprobmax^D :', xprobmax^D 
   write(unitterm,*) 'block_nx^D :', block_nx^D 
   \}
   write(*,*) '==================End of the Grid parameters==================='
endif

read(unitpar,paramlist)

if (dtpar>zero) time_accurate=.true.

if(.not.time_accurate) then
  if(residmin<=smalldouble .or. residmax==bigdouble) then
   if (mype==0) write(unitterm,*)"Non time_accurate SS computation needs values residmin and residmax"
   call mpistop("Provide values for residual bounds in stoplist")
  end if
end if

! Warn when too few blocks at start of simulation 
! -----------No need-----------! load balance can do better
!if (mype.eq.0 .and. snapshotini.eq.-1 .and. {^D& floor(dble(nxlone^D)/(dble(ixGhi^D)-2.0d0*dble(dixB))) |*} .lt. npe) then
!code test
!   call mpistop('Need at least as many blocks on level 1 as cores to initialize!')
!end if


read(unitpar,cfclist)
cfc_tol(1) = cfc_tol1
cfc_tol(2) = cfc_tol2
cfc_tol(3) = cfc_tol3
cfc_tol_evolve(1) = cfc_tol_evolve1
cfc_tol_evolve(2) = cfc_tol_evolve2
cfc_tol_evolve(3) = cfc_tol_evolve3



{#IFDEF DY_SP
select case (cfc_coordinate)
case ('cartesian')
coordinate = cartesian
case ('spherical')
coordinate = spherical
case ('cylindrical')
coordinate = cylindrical
case default
stop 'you chosen CFC, but you need to select coordinate'
end select
}

if (cfc_evolve) then
   if (.not. useprimitiveRel) call mpistop('CFC must use useprimitiveRel')
   call cfc_solver_activate()
   if (mype == 0) then
    write(*,*) 'CFC EVOLVE = .TRUE.'
   endif
else
   if (mype == 0) then
    write(*,*) 'CFC EVOLVE = .FALSE.'
   endif
endif

if (initialize_metric .or. cfc_evolve) then
   use_multigrid = .true.
   !metric_vars_interpolation = .true.
else
   use_multigrid = .false.
   !metric_vars_interpolation = .false.
endif

if (.not. cfc_evolve .and. use_gw_br) call mpistop('Cannot use GW backreaction without evolution xCFC solver')
if (use_gw_br) then 
   call gw_br_activate()
   !call mpistop('not finish implementation of gw br')
   gw_br_tol(1) = gw_br_tol1
   gw_br_tol(2) = gw_br_tol2
   gw_br_tol(3) = gw_br_tol3
   gw_br_tol_evolve(1) = gw_br_tol_evolve1
   gw_br_tol_evolve(2) = gw_br_tol_evolve2
   gw_br_tol_evolve(3) = gw_br_tol_evolve3

   if (.not. gw_br_use_I3_ij .and. gw_br_use_wi) &
       call mpistop('Do not turn on gw_br_use_wi when not using I3_ij for gw_br')

endif
{#IFNDEF GW_BR
if (use_gw_br) call mpistop('Please define GW_BR to use gw_br')
}
{#IFDEF GW_BR
!if (.not. use_gw_br) call mpistop('If you do not use_gw_br, turn off define GW_BR')
}

read(unitpar,eoslist)

select case(eos_type_input)
case("tabulated")
  call eos_tabulated_activate()
case("idealgas")
  call eos_idealgas_activate()
case("polytrope")
  call eos_polytrope_activate()
case("hybrid")
  call eos_hybrid_activate()
case default
  call mpistop('pls specify eos type or you missed the eos type input in para file')
end select
{#IFNDEF TABEOS
  if (eos_type == tabulated) then
     call mpistop('using tabulated eos, please define TABEOS in defintion file first')
  endif
}
{#IFDEF TABEOS
  if (eos_type /= tabulated) then
     call mpistop('You defined TABEOS in definition, if you dont use tabeos, pls undefine it')
  endif
}

read(unitpar, outflowlist)

read(unitpar, m1list)

close(unitpar)

if (mype==0) then
   print*,'Reading from inifile: ', trim(inifile)
   print*,'snapshotini         : ', snapshotini
   print*,'slicenext           : ', slicenext
   print*,'collapsenext        : ', collapsenext
   print*,'shellnext           : ', shellnext
   print*,'Filenameini         : ', trim(filenameini)
   print*,'Filenameout         : ', trim(filenameout)
   print*,'Converting?         : ', convert
   print*,'                                                                '
endif

if (mype.eq.0) print*,'-----------------------------------------------------------------------------'

if(mype==0) write(*,*) 'typelimiter1 at highest amr level', typelimiter1(mxnest) 
end subroutine readparameters
!=============================================================================
subroutine saveamrfile(ifile)

! following specific for Intel compiler and use on VIC3 with MPT
!DEC$ ATTRIBUTES NOINLINE :: write_snapshot

  use mod_slice, only: write_slice
  use mod_collapse, only: write_collapsed
  {^IFTHREED
    use mod_healpix_det, only: save_outflow_snap_to_hdf5, calculate_all_outflow
  }
  {#IFNDEF D1
  use mod_shell, only: write_shell }
include 'amrvacdef.f'
integer:: ifile
!-----------------------------------------------------------------------------
select case (ifile)
case (fileout_)
{#IFDEF HDF5
   call write_snapshot_hdf5
\}
{#IFNDEF HDF5
   if(endian_swap) typeparIO=-1
   select case (typeparIO)
      case (1) ! Parallel read and write
         call write_snapshot
      case (0,-2) ! Parallel read, serial write (-2)
         call write_snapshot_nopar
      case (-1) ! Serial read and write
         call write_snapshot_noparf
   end select
\}
{#IFDEF BOUNDARYDRIVER
   call write_boundary
}
!opedit: now we can also convert directly and will when autoconvert is set in inifile: 
   if (autoconvert) call generate_plotfile
{#IFDEF PARTICLES
   call write_particles_snapshot
}
case (fileslice_)
   call write_slice
case (filecollapse_)
   call write_collapsed
case (fileshell_)
{#IFNDEF D1
    call write_shell
}
{#IFDEF D1
   write(*,*) 'Shell output is not defined for 1D problems'  
}
case (filelog_)
   select case (typefilelog)
   case ('default')
      call printlog_default
   case ('special')
      call printlog_special
   case default
      call mpistop("Error in SaveFile: Unknown typefilelog")
   end select
case (fileanalysis_)
  call write_analysis
case (fileofdet_)
  {#IFDEF D3
    if (healpix_det_open) then
        call calculate_all_outflow(it)
        call save_outflow_snap_to_hdf5(it)
    endif
  }
case default
   write(*,*) 'No save method is defined for ifile=',ifile
   call mpistop("")
end select

! opedit: Flush stdout and stderr from time to time.
call flush(unitterm)

end subroutine saveamrfile
!=============================================================================
subroutine write_snapshot
use mod_forest
include 'amrvacdef.f'

integer :: file_handle, amode, igrid, Morton_no, iwrite
integer :: nx^D, count_forest
integer(kind=MPI_OFFSET_KIND) :: offset, init_fsize
integer, dimension(ngridshi) :: iorequest {#IFDEF STAGGERED , iorequest_stg}
integer, dimension(MPI_STATUS_SIZE,ngridshi) :: iostatus {#IFDEF STAGGERED , iostatus_stg}
integer, dimension(MPI_STATUS_SIZE) :: status
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------
if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"

init_fsize = 0
amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
! create if not exist
call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)
! truncate the file to be size zero
call MPI_FILE_SET_SIZE(file_handle, init_fsize, ierrmpi)

iorequest=MPI_REQUEST_NULL
{#IFDEF STAGGERED
iorequest_stg=MPI_REQUEST_NULL }

iwrite=0
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      call set_tmpGlobals(igrid)
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
   endif
   iwrite=iwrite+1
   offset=int(size_block_io {#IFDEF STAGGERED + size_block_stg_io} & 
        ,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IWRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
        iorequest(iwrite),ierrmpi)
{#IFDEF STAGGERED   
   offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IWRITE_AT(file_handle,offset,pws(igrid)%w,1,type_block_staggered_io, &
        iorequest_stg(iwrite),ierrmpi)
}
end do

if (iwrite>0) call MPI_WAITALL(iwrite,iorequest,iostatus,ierrmpi)
{#IFDEF STAGGERED
if (iwrite>0) call MPI_WAITALL(iwrite,iorequest_stg,iostatus_stg,ierrmpi)
}

call MPI_FILE_CLOSE(file_handle,ierrmpi)
if (mype==0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, &
                      file_handle,ierrmpi)

   call write_forest(file_handle, count_forest)

   {nx^D=ixMhi^D-ixMlo^D+1
   call MPI_FILE_WRITE(file_handle,nx^D,1,MPI_INTEGER,status,ierrmpi)\}
   call MPI_FILE_WRITE(file_handle,eqpar,neqpar+nspecialpar, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,count_forest,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,levmax,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndim,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nw,1,MPI_INTEGER,status,ierrmpi)
!{#IFDEF STAGGERED
   call MPI_FILE_WRITE(file_handle,nws,1,MPI_INTEGER,status,ierrmpi)
!}
   call MPI_FILE_WRITE(file_handle,neqpar+nspecialpar,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

   call MPI_FILE_CLOSE(file_handle,ierrmpi)
end if
snapshot=snapshot+1

end subroutine write_snapshot
!=============================================================================
subroutine write_snapshot_nopar
use mod_forest
include 'amrvacdef.f'


integer :: file_handle, amode, igrid, Morton_no, iwrite
integer :: nx^D, count_forest

integer(kind=MPI_OFFSET_KIND) :: offset, init_fsize

integer, allocatable :: iorecvstatus(:,:),ioastatus(:,:)
integer, allocatable :: igrecvstatus(:,:)
integer, allocatable :: igrid_recv(:) 
{#IFDEF STAGGERED
integer, allocatable :: iorecvstatus_stg(:,:),ioastatus_stg(:,:)
}
integer, dimension(MPI_STATUS_SIZE) :: status

integer  :: ipe,insend,inrecv,nrecv,nwrite
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------

call MPI_BARRIER(icomm,ierrmpi)

if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

nrecv=0
inrecv=0
nwrite=0
insend=0
iwrite=0

if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    itag=Morton_no
{#IFDEF STAGGERED
    itag_stg=Morton_no+ngridshi
}
    insend=insend+1
    if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      mygeo =>pgeo(igrid)
      {#IFNDEF DY_SP
      if(covariant)myM => mygeo%m
      }
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
    endif
    call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
    call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
    {#IFDEF STAGGERED
    call MPI_SEND(pws(igrid)%w,1,type_block_staggered_io, 0,itag_stg,icomm,ierrmpi)
    }
 end do
else 
 ! mype==0
 nwrite=(Morton_stop(0)-Morton_start(0)+1)

 ! master processor writes out
 write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"

 init_fsize = 0
 amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
 ! create if not exist
 call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)
 ! truncate the file to be size zero
 call MPI_FILE_SET_SIZE(file_handle, init_fsize, ierrmpi)


 ! writing his local data first
 do Morton_no=Morton_start(0),Morton_stop(0)
   igrid=sfc_to_igrid(Morton_no)
   iwrite=iwrite+1
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine,
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      mygeo =>pgeo(igrid)
      {#IFNDEF DY_SP
      if(covariant)myM => mygeo%m
      }
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
   endif

   offset=int(size_block_io {#IFDEF STAGGERED + size_block_stg_io} & 
         ,kind=MPI_OFFSET_KIND) &
           *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
                           MPI_STATUSES_IGNORE,ierrmpi)
{#IFDEF STAGGERED   
   offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
   call MPI_FILE_WRITE_AT(file_handle,offset,pws(igrid)%w,1,type_block_staggered_io, &
        MPI_STATUSES_IGNORE,ierrmpi)
   }
   

 end do
 ! write data communicated from other processors
 if(npe>1)then
  nrecv=(Morton_stop(npe-1)-Morton_start(1)+1)
  inrecv=0
  allocate(igrid_recv(nrecv))
  allocate(igrecvstatus(MPI_STATUS_SIZE,nrecv),iorecvstatus(MPI_STATUS_SIZE,nrecv))
  allocate(ioastatus(MPI_STATUS_SIZE,nrecv))
  {#IFDEF STAGGERED
  allocate(iorecvstatus_stg(MPI_STATUS_SIZE,nrecv))
  allocate(ioastatus_stg(MPI_STATUS_SIZE,nrecv))
  }

  do ipe =1, npe-1
   do Morton_no=Morton_start(ipe),Morton_stop(ipe)
     iwrite=iwrite+1
     itag=Morton_no
{#IFDEF STAGGERED
     itag_stg=Morton_no+ngridshi
}
     inrecv=inrecv+1
     call MPI_RECV(igrid_recv(inrecv),1,MPI_INTEGER, ipe,itag,icomm,&
                   igrecvstatus(:,inrecv),ierrmpi)

     allocate(pwio(igrid_recv(inrecv))%w(ixG^T,1:nw))
     call MPI_RECV(pwio(igrid_recv(inrecv))%w,1,type_block_io,ipe,itag,icomm,&
                   iorecvstatus(:,inrecv),ierrmpi)
{#IFDEF STAGGERED
     allocate(pwsio(igrid_recv(inrecv))%w(ixGlo^D-1:ixGhi^D,1:ndim))
     call MPI_RECV(pwsio(igrid_recv(inrecv))%w,1,type_block_staggered_io,ipe,itag_stg,icomm,&
                   iorecvstatus_stg(:,inrecv),ierrmpi)
}
     offset=int(size_block_io {#IFDEF STAGGERED + size_block_stg_io} & 
           ,kind=MPI_OFFSET_KIND) &
             *int(Morton_no-1,kind=MPI_OFFSET_KIND)
     call MPI_FILE_WRITE_AT(file_handle,offset,pwio(igrid_recv(inrecv))%w,1,&
                          type_block_io,ioastatus(:,inrecv),ierrmpi)
{#IFDEF STAGGERED
     offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
     call MPI_FILE_WRITE_AT(file_handle,offset,pwsio(igrid_recv(inrecv))%w,1,&
                          type_block_staggered_io,ioastatus_stg(:,inrecv),ierrmpi)
}
     deallocate(pwio(igrid_recv(inrecv))%w)
{#IFDEF STAGGERED
     deallocate(pwsio(igrid_recv(inrecv))%w)
}
   end do
  end do
  deallocate(igrecvstatus,iorecvstatus,ioastatus,igrid_recv)
  {#IFDEF STAGGERED
  deallocate(iorecvstatus_stg,ioastatus_stg)
  }
 end if
end if


call MPI_FILE_CLOSE(file_handle,ierrmpi)


if (mype==0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, &
                      file_handle,ierrmpi)

   call write_forest(file_handle, count_forest)
   {nx^D=ixMhi^D-ixMlo^D+1
   call MPI_FILE_WRITE(file_handle,nx^D,1,MPI_INTEGER,status,ierrmpi)\}
   call MPI_FILE_WRITE(file_handle,eqpar,neqpar+nspecialpar, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,count_forest,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,levmax,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndim,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nw,1,MPI_INTEGER,status,ierrmpi)
!{#IFDEF STAGGERED
   call MPI_FILE_WRITE(file_handle,nws,1,MPI_INTEGER,status,ierrmpi)
!}
   call MPI_FILE_WRITE(file_handle,neqpar+nspecialpar,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

   call MPI_FILE_CLOSE(file_handle,ierrmpi)
end if

snapshot=snapshot+1

call MPI_BARRIER(icomm,ierrmpi)

end subroutine write_snapshot_nopar
!=============================================================================
subroutine write_snapshot_noparf
use mod_forest
include 'amrvacdef.f'

integer :: igrid, Morton_no
integer :: nx^D

integer, allocatable :: iostatus(:,:),iorecvstatus(:,:),ioastatus(:,:)
integer, allocatable :: igrecvstatus(:,:)
integer, allocatable :: iorequest(:),igrid_recv(:) 
{#IFDEF STAGGERED
integer, allocatable :: iostatus_stg(:,:),iorecvstatus_stg(:,:),ioastatus_stg(:,:)
integer, allocatable :: iorequest_stg(:)
}
integer  :: ipe,insend,inrecv,nrecv,nwrite,count_forest
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------
call MPI_BARRIER(icomm,ierrmpi)

if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

nrecv=0
inrecv=0
nwrite=0
insend=0

if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    itag=Morton_no
{#IFDEF STAGGERED
    itag_stg=Morton_no+ngridshi
}
    insend=insend+1
    if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      mygeo =>pgeo(igrid)
      {#IFNDEF DY_SP
      if(covariant)myM => mygeo%m
      }
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
    endif
    call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
    call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
    {#IFDEF STAGGERED
    call MPI_SEND(pws(igrid)%w,1,type_block_staggered_io, 0,itag_stg,icomm,ierrmpi)
    }
 end do
else 
 ! mype==0
 nwrite=(Morton_stop(0)-Morton_start(0)+1)
 allocate(iorequest(nwrite),iostatus(MPI_STATUS_SIZE,nwrite))
{#IFDEF STAGGERED
 allocate(iorequest_stg(nwrite),iostatus_stg(MPI_STATUS_SIZE,nwrite))
}
 iorequest=MPI_REQUEST_NULL
{#IFDEF STAGGERED
 iorequest_stg=MPI_REQUEST_NULL}

 ! master processor writes out
 write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"
 if(endian_swap) then
  {#IFNDEF BIGENDIAN
   open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
        status='replace',convert='BIG_ENDIAN')
  }
  {#IFDEF BIGENDIAN
   open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
        status='replace',convert='LITTLE_ENDIAN')
  }
 else
   open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
        status='replace')
 end if
 ! writing his local data first
 do Morton_no=Morton_start(0),Morton_stop(0)
   igrid=sfc_to_igrid(Morton_no)
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine,
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      mygeo =>pgeo(igrid)
      {#IFNDEF DY_SP
      if(covariant)myM => mygeo%m
      }
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
   endif
   write(unitsnapshot) pw(igrid)%w(ixM^T,1:nw)
{#IFDEF STAGGERED
   write(unitsnapshot) pws(igrid)%w(ixMlo^D-1:ixMhi^D,1:nws)
}
 end do
 ! write data communicated from other processors
 if(npe>1)then
  nrecv=(Morton_stop(npe-1)-Morton_start(1)+1)
  inrecv=0
  allocate(igrid_recv(nrecv))
  allocate(igrecvstatus(MPI_STATUS_SIZE,nrecv),iorecvstatus(MPI_STATUS_SIZE,nrecv))
  allocate(ioastatus(MPI_STATUS_SIZE,nrecv))
  {#IFDEF STAGGERED
  allocate(iorecvstatus_stg(MPI_STATUS_SIZE,nrecv))
  allocate(ioastatus_stg(MPI_STATUS_SIZE,nrecv))
  }

  do ipe =1, npe-1
   do Morton_no=Morton_start(ipe),Morton_stop(ipe)
     itag=Morton_no
{#IFDEF STAGGERED
     itag_stg=Morton_no+ngridshi
}
     inrecv=inrecv+1
     call MPI_RECV(igrid_recv(inrecv),1,MPI_INTEGER, ipe,itag,icomm,&
                   igrecvstatus(:,inrecv),ierrmpi)

     allocate(pwio(igrid_recv(inrecv))%w(ixG^T,1:nw))
     call MPI_RECV(pwio(igrid_recv(inrecv))%w,1,type_block_io,ipe,itag,icomm,&
                   iorecvstatus(:,inrecv),ierrmpi)
{#IFDEF STAGGERED
     allocate(pwsio(igrid_recv(inrecv))%w(ixGlo^D-1:ixGhi^D,1:ndim))
     call MPI_RECV(pwsio(igrid_recv(inrecv))%w,1,type_block_staggered_io,ipe,itag_stg,icomm,&
                   iorecvstatus_stg(:,inrecv),ierrmpi)
}
     write(unitsnapshot) pwio(igrid_recv(inrecv))%w(ixM^T,1:nw)
{#IFDEF STAGGERED
     write(unitsnapshot) pwsio(igrid_recv(inrecv))%w(ixMlo^D-1:ixMhi^D,1:nws)
}
     deallocate(pwio(igrid_recv(inrecv))%w)
{#IFDEF STAGGERED
     deallocate(pwsio(igrid_recv(inrecv))%w)
}
   end do
  end do
  deallocate(igrecvstatus,iorecvstatus,ioastatus,igrid_recv)
  {#IFDEF STAGGERED
  deallocate(iorecvstatus_stg,ioastatus_stg)
  }
 end if
end if

if(nwrite>0)then
  call MPI_WAITALL(nwrite,iorequest,iostatus,ierrmpi) 
{#IFDEF STAGGERED
  call MPI_WAITALL(nwrite,iorequest_stg,iostatus_stg,ierrmpi)
}
  if(mype==0) then
    deallocate(iorequest,iostatus)
    {#IFDEF STAGGERED
    deallocate(iorequest_stg, iostatus_stg)}
  end if
end if

if(mype==0) then
  call write_forest(unitsnapshot, count_forest)
  {nx^D=ixMhi^D-ixMlo^D+1
  write(unitsnapshot) nx^D\}
  write(unitsnapshot) eqpar
  write(unitsnapshot) count_forest
  write(unitsnapshot) nleafs
  write(unitsnapshot) levmax 
  write(unitsnapshot) ndim 
  write(unitsnapshot) ndir
  write(unitsnapshot) nw
  write(unitsnapshot) nws !! Staggered variables
  write(unitsnapshot) neqpar+nspecialpar
  write(unitsnapshot) it
  write(unitsnapshot) t
  close(unitsnapshot)
end if

snapshot=snapshot+1

call MPI_BARRIER(icomm,ierrmpi)

end subroutine write_snapshot_noparf
!=============================================================================
subroutine read_snapshot
use mod_forest
include 'amrvacdef.f'

integer :: file_handle, amode, igrid, Morton_no, iread, count_forest
integer :: levmaxini, ndimini, ndirini, nwini, nwsini, neqparini, nxini^D
integer :: global_integer(1:9)

integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(ngridshi) :: iorequest {#IFDEF STAGGERED , iorequest_stg}
integer, dimension(MPI_STATUS_SIZE,ngridshi) :: iostatus {#IFDEF STAGGERED , iostatus_stg}
integer, dimension(MPI_STATUS_SIZE) :: status
character(len=80) :: filename
logical :: fexist
!-----------------------------------------------------------------------------
!!!call MPI_BARRIER(icomm,ierrmpi)
! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"

if(mype==0) then
  inquire(file=filename,exist=fexist)
  if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found!")
endif

amode=MPI_MODE_RDONLY
call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)


!call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
!call MPI_TYPE_EXTENT(MPI_INTEGER,size_int,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size_int,ierrmpi)

offset=-int( 9*size_int+size_double,kind=MPI_OFFSET_KIND)
!call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

call MPI_FILE_READ_ALL(file_handle,global_integer,9,MPI_INTEGER,status,ierrmpi)
count_forest = global_integer(1)
nleafs = global_integer(2)
nleafs_active = nleafs
levmaxini = global_integer(3)
ndimini = global_integer(4)
ndirini = global_integer(5)
nwini = global_integer(6)
nwsini = global_integer(7)
neqparini = global_integer(8)
it = global_integer(9)
call MPI_FILE_READ_ALL(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

! check if settings are suitable for restart
if (levmaxini>mxnest) then
   if (mype==0) write(*,*) "number of levels in restart file = ",levmaxini
   if (mype==0) write(*,*) "mxnest = ",mxnest
   call mpistop("mxnest should be at least number of levels in restart file")
end if
if (ndimini/=ndim) then
   if (mype==0) write(*,*) "ndim in restart file = ",ndimini
   if (mype==0) write(*,*) "ndim = ",ndim
   call mpistop("reset ndim to ndim in restart file")
end if
if (ndirini/=ndir) then
   if (mype==0) write(*,*) "ndir in restart file = ",ndirini
   if (mype==0) write(*,*) "ndir = ",ndir
   call mpistop("reset ndir to ndir in restart file")
end if
if (nw/=nwini) then
   if (mype==0) write(*,*) "nw=",nw," and nw in restart file=",nwini
   call mpistop("currently, changing nw at restart is not allowed")
end if
{#IFDEF STAGGERED
if (nws/=nwsini) then
   if (mype==0) write(*,*) "nws=",nws," and nws in restart file=",nwsini
   call mpistop("currently, changing nws at restart is not allowed")
end if
}
offset=offset-int(ndimini*size_int+neqparini*size_double,kind=MPI_OFFSET_KIND)
!call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

{call MPI_FILE_READ_ALL(file_handle,nxini^D,1,MPI_INTEGER,status,ierrmpi)\}
if (ixGhi^D/=nxini^D+2*dixB|.or.) then
   if (mype==0) write(*,*) "Error: reset resolution to ",nxini^D+2*dixB
   call mpistop("change with setamrvac")
end if
neqparini=min(neqparini,neqpar+nspecialpar)
call MPI_FILE_READ_ALL(file_handle,eqpar,neqparini, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)


call read_forest(file_handle, count_forest)

iorequest=MPI_REQUEST_NULL
{#IFDEF STAGGERED
iorequest_stg=MPI_REQUEST_NULL }

iread=0
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   call alloc_node(igrid)
   iread=iread+1
   offset=int(size_block_io {#IFDEF STAGGERED + size_block_stg_io} &
        ,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IREAD_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
                          iorequest(iread),ierrmpi)
{#IFDEF STAGGERED   
   offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IREAD_AT(file_handle,offset,pws(igrid)%w,1,type_block_staggered_io, &
        iorequest_stg(iread),ierrmpi)
}
end do

if (iread>0) call MPI_WAITALL(iread,iorequest,iostatus,ierrmpi)
{#IFDEF STAGGERED
if (iread>0) call MPI_WAITALL(iread,iorequest_stg,iostatus_stg,ierrmpi)
}

call MPI_FILE_CLOSE(file_handle,ierrmpi)

!!!call MPI_BARRIER(icomm,ierrmpi)
end subroutine read_snapshot
!=============================================================================
subroutine read_snapshot_nopar
use mod_forest
include 'amrvacdef.f'

double precision :: wio(ixG^T,1:nw) {#IFDEF STAGGERED ,wsio(ixGlo^D-1:ixGhi^D,1:nws)}
integer :: file_handle, amode, igrid, Morton_no, iread
integer :: levmaxini, ndimini, ndirini, nwini, nwsini, neqparini, nxini^D, count_forest
integer :: global_integer(1:9)

integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(ngridshi) :: iorequest {#IFDEF STAGGERED , iorequest_stg}
integer, dimension(MPI_STATUS_SIZE) :: status
integer, dimension(MPI_STATUS_SIZE) :: iostatus {#IFDEF STAGGERED , iostatus_stg}

integer, allocatable :: iorecvstatus(:,:) {#IFDEF STAGGERED , iorecvstatus_stg(:,:)}
integer :: ipe,inrecv,nrecv
integer :: sendini(9+^ND)
character(len=80) :: filename
logical :: fexist
!-----------------------------------------------------------------------------
!!!call MPI_BARRIER(icomm,ierrmpi)
! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"
if (mype==0) then
 inquire(file=filename,exist=fexist)
 if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found!")
 amode=MPI_MODE_RDONLY
 call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)

 !call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
 !call MPI_TYPE_EXTENT(MPI_INTEGER,size_int,ierrmpi)
 call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
 call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size_int,ierrmpi)

 offset=-int( 9*size_int+size_double,kind=MPI_OFFSET_KIND)
 !call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
 call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

 call MPI_FILE_READ_ALL(file_handle,global_integer,9,MPI_INTEGER,status,ierrmpi)
 count_forest = global_integer(1)
 nleafs = global_integer(2)
 nleafs_active = nleafs
 levmaxini = global_integer(3)
 ndimini = global_integer(4)
 ndirini = global_integer(5)
 nwini = global_integer(6)
 nwsini = global_integer(7)
 neqparini = global_integer(8)
 it = global_integer(9)
 call MPI_FILE_READ(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)
 ! check if settings are suitable for restart
 if (levmaxini>mxnest) then
    write(*,*) "number of levels in restart file = ",levmaxini
    write(*,*) "mxnest = ",mxnest
    call mpistop("mxnest should be at least number of levels in restart file")
 end if
 if (ndimini/=ndim) then
    write(*,*) "ndim in restart file = ",ndimini
    write(*,*) "ndim = ",ndim
    call mpistop("reset ndim to ndim in restart file")
 end if
 if (ndirini/=ndir) then
    write(*,*) "ndir in restart file = ",ndirini
    write(*,*) "ndir = ",ndir
    call mpistop("reset ndir to ndir in restart file")
 end if
 if (nw/=nwini) then
    write(*,*) "nw=",nw," and nw in restart file=",nwini
    call mpistop("currently, changing nw at restart is not allowed")
 end if
{#IFDEF STAGGERED
 if (nws/=nwsini) then
   if (mype==0) write(*,*) "nws=",nws," and nws in restart file=",nwsini
   call mpistop("currently, changing nws at restart is not allowed")
 end if
}
 offset=offset-int(ndimini*size_int+neqparini*size_double,kind=MPI_OFFSET_KIND)
 !call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
 call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

 {call MPI_FILE_READ(file_handle,nxini^D,1,MPI_INTEGER,status,ierrmpi)\}
 if (ixGhi^D/=nxini^D+2*dixB|.or.) then
    write(*,*) "Error: reset resolution to ",nxini^D+2*dixB
    call mpistop("change with setamrvac")
 end if
 neqparini=min(neqparini,neqpar+nspecialpar)
 call MPI_FILE_READ(file_handle,eqpar,neqparini, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)
end if

! broadcast the global parameters first
if (npe>1) then
  if (mype==0) then
     sendini=(/count_forest,nleafs,levmaxini,ndimini,ndirini,nwini, &
               nwsini,neqparini,it,^D&nxini^D /)
  end if
  call MPI_BCAST(sendini,9+^ND,MPI_INTEGER,0,icomm,ierrmpi)
  count_forest=sendini(1);nleafs=sendini(2);levmaxini=sendini(3);
  ndimini=sendini(4);ndirini=sendini(5);nwini=sendini(6);
  nwsini=sendini(7);neqparini=sendini(8);it=sendini(9);
  nxini^D=sendini(9+^D);
  nleafs_active = nleafs
  call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
  call MPI_BCAST(eqpar,neqparini,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if

call read_forest(file_handle, count_forest)

if (mype==0)then
   iread=0
!!! HO: What is this for? (iorequest)
   iorequest=MPI_REQUEST_NULL
   {#IFDEF STAGGERED
   iorequest_stg=MPI_REQUEST_NULL}
   do Morton_no=Morton_start(0),Morton_stop(0)
      igrid=sfc_to_igrid(Morton_no)
      call alloc_node(igrid)
      iread=iread+1
      offset=int(size_block_io {#IFDEF STAGGERED + size_block_stg_io} &
            ,kind=MPI_OFFSET_KIND) &
             *int(Morton_no-1,kind=MPI_OFFSET_KIND)
      call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
                            iostatus,ierrmpi)
      {#IFDEF STAGGERED
      offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
      call MPI_FILE_READ_AT(file_handle,offset,pws(igrid)%w,1,type_block_staggered_io, &
                            iostatus_stg,ierrmpi)
       }
   end do
   if (npe>1) then
    do ipe=1,npe-1
     do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       iread=iread+1
       itag=Morton_no
       {#IFDEF STAGGERED
       itag_stg=Morton_no+ngridshi}
       offset=int(size_block_io {#IFDEF STAGGERED + size_block_stg_io} &
             ,kind=MPI_OFFSET_KIND) &
              *int(Morton_no-1,kind=MPI_OFFSET_KIND)
       call MPI_FILE_READ_AT(file_handle,offset,wio,1,type_block_io, &
                            iostatus,ierrmpi)
       call MPI_SEND(wio,1,type_block_io, ipe,itag,icomm,ierrmpi)
       {#IFDEF STAGGERED
       offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
       call MPI_FILE_READ_AT(file_handle,offset,wsio,1,type_block_staggered_io, &
                            iostatus_stg,ierrmpi)
       call MPI_SEND(wsio,1,type_block_staggered_io, ipe,itag_stg,icomm,ierrmpi)
       }
     end do
    end do
   end if
   call MPI_FILE_CLOSE(file_handle,ierrmpi)
else
   nrecv=(Morton_stop(mype)-Morton_start(mype)+1)
   allocate(iorecvstatus(MPI_STATUS_SIZE,nrecv))
{#IFDEF STAGGERED
   allocate(iorecvstatus_stg(MPI_STATUS_SIZE,nrecv))}
   inrecv=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid=sfc_to_igrid(Morton_no)
      itag=Morton_no
      {#IFDEF STAGGERED
      itag_stg=Morton_no+ngridshi}
      call alloc_node(igrid)
      inrecv=inrecv+1
      call MPI_RECV(pw(igrid)%w,1,type_block_io,0,itag,icomm,&
                    iorecvstatus(:,inrecv),ierrmpi)
      {#IFDEF STAGGERED
      call MPI_RECV(pws(igrid)%w,1,type_block_staggered_io,0,itag_stg,icomm,&
                    iorecvstatus_stg(:,inrecv),ierrmpi)
      }
   end do
   deallocate(iorecvstatus {#IFDEF STAGGERED , iorecvstatus_stg})
end if

call MPI_BARRIER(icomm,ierrmpi)

end subroutine read_snapshot_nopar
!=============================================================================
subroutine read_snapshot_noparf
use mod_forest
include 'amrvacdef.f'

double precision :: wio(ixG^T,1:nw) {#IFDEF STAGGERED ,wsio(ixGlo^D-1:ixGhi^D,1:nws)}
integer :: file_handle, amode, igrid, Morton_no, iread
integer :: levmaxini, ndimini, ndirini, nwini, nwsini, neqparini, nxini^D, count_forest

integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(ngridshi) :: iorequest {#IFDEF STAGGERED , iorequest_stg}
integer, dimension(MPI_STATUS_SIZE) :: status
integer, dimension(MPI_STATUS_SIZE) :: iostatus {#IFDEF STAGGERED , iostatus_stg}

integer, allocatable :: iorecvstatus(:,:) {#IFDEF STAGGERED , iorecvstatus_stg(:,:)}
integer :: ipe,inrecv,nrecv
integer :: sendini(9+^ND)
character(len=80) :: filename
logical :: fexist
integer :: file_size, idx
!-----------------------------------------------------------------------------
call MPI_BARRIER(icomm,ierrmpi)
! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"
if (mype==0) then
 inquire(file=filename,exist=fexist,size=file_size)
 if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found!")
 ! Open and read file
 size_int=sizeof(4)
 size_double=sizeof(1.0d0)

 offset=file_size+1-(9*size_int+size_double)
 open(unitsnapshot,file=filename,form='unformatted',access='stream',&
      status='old',action='read')

 read(unitsnapshot,pos=offset) count_forest
 read(unitsnapshot) nleafs
 read(unitsnapshot) levmaxini
 read(unitsnapshot) ndimini
 read(unitsnapshot) ndirini
 read(unitsnapshot) nwini
 read(unitsnapshot) nwsini
 read(unitsnapshot) neqparini

 read(unitsnapshot) it

 read(unitsnapshot) t

 ! check if settings are suitable for restart
 if (levmaxini>mxnest) then
    write(*,*) "number of levels in restart file = ",levmaxini
    write(*,*) "mxnest = ",mxnest
    call mpistop("mxnest should be at least number of levels in restart file")
 end if
 if (ndimini/=ndim) then
    write(*,*) "ndim in restart file = ",ndimini
    write(*,*) "ndim = ",ndim
    call mpistop("reset ndim to ndim in restart file")
 end if
 if (ndirini/=ndir) then
    write(*,*) "ndir in restart file = ",ndirini
    write(*,*) "ndir = ",ndir
    call mpistop("reset ndir to ndir in restart file")
 end if
 if (nw/=nwini) then
    write(*,*) "nw=",nw," and nw in restart file=",nwini
    call mpistop("currently, changing nw at restart is not allowed")
 end if
{#IFDEF STAGGERED
 if (nws/=nwsini) then
   if (mype==0) write(*,*) "nws=",nws," and nws in restart file=",nwsini
   call mpistop("currently, changing nws at restart is not allowed")
 end if
}

 ! ---
 offset=offset-(ndimini*size_int+neqparini*size_double)
 {read(unitsnapshot,pos=offset+(^D-1)*size_int)nxini^D\}
 if (ixGhi^D/=nxini^D+2*dixB|.or.) then
    write(*,*) "Error: reset resolution to ",nxini^D+2*dixB
    call mpistop("change with setamrvac")
 end if
 neqparini=min(neqparini,neqpar+nspecialpar)
 do idx=1,neqparini
     read(unitsnapshot) eqpar(idx)
 end do

 close(unitsnapshot) 
end if ! mype=0

! Broadcast global parameters
if (npe>1) then

  if (mype==0) then
     sendini=(/count_forest,nleafs,levmaxini,ndimini,ndirini,nwini, &
               nwsini,neqparini,it,^D&nxini^D /)
  end if
  call MPI_BCAST(sendini,9+^ND,MPI_INTEGER,0,icomm,ierrmpi)
  count_forest=sendini(1);nleafs=sendini(2);levmaxini=sendini(3);
  ndimini=sendini(4);ndirini=sendini(5);nwini=sendini(6);
  nwsini=sendini(7);neqparini=sendini(8);it=sendini(9);
  nxini^D=sendini(9+^D);
  nleafs_active = nleafs
  call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
  call MPI_BCAST(eqpar,neqparini,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if

open(unitsnapshot,file=filename,form='unformatted',access='stream',&
     status='old',action='read')

! HO: We pass an MPI file handle to avoid changing the interface,
! even though it is not used
call read_forest(file_handle, count_forest)

! Read and communicate arrays
if (mype==0)then
   iread=0

   do Morton_no=Morton_start(0),Morton_stop(0)
      igrid=sfc_to_igrid(Morton_no)
      call alloc_node(igrid)
      iread=iread+1
      offset=1+(size_block_io {#IFDEF STAGGERED + size_block_stg_io}) &
              *(Morton_no-1)

      read(unitsnapshot,pos=offset) pw(igrid)%w(ixM^T,:)
      {#IFDEF STAGGERED
      offset=offset + size_block_io
      read(unitsnapshot,pos=offset) pws(igrid)%w(ixMlo^D-1:ixMhi^D,:) }
   end do
   if (npe>1) then
    do ipe=1,npe-1
     do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       iread=iread+1
       itag=Morton_no
       {#IFDEF STAGGERED
       itag_stg=Morton_no+ngridshi}
       offset=1+(size_block_io {#IFDEF STAGGERED + size_block_stg_io}) &
              *(Morton_no-1)
       read(unitsnapshot,pos=offset) wio(ixM^T,:)
       call MPI_SEND(wio,1,type_block_io, ipe,itag,icomm,ierrmpi)
       {#IFDEF STAGGERED
       offset=offset + size_block_io
       read(unitsnapshot,pos=offset) wsio(ixMlo^D-1:ixMhi^D,:) 
       call MPI_SEND(wsio,1,type_block_staggered_io,ipe,itag_stg,icomm,ierrmpi)
       }
     end do
    end do
   end if

   close(unitsnapshot)

else
   ! mype /= 0
   nrecv=(Morton_stop(mype)-Morton_start(mype)+1)

   allocate(iorecvstatus(MPI_STATUS_SIZE,nrecv))
{#IFDEF STAGGERED
   allocate(iorecvstatus_stg(MPI_STATUS_SIZE,nrecv))}

   iorequest=MPI_REQUEST_NULL
   {#IFDEF STAGGERED
   iorequest_stg=MPI_REQUEST_NULL}

   inrecv=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid=sfc_to_igrid(Morton_no)
      itag=Morton_no
      {#IFDEF STAGGERED
      itag_stg=Morton_no+ngridshi}
      call alloc_node(igrid)
      inrecv=inrecv+1
      call MPI_RECV(pw(igrid)%w,1,type_block_io,0,itag,icomm,&
                    iorecvstatus(:,inrecv),ierrmpi)
      {#IFDEF STAGGERED
      call MPI_RECV(pws(igrid)%w,1,type_block_staggered_io,0,itag_stg,icomm,&
                    iorecvstatus_stg(:,inrecv),ierrmpi)
      }
   end do
   deallocate(iorecvstatus {#IFDEF STAGGERED , iorecvstatus_stg})
end if

if(nrecv>0)then
  call MPI_WAITALL(nrecv,iorequest,iostatus,ierrmpi) 
{#IFDEF STAGGERED
  call MPI_WAITALL(nrecv,iorequest_stg,iostatus_stg,ierrmpi)
}

end if

end subroutine read_snapshot_noparf
!=============================================================================
subroutine printlog_default

! printlog: calculates volume averaged mean values 
use mod_timing
use mod_forest,only:nleafs,nleafs_active,nleafs_level
include 'amrvacdef.f'

logical          :: fileopen
integer          :: iigrid, igrid, level, iw, i
double precision :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
double precision :: dvolume(ixG^T), volumeflat(1:nlevelshi)
integer          :: numlevels, nx^D, nc, ncells, dit
double precision :: dtTimeLast, now, cellupdatesPerSecond, activeBlocksPerCore, wctPerCodeTime, timeToFinish
integer, dimension(1:nlevelshi) :: isum_send, isum_recv
double precision, dimension(1:nw+1+nlevelshi) :: dsum_send, dsum_recv
character(len=80) :: filename
character(len=2048) :: line
logical, save :: opened=.false.
integer :: amode, status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------

volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
wmean(1:nw)= zero

do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   volumeflat(level)=volumeflat(level)+ &
          {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*}
   if (slab) then
      dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
   else
      dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
      volume(level)=volume(level)+sum(dvolume(ixM^T))
   end if
   do iw=1,nw
      wmean(iw)=wmean(iw)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,iw))
   end do
end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))

numlevels=levmax-levmin+1
dsum_send(1:nw)=wmean(1:nw)
dsum_send(nw+1)=voltotal
dsum_send(nw+2:nw+1+numlevels)=volumeflat(levmin:levmax)
call MPI_REDUCE(dsum_send,dsum_recv,nw+1+numlevels,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,icomm,ierrmpi)

if (mype==0) then

! To compute cell updates per second, we do the following:
nx^D=ixMhi^D-ixMlo^D+1;
nc={nx^D*}
ncells = nc * nleafs_active
! assumes the number of active leafs haven't changed since last compute.
now        = MPI_WTIME()
dit        = it - itTimeLast
dtTimeLast = now - timeLast
itTimeLast = it
timeLast   = now
cellupdatesPerSecond = dble(ncells) * dble(nstep) * dble(dit) / (dtTimeLast * dble(npe))
! blocks per core:
activeBlocksPerCore = dble(nleafs_active) / dble(npe)
! Wall clock time per code time unit in seconds:
if (dt .gt. zero) then 
   wctPerCodeTime = dtTimeLast / (dble(dit) * dt)
else
   wctPerCodeTime = zero
end if

! Wall clock time to finish in hours:
timeToFinish = (tmax - t) * wctPerCodeTime / 3600.0d0

   wmean(1:nw)=dsum_recv(1:nw)
   voltotal=dsum_recv(nw+1)
   volumeflat(levmin:levmax)=dsum_recv(nw+2:nw+1+numlevels)

   wmean=wmean/voltotal

   ! determine coverage in coordinate space
   volprob={(xprobmax^D-xprobmin^D)|*}
   volumeflat(levmin:levmax)=volumeflat(levmin:levmax)/volprob

   if (.not.opened) then
      ! generate filename
      write(filename,"(a,a)") TRIM(filenamelog),".log"

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      amode=ior(amode,MPI_MODE_APPEND)
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, &
                         MPI_INFO_NULL,log_fh,ierrmpi)
      opened=.true.

      if (snapshotini==-1) then
          call MPI_FILE_WRITE(log_fh,fileheadout,len_trim(fileheadout), &
                              MPI_CHARACTER,status,ierrmpi)
          !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
          call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)
      endif

      i=len_trim(wnameslog)-1
      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<1024) write(wnameslog(i:i+1),"(a,i1)") "c",level
          else
            if (i+2<1024) write(wnameslog(i:i+2),"(a,i2)") "c",level
          endif
      end do

      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<1024) write(wnameslog(i:i+1),"(a,i1)") "n",level
          else
            if (i+2<1024) write(wnameslog(i:i+2),"(a,i2)") "n",level
          endif
      end do
      if (time_accurate) then
         if(residmin>smalldouble) then
           write(line,'(a15,a1024)')"it   t  dt res ",wnameslog
         else
           write(line,'(a15,a1024)')"it   t   dt    ",wnameslog
         endif
      else
         if(residmin>smalldouble) then
           write(line,'(a7,a1024)')"it res ",wnameslog
         else
           write(line,'(a7,a1024)')"it     ",wnameslog
         endif
      end if

      line=trim(line)//"| Xload Xmemory 'Cell_Updates /second/core'"
      line=trim(line)//" 'Active_Blocks/Core' 'Wct Per Code Time [s]' 'TimeToFinish [hrs]'"

      
      if (snapshotini==-1) then
          call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, &
                              status,ierrmpi)
      endif
   end if
   !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   if (time_accurate) then
      if(residmin>smalldouble) then
         write(line,'(i7,3(es12.4))')it,t,dt,residual
      else
         write(line,'(i7,2(es12.4))')it,t,dt
      endif
   else
      if(residmin>smalldouble) then
         write(line,'(i7,1(es12.4))')it,residual
      else
         write(line,'(i7)')it
      endif
   end if
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                       MPI_CHARACTER,status,ierrmpi)
   do iw=1,nw
      write(line,'(es12.4)')wmean(iw)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   do level=1,mxnest
      write(line,'(es12.4)')volumeflat(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do

   do level=1,mxnest
      write(line,'(i8)') nleafs_level(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do

   write(line,'(a3,6(es10.2))') ' | ', Xload, Xmemory, cellupdatesPerSecond, &
        activeBlocksPerCore, wctPerCodeTime, timeToFinish
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)

end if

end subroutine printlog_default
!=============================================================================
{#IFDEF HDF5
subroutine write_snapshot_hdf5
use mod_forest
use HDF5
use mod_metric, only: CoordToCart
include 'amrvacdef.f'

integer :: igrid, Morton_no,error,i,j,k,d,counter,offset,nx^D,&
rank,alloc_status,num_blocks_global,num_blocks_local,ierr,naux
integer, dimension(npe) :: morton_counts
integer :: num_datasets=1
integer, dimension(:), allocatable :: num_variables
character(len=80) :: filename,attributename,datasetname
character(len=80), dimension(:), allocatable :: dataset_names
logical, save :: firstsnapshot=.true.

integer :: ix^L

!Buffers for cartesian data, they will only be allocated if the
!grid is non-cartesian.
double precision, dimension(:^D&,:), allocatable :: xCart, wCart
double precision, dimension(:), allocatable :: normconv

!Buffers for aggregating the data. They will only be allocated
!if fastIO is set to true.
integer, dimension(:), allocatable :: buffer_levels
double precision, dimension(:,:,:), allocatable :: buffer_coords
double precision, dimension(:^D&,:,:), allocatable :: buffer_cc_data
{#IFDEF STAGGERED
double precision, dimension(:^D&,:,:), allocatable :: buffer_st_data
\}
integer :: block_idx, ntot, save_GZ_HDF5

!HDF5 variables
integer(hid_t) :: plist_acc_id,plist_wr_id,file_id,atype_id,attribute
integer(hid_t) :: dspace_scalar,dspace_eqpar
integer(hid_t) :: dspace_blocks,dspace_blocks_grid,dspace_blocks_tmp
integer(hid_t), dimension(:), allocatable :: dspace_blocks_grid_tmp
integer(hid_t), dimension(2) :: dspaces_vars_blocks_nx,dsets_celldata,&
dspaces_vars_blocks_nx_tmp,memspaces_vars_blocks_nx
integer(hid_t) :: memspace_blocks, memspace_blocks_grid
integer(hid_t) :: dset_levels,dset_grid,dset_test
integer(hsize_t), dimension(5) :: stride = (/1,1,1,1,1/)
integer(hsize_t), dimension(5) :: dims_count = (/1,1,1,1,1/)
integer(hsize_t), dimension(1) :: dims
integer(hsize_t), dimension(5) :: block,dims_start
double precision, dimension(12,12,12) :: carts1,carts2,carts3
!-----------------------------------------------------------------------------
if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

if (convert) snapshot = snapshot-1

num_blocks_global=nleafs
num_blocks_local=Morton_stop(mype)-Morton_start(mype)+1

{#IFDEF STAGGERED
num_datasets=2
\}

allocate(num_variables(1:num_datasets), stat=alloc_status)
allocate(dataset_names(1:num_datasets), stat=alloc_status)
num_variables(1) = nw
dataset_names(1) = "Cell-centered Variables"
{#IFDEF STAGGERED
num_variables(2) = nws
dataset_names(2) = "Staggered Variables"
\}
naux = nwauxio
num_variables(1) = num_variables(1) + naux

allocate(dspace_blocks_grid_tmp(1:ndim), stat=alloc_status)

! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".hdf5"

call H5open_f(error)

!Create new file
  call H5Pcreate_f(H5P_FILE_ACCESS_F,plist_acc_id,error)
  call H5Pset_fapl_mpio_f(plist_acc_id,MPI_COMM_WORLD,MPI_INFO_NULL,error)
  call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error,access_prp=plist_acc_id)
  call H5Pclose_f(plist_acc_id,error)

!First write metadata as attributes in root-directory
!-----------------------------------------------------------------------
!Create the dataspaces for the attribute
call H5Screate_f(H5S_SCALAR_f, dspace_scalar, error)
dims_count(1)=neqpar+nspecialpar
call H5Screate_simple_f(1,dims_count, dspace_eqpar, error)

!Write iteration number
dims=1
call H5Acreate_f(file_id, "Iteration", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, it, dims, error)
call H5Aclose_f(attribute, error)

!Write simulation time
call H5Acreate_f(file_id, "Time", H5T_IEEE_F64BE, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_DOUBLE, t, dims, error)
call H5Aclose_f(attribute, error)

!Write extent of grid in x^D-direction
{nx^D=ixGhi^D-ixGlo^D+1
write(attributename,"(a,i1.0)") "ExtentX",^D
call H5Acreate_f(file_id,attributename,H5T_NATIVE_INTEGER,dspace_scalar,&
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, nx^D, dims, error)
call H5Aclose_f(attribute, error)\}

!Write equation parameters
dims=neqpar+nspecialpar
call H5Acreate_f(file_id, "EqPar", H5T_IEEE_F64BE, dspace_eqpar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_DOUBLE, eqpar, dims, error)
call H5Aclose_f(attribute, error)

!Write number of active grids
dims=1
call H5Acreate_f(file_id, "Nleafs", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, nleafs, dims, error)
call H5Aclose_f(attribute, error)

!Write max refinement level
call H5Acreate_f(file_id, "LevMax", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, levmax, dims, error)
call H5Aclose_f(attribute, error)

!Write number of ghost-zones
call H5Acreate_f(file_id, "GhostZones", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, dixB, dims, error)
call H5Aclose_f(attribute, error)

!Write if ghost-zones are included in output
!Note that HDF5 doesnt support Bools. We use 0=False, 1=True
if(save_GZ) then
  save_GZ_HDF5=1
else
  save_GZ_HDF5=0
end if
call H5Acreate_f(file_id, "HasGhostZones",H5T_NATIVE_INTEGER,dspace_scalar,&
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, save_GZ_HDF5, dims, error)
call H5Aclose_f(attribute, error)

!Write number of dimensions
call H5Acreate_f(file_id, "Dims", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, ndim, dims, error)
call H5Aclose_f(attribute, error)

!Write number of vector components
call H5Acreate_f(file_id, "VectorComponents", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, ndir, dims, error)
call H5Aclose_f(attribute, error)

!Write number of cell-centered variables
call H5Acreate_f(file_id, "nw", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, nw, dims, error)
call H5Aclose_f(attribute, error)

!Write number of staggered variables
call H5Acreate_f(file_id, "nws", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, nws, dims, error)
call H5Aclose_f(attribute, error)

!Write number of equation parameters
call H5Acreate_f(file_id, "NEqPar", H5T_NATIVE_INTEGER, dspace_scalar, &
                 attribute, error)
call H5Awrite_f(attribute, H5T_NATIVE_INTEGER, neqpar+nspecialpar, dims, error)
call H5Aclose_f(attribute, error)

!Close dataspaces 
call H5Sclose_f(dspace_scalar,error)
call H5Sclose_f(dspace_eqpar, error)
!-----------------------------------------------------------------------------------

!Create dataspaces for datasets
dims_count(1) = num_blocks_global
call H5Screate_simple_f(1,dims_count, dspace_blocks, error)
dims_count(1) = ndim
if(save_GZ .eqv. .True.) then
  dims_count(2) = nx^D*
else
  dims_count(2) = {(nx^D - 2*dixB+1)*&\}
  1 !This is a terrrible hack. Check the compiled F90source to understand this!
end if
dims_count(3) = num_blocks_global
call H5Screate_simple_f(3, dims_count, dspace_blocks_grid, error)
dims_count(ndim+1) = num_variables(1)
dims_count(ndim+2) = num_blocks_global
if(save_GZ .eqv. .True.) then
   {dims_count(^D) = nx^D\}
else
   {dims_count(^D) = (nx^D - 2*dixB+1)\}
end if
rank = ndim+2
call H5Screate_simple_f(rank,dims_count,dspaces_vars_blocks_nx(1),error)
{#IFDEF STAGGERED
dims_count(ndim+1) = num_variables(2)
dims_count(ndim+2) = num_blocks_global
{dims_count(^D) = 1+dims_count(^D)\}
call H5Screate_simple_f(rank,dims_count,dspaces_vars_blocks_nx(2),error)
\}


!Create the datasets
call H5Dcreate_f(file_id, "Levels", H5T_STD_I32BE, dspace_blocks, &
            dset_levels, error)
call H5Dcreate_f(file_id,"Grid", H5T_IEEE_F64BE,dspace_blocks_grid,&
            dset_grid,error)
call H5Dcreate_f(file_id, dataset_names(1), H5T_IEEE_F64BE, &
                 dspaces_vars_blocks_nx(1), dsets_celldata(1), error)
{#IFDEF STAGGERED
call H5Dcreate_f(file_id, dataset_names(2), H5T_IEEE_F64BE, &
                 dspaces_vars_blocks_nx(2), dsets_celldata(2), error)
\}

!Create MPI property list.
call H5Pcreate_f(H5P_DATASET_XFER_F,plist_wr_id,error)
call H5Pset_dxpl_mpio_f(plist_wr_id, H5FD_MPIO_COLLECTIVE_F, error)
 
if(fastIO .eqv. .true.) then

  !Prepare dataspaces for describing data in memory
  igrid = sfc_to_igrid(Morton_start(mype))
  dims=num_blocks_local
  call H5Screate_simple_f(1,dims,memspace_blocks,error)
  dims_count(1) = num_blocks_local
  dims_count(3) = ndim
  if(save_GZ .eqv. .True.) then
    ntot = nx^D*
    {ixmin^D = 1\}
    {ixmax^D = nx^D\}
  else
    ntot = {(nx^D - 2*dixB + 1)*& \}
    1
    {ixmin^D = 1+dixB\}
    {ixmax^D = nx^D-dixB+1\}
  end if
  dims_count(2) = ntot
  rank=3
  call H5Screate_simple_f(rank, dims_count, memspace_blocks_grid, error)
  {dims_count(^D+2)=ixmax^D-ixmin^D+1\}
  dims_count(2)=num_variables(1)
  dims_count(1)=num_blocks_local
  rank=ndim+2
  call H5Screate_simple_f(rank, dims_count, memspaces_vars_blocks_nx(1), error)
  {#IFDEF STAGGERED
  {dims_count(^D+2)=dims_count(^D+2)+1\}
  dims_count(2)=num_variables(2)
  call H5Screate_simple_f(rank, dims_count, memspaces_vars_blocks_nx(2), error)
  \}

  !Allocate data_buffers
  allocate(buffer_levels(1:num_blocks_local), stat=alloc_status)
  allocate(buffer_coords(1:ndim,1:ntot,1:num_blocks_local), stat=alloc_status)
  allocate(buffer_cc_data(pw(igrid)%ixG^S,1:num_variables(1),1:num_blocks_local), stat=alloc_status)
  {#IFDEF STAGGERED
  allocate(buffer_st_data(pws(igrid)%ixG^S,1:nws,1:num_blocks_local), stat=alloc_status)
  \}

  do Morton_no=Morton_start(mype), Morton_stop(mype)
  
    igrid=sfc_to_igrid(Morton_no)
    call set_tmpGlobals(igrid)
    block_idx = Morton_no - Morton_start(mype) + 1

    allocate(wCart(px(igrid)%ixG^S,1:(num_variables(1))), stat = alloc_status) 
    allocate(xCart(px(igrid)%ixG^S,1:ndim), stat = alloc_status) 
    allocate(normconv(0:num_variables(1)), stat = alloc_status)

    wCart(:^D&,1:nw) = pw(igrid)%w
    !call specialvar_output(px(igrid)%ixG^L,ixM^LL^LADD1,num_variables(1),wCart,ps(igrid),normconv)
      ! call specialvar_output(px(igrid)%ixG^L,ixM^LL^LADD1,num_variables(1),wCart,ps(igrid),normconv, &
      !                       !pgeo(igrid)%dx, node(plevel_,igrid))
      !                       pgeo(igrid)%dx, node(plevel_,igrid), psold(igrid))
    !M1_test
      call specialvar_output(px(igrid)%ixG^L,ixM^LL^LADD1,num_variables(1),wCart,ps(igrid),normconv, &
      !pgeo(igrid)%dx, node(plevel_,igrid))
      pgeo(igrid)%dx, node(plevel_,igrid), psold(igrid), psm1(igrid))

    do d=1,ndim
      xCart(:^D&,d) = px(igrid)%x(:^D&,d) - dx(d,node(plevel_,igrid))/2.
    end do

    deallocate(normconv, stat=alloc_status)

    !With HDF5 we always save prim in accordance with the input routine

!stop 'why call c2p'
!    call primitive(px(igrid)%ixG^L,px(igrid)%ixG^L,wCart, &
!                   px(igrid)%x)

    !First coordinates and vectors need to be converted to cartesian for
    !data visualization to work.
    if(.not. convert_nocartesian) then
      call CoordToCart(px(igrid)%ixG^L, px(igrid)%ixG^L, xCart, xCart)
    end if

    !Aggregate data
    buffer_levels(block_idx) = node(plevel_,igrid)
    do d=1,ndim
      buffer_coords(d,:,block_idx) = reshape(xCart(ix^S,d), (/ ntot /))
    end do
    buffer_cc_data(:^D&,:,block_idx) = wCart(:^D&,1:num_variables(1))
    {#IFDEF STAGGERED
    buffer_st_data(:^D&,:,block_idx) = pws(igrid)%w
    \}

    deallocate(xCart, stat = alloc_status) 
    deallocate(wCart, stat = alloc_status) 
 
  end do
  
  !Prepare local hyperslabbed dataspaces for writing datasets to file
  dims_count=(/1,1,1,1,1/)
  dims_start(1)=Morton_start(mype)-1
  block(1)=num_blocks_local
  call H5Sselect_hyperslab_f(dspace_blocks, H5S_SELECT_SET_F,&
                             dims_start,dims_count,error,stride,block)
  dims_start(2)=0
  dims_start(3)=Morton_start(mype)-1
  block(3)=num_blocks_local
  block(2)=ntot
  block(1)=ndim
  dims_start(1)=0
  call H5dget_space_f(dset_grid, dspace_blocks_grid_tmp(1), error)
  call H5Sselect_hyperslab_f(dspace_blocks_grid_tmp(1), H5S_SELECT_SET_F,&
                             dims_start,dims_count,error,stride,block)
  dims_start(ndim+2)=Morton_start(mype)-1
  block(ndim+2)=num_blocks_local
  if(save_GZ .eqv. .True.) then
    {dims_start(^D) = 0
    block(^D) = nx^D\}
  else
    {dims_start(^D) = 0
    block(^D) = nx^D-2*dixB+1\}
  end if
  dims_start(ndim+1) = 0
  block(ndim+1) = num_variables(1)
  call H5Sselect_hyperslab_f(dspaces_vars_blocks_nx(1),&
                   H5S_SELECT_SET_F,dims_start,dims_count, error, stride, block)
  {#IFDEF STAGGERED
  {block(^D) = 1+block(^D)\}
  block(ndim+1) = num_variables(2)
  call H5Sselect_hyperslab_f(dspaces_vars_blocks_nx(2),&
                   H5S_SELECT_SET_F,dims_start,dims_count, error, stride, block)
  \}

  !Dump data
  dims=num_blocks_local
  call H5Dwrite_f(dset_levels,H5T_NATIVE_INTEGER,buffer_levels,dims,error,&
                  memspace_blocks,dspace_blocks,plist_wr_id)

  dims_count(1)=num_blocks_local
  dims_count(3)=ndim
  dims_count(2) = ntot

  call H5Dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, buffer_coords, dims_count, &
                  error,memspace_blocks_grid, dspace_blocks_grid_tmp(1), plist_wr_id)

  
  {dims_count(^D+2) = ixmax^D-ixmin^D+1\}
  dims_count(2) = num_variables(1)
  dims_count(1)=num_blocks_local
  call H5Dwrite_f(dsets_celldata(1), H5T_NATIVE_DOUBLE, buffer_cc_data(ix^S,:,:), dims_count, error, &
                  memspaces_vars_blocks_nx(1),dspaces_vars_blocks_nx(1), plist_wr_id)
  {#IFDEF STAGGERED
  {dims_count(^D+2) = 1+dims_count(^D+2)\}
  {ixmax^D=ixmax^D+1\}
  dims_count(2) = num_variables(2)
  call H5Dwrite_f(dsets_celldata(2), H5T_NATIVE_DOUBLE, buffer_st_data(ix^S,:,:), dims_count, error, &
                  memspaces_vars_blocks_nx(2),dspaces_vars_blocks_nx(2), plist_wr_id)
  \}
  
  !Deallocate data buffers
  deallocate(buffer_levels, stat=alloc_status)
  deallocate(buffer_coords, stat=alloc_status)
  deallocate(buffer_cc_data, stat=alloc_status)
  {#IFDEF STAGGERED
  deallocate(buffer_st_data, stat=alloc_status)
  \}

  !Close temporary dataspaces
  call H5Sclose_f(dspace_blocks_grid_tmp(1), error)
 
else
 
  !Prepare dataspaces for describing data in memory
  dims=1
  call H5Screate_simple_f(1,dims,memspace_blocks,error)
  if(save_GZ .eqv. .True.) then
    {dims_count(^D) = nx^D\}
  else
    {dims_count(^D) = nx^D - 2*dixB + 1\}
  end if
  rank=ndim
  call H5Screate_simple_f(rank, dims_count, memspace_blocks_grid, error)
  dims_count(ndim+1)=num_variables(1)
  rank=ndim+1
  call H5Screate_simple_f(rank, dims_count, memspaces_vars_blocks_nx(1), error)
  {#IFDEF STAGGERED
  {dims_count(^D)=1+dims_count(^D)\}
  dims_count(ndim+1)=num_variables(2)
  rank=ndim+1
  call H5Screate_simple_f(rank, dims_count, memspaces_vars_blocks_nx(2), error)
  \}
  
  !This is neccessary so that each process knows when to write collectively and when not
  call MPI_ALLGATHER(Morton_stop(mype)-Morton_start(mype),1,MPI_INT,morton_counts,1,MPI_INT,MPI_COMM_WORLD,ierr)
  
  !Loop over all local blocks
  do Morton_no=Morton_start(mype), Morton_stop(mype)

    igrid=sfc_to_igrid(Morton_no)
    call set_tmpGlobals(igrid)

    allocate(wCart(px(igrid)%ixG^S,1:(num_variables(1))), stat = alloc_status) 
    allocate(xCart(px(igrid)%ixG^S,1:ndim), stat = alloc_status) 
    allocate(normconv(0:num_variables(1)), stat = alloc_status)
    

    wCart(:^D&,1:nw) = pw(igrid)%w 
      ! M1_test
      !call specialvar_output(px(igrid)%ixG^L,ixM^LL^LADD1,num_variables(1),wCart,ps(igrid),normconv, &
      !                       pgeo(igrid)%dx, node(plevel_,igrid), psold(igrid))
      call specialvar_output(px(igrid)%ixG^L,ixM^LL^LADD1,num_variables(1),wCart,ps(igrid),normconv, &
                             pgeo(igrid)%dx, node(plevel_,igrid), psold(igrid), psm1(igrid))
    do d=1,ndim
      xCart(:^D&,d) = px(igrid)%x(:^D&,d) - dx(d,node(plevel_,igrid))/2.
    end do

    deallocate(normconv, stat=alloc_status)
 
    !With HDF5 we always save prim in accordance with the input routine
    !call primitive(px(igrid)%ixG^L,px(igrid)%ixG^L,wCart(:^D&,1:nw), &
    !               px(igrid)%x)

    !First coordinates and vectors need to be converted to cartesian for
    !data visualization to work.
    if(.not. convert_nocartesian) then
      call CoordToCart(px(igrid)%ixG^L, px(igrid)%ixG^L, xCart, xCart)
    end if
 
    !Copy Dataspace for hyperslabbing
    call H5dget_space_f(dset_levels, dspace_blocks_tmp, error)
    do d=1,ndim
      call H5dget_space_f(dset_grid, dspace_blocks_grid_tmp(d), error)
    end do
    call H5dget_space_f(dsets_celldata(1), dspaces_vars_blocks_nx_tmp(1), error)
    {#IFDEF STAGGERED
    call H5dget_space_f(dsets_celldata(2), dspaces_vars_blocks_nx_tmp(2), error)
    \}
  
    !Prepare local hyperslabbed dataspaces for writing datasets to file
    dims_count=(/1,1,1,1,1/)
    dims_start(1) = Morton_no-1
    block(1) = 1
    call H5Sselect_hyperslab_f(dspace_blocks_tmp, H5S_SELECT_SET_F,&
                               dims_start,dims_count,error,stride,block)
    dims_start(2)=0
    dims_start(3)=Morton_no-1
    block(3) = 1
    if(save_GZ .eqv. .True.) then
      block(2) = nx^D*
    else
      block(2) = {(nx^D - 2*dixB+1)*&\}
      1. !This is a terrible hack, but so be it. Check the compiled source file to understand this!
    end if
    block(1) = 1
    do d=1,ndim
      dims_start(1) = d-1
      call H5Sselect_hyperslab_f(dspace_blocks_grid_tmp(d), H5S_SELECT_SET_F,&
                                 dims_start,dims_count, error, stride, block)
    end do
    dims_start(ndim+2)=Morton_no-1
    block(ndim+2)=1
    if(save_GZ .eqv. .True.) then
      {dims_start(^D) = 0
      block(^D) = nx^D\}
    else
      {dims_start(^D) = 0
      block(^D) = nx^D-2*dixB+1\}
    end if
    dims_start(ndim+1) = 0
    block(ndim+1) = num_variables(1)
    call H5Sselect_hyperslab_f(dspaces_vars_blocks_nx_tmp(1),&
                     H5S_SELECT_SET_F,dims_start,dims_count, error, stride, block)
    {#IFDEF STAGGERED
    {dims_start(^D) = 0
    block(^D) = 1+block(^D)\}
    block(ndim+1) = num_variables(2)
    call H5Sselect_hyperslab_f(dspaces_vars_blocks_nx_tmp(2),&
                     H5S_SELECT_SET_F,dims_start,dims_count, error, stride, block)
    \}
  
    !CHeck if independent write-calls will be necessary. If so, change plist_wr_id
    if(maxval(morton_counts)==(Morton_no-Morton_start(mype))&
       .and. maxval(morton_counts)/=minval(morton_counts))then
      !Change data_transfer_list to independent
      call H5Pcreate_f(H5P_DATASET_XFER_F,plist_wr_id,error)
      call H5Pset_dxpl_mpio_f(plist_wr_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if
  
    !Dump data
    dims=1
    call H5Dwrite_f(dset_levels, H5T_NATIVE_INTEGER, node(plevel_,igrid), dims, error,&
                    memspace_blocks, dspace_blocks_tmp, plist_wr_id)
  
    if(save_GZ .eqv. .True.) then
      {dims_count(^D) = nx^D\}
      dims_count(ndim+1)=1
      {ixmin^D = 1\}
      {ixmax^D = nx^D\}
    else
      {dims_count(^D) = nx^D-2*dixB+1\}
      dims_count(ndim+1)=1
      {ixmin^D = 1+dixB\}
      {ixmax^D = nx^D-dixB+1\}
    end if

    do d=1,ndim
      call H5Dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xCart(ix^S,d), dims_count, &
                      error,memspace_blocks_grid, dspace_blocks_grid_tmp(d), plist_wr_id)
    end do
  
  
    dims_count(ndim+1) = num_variables(1)
    call H5Dwrite_f(dsets_celldata(1), H5T_NATIVE_DOUBLE, wCart(ix^S,:), dims_count, error, &
                    memspaces_vars_blocks_nx(1),dspaces_vars_blocks_nx_tmp(1), plist_wr_id)
    
    {#IFDEF STAGGERED
    {dims_count(^D) = 1+dims_count(^D)\}
    {ixmax^D=ixmax^D+1\}
    dims_count(ndim+1) = num_variables(2)
    call H5Dwrite_f(dsets_celldata(2), H5T_NATIVE_DOUBLE, pws(igrid)%w(ix^S,:), dims_count, error, &
                    memspaces_vars_blocks_nx(2),dspaces_vars_blocks_nx_tmp(2), plist_wr_id)
    \}
  
    call H5Sclose_f(dspace_blocks_tmp, error)
    do d=1,ndim
      call H5Sclose_f(dspace_blocks_grid_tmp(d), error)
    end do
    call H5Sclose_f(dspaces_vars_blocks_nx_tmp(1), error)
    {#IFDEF STAGGERED 
    call H5Sclose_f(dspaces_vars_blocks_nx_tmp(2), error)
    \}
  
    deallocate(xCart, stat = alloc_status) 
    deallocate(wCart, stat = alloc_status) 
  end do
end if

!Close property list
call H5Pclose_f(plist_wr_id, error)

!Close dataspaces for describing the memory
call H5Sclose_f(memspace_blocks, error)
call H5Sclose_f(memspace_blocks_grid, error)
call H5Sclose_f(memspaces_vars_blocks_nx(1), error)
{#IFDEF STAGGERED
call H5Sclose_f(memspaces_vars_blocks_nx(2), error)
\}

!Close dataspaces for describing the file
call H5Sclose_f(dspace_blocks, error)
call H5Sclose_f(dspace_blocks_grid, error)
call H5Sclose_f(dspaces_vars_blocks_nx(1), error)
{#IFDEF STAGGERED
call H5Sclose_f(dspaces_vars_blocks_nx(2), error)
\}

!Close datasets
call H5Dclose_f(dset_levels, error)
call H5Dclose_f(dset_grid, error)
call H5Dclose_f(dsets_celldata(1), error)
{#IFDEF STAGGERED
call H5Dclose_f(dsets_celldata(2), error)
\}

!-----------------------------------------------------------------------------------
call H5Fclose_f(file_id,error)
call H5close_f(error)
deallocate(dspace_blocks_grid_tmp, stat=alloc_status)

snapshot=snapshot+1

if(mype==0 .and. write_xdmf .eqv. .true.) then
  call MakeXDMF(filename,num_variables,num_datasets,dataset_names)
end if

end subroutine write_snapshot_hdf5
!=============================================================================
subroutine MakeXDMF(fname,num_variables,num_datasets,dataset_names)
use mod_forest
include 'amrvacdef.f'

integer, intent(in) :: num_datasets
character(len=80), dimension(num_datasets) :: dataset_names
integer, dimension(num_datasets) :: num_variables
character(len=10), dimension(nw+nwauxio) :: variable_names
character(len=80) :: fname
character(len=80) :: fname2
integer, parameter :: out_unit=20
integer :: igrid, Morton_no, n_block, nx1, nx2, nx3, n_quantity, n_variable,index, ixfile
integer :: ixmin1, ixmin2, ixmin3, ixmax1, ixmax2, ixmax3
character(len=10)   :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead,varnames_tmp
!-----------------------------------------------------------------------------
! Initialize:
nx1=0;nx2=0;nx3=0
ixmin1=0;ixmin2=0;ixmin3=0

!Get variable names
call getheadernames(wnamei,xandwnamei,outfilehead) ! appends to primnames
varnames_tmp=primnames

n_quantity=1
do while (varnames_tmp /= "")
  index=SCAN(varnames_tmp,' ')
  variable_names(n_quantity)=varnames_tmp(1:index-1)
  varnames_tmp=varnames_tmp(index+1:)
  n_quantity=n_quantity+1
end do


{nx^D=ixGhi^D-ixGlo^D+1\}

if(save_GZ .eqv. .False.) then
  {nx^D = nx^D - 2 * dixB+1\}
  {ixmin^D = 0\} 
  {ixmax^D = nx^D\} 
else
  {ixmin^D=0\}
  {ixmax^D=nx^D\}
end if

if(xdmf_type=='Cell' .and. (.not. save_GZ)) then
  {ixmax^D = ixmax^D -1 \}
end if

do ixfile=len_trim(fname),1,-1
   if (fname(ixfile:ixfile) .eq. '/') exit
end do
ixfile=ixfile+1

write(fname2,'(a,i4.4,a)') TRIM(filenameout), snapshot-1, ".xdmf"

open(unit=out_unit,file=fname2,action="write",status="replace")

!Write header
write(out_unit,"(a)") '<?xml version="1.0" ?>'
write(out_unit,"(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
write(out_unit,"(a)") '<Xdmf Version="2.0">'
write(out_unit,"(a)") '<Domain>'
write(out_unit,"(a)") '<Grid Name="Mesh" GridType="Collection">'

do n_block=0, nleafs-1
  
  write(out_unit,"(a,I0,a)") '  <Grid Name="MeshBlock', n_block, '" GridType="Uniform">'

  !Write topology
  select case(ndim)
    case(3)
      write(out_unit,"(a,I0,a,I0,a,I0,a)")&
      '    <Topology TopologyType="3DSMesh" NumberOfElements="', nx3,' ',nx2,' ',nx1, '"/>'
    case(2)
      write(out_unit,"(a,I0,a,I0,a)")&
      '    <Topology TopologyType="2DSMesh" NumberOfElements="', nx2,' ',nx1, '"/>'
    case default
      write(*,*) "XDMF NOT YET IMPLEMENTED FOR 1D."
    end select

  !Write geometry
  select case(ndim)
    case(3)
      write(out_unit,"(a)") '    <Geometry GeometryType="XYZ">'
   case(2)
      write(out_unit,"(a)") '    <Geometry GeometryType="XY">'
    case default
      write(*,*) "XDMF NOT YET IMPLEMENTED FOR 1D."
    end select

  write(out_unit,"(a,I0,a,I0,a)")&  
  '      <DataItem ItemType="HyperSlab" Dimensions="',nx^D*,&
  ' ',ndim,'">'
  write(out_unit,"(a,I0,a,I0,a,I0,a)")& 
  '      <DataItem Dimensions="3 3" NumberType="Int"> ',&
             n_block,' 0 0 1 1 1 1 ',nx^D*,' ',ndim,' </DataItem>'
  write(out_unit,"(a,I0,a,I0,a,I0,a,a,a)") '        <DataItem Dimensions="',nleafs,&
                         ' ',nx^D*,' ',ndim,'" Format="HDF"> ',TRIM(fname(ixfile:)),':/Grid </DataItem>'
  write(out_unit,"(a)") '      </DataItem>'

  
  write(out_unit,"(a)") '    </Geometry>'

  !Write description of cell-centered data
  n_quantity=1
  do n_variable=0,num_variables(1)-1

    if(xdmf_type=='Cell') then 
      write(out_unit,"(a,a,a)") '    <Attribute Name="',TRIM(variable_names(n_quantity)),&
                                '" Center="Cell">'
    else
      write(out_unit,"(a,a,a)") '    <Attribute Name="',TRIM(variable_names(n_quantity)),&
                                '" Center="Node">'
    end if
    n_quantity=n_quantity+1 
    
    if(nx3>1) then
      write(out_unit,"(a,I0,a,I0,a,I0,a)") '      <DataItem ItemType="HyperSlab" Dimensions="',&
                                           ixmax3,' ',ixmax2,' ',ixmax1,'">'
      write(out_unit,"(a,I0,a,I0,a,I0,a,I0,a,I0,a,I0,a,I0,a,I0,a)")&
      '        <DataItem Dimensions="3 5" NumberType="Int"> ',n_block,' ',n_variable,&
      ' ',ixmin3,' ',ixmin2,' ',ixmin1 ,' 1 1 1 1 1 1 1 ',ixmax3,' ',ixmax2,' ',ixmax1,' </DataItem>'
    else
      write(out_unit,"(a,I0,a,I0,a)") '      <DataItem ItemType="HyperSlab" Dimensions="',&
                                           ixmax2,' ',ixmax1,'">'
      write(out_unit,"(a,I0,a,I0,a,I0,a,I0,a,I0,a,I0,a)")&
      '        <DataItem Dimensions="3 4" NumberType="Float"> ',n_block,' ',n_variable,&
      ' ',ixmin2,' ',ixmin1,' 1 1 1 1 1 1 ',ixmax2,' ',ixmax1,' </DataItem>'
    end if
    
    if(nx3>1) then
        write(out_unit,"(a,I0,a,I0,a,I0,a,I0,a,I0,a,a,a,a,a)")& 
        '        <DataItem Dimensions="',nleafs,' ',num_variables(1),' ',ixmax3,' ',ixmax2,' ',ixmax1,&
        '" Format="HDF"> ',TRIM(fname(ixfile:)),':/',TRIM(dataset_names(1)),' </DataItem>'
        write(out_unit,"(a)") '      </DataItem>'
        write(out_unit,"(a)") '    </Attribute>'
    else
        write(out_unit,"(a,I0,a,I0,a,I0,a,I0,a,a,a,a,a)")& 
        '        <DataItem Dimensions="',nleafs,' ',num_variables(1),' ',ixmax2,' ',ixmax1,&
        '" Format="HDF"> ',TRIM(fname(ixfile:)),':/',TRIM(dataset_names(1)),' </DataItem>'
        write(out_unit,"(a)") '      </DataItem>'
        write(out_unit,"(a)") '    </Attribute>'
    end if 

  end do
  
  !End block
  write(out_unit,"(a)") '  </Grid>'
end do

!Complete header
write(out_unit,"(a)") '</Grid>' 
write(out_unit,"(a)") '</Domain>' 
write(out_unit,"(a)") '</Xdmf>' 

close(out_unit)

end subroutine MakeXDMF
!=================================================================================================
subroutine read_snapshot_hdf5
use mod_forest
use HDF5
include 'amrvacdef.f'

integer :: error, igrid, Morton_no,levmaxini, ndimini, ndirini, nwini,rank,&
{#IFDEF STAGGERED nwsini,} neqparini, nxini^D, ierr, alloc_stat, HasGZ, GZ
integer :: ix^L
integer(hsize_t), dimension(3) :: dims
integer(hsize_t), dimension(4) :: dims4
integer(hsize_t), dimension(5) :: dims_count,dims_start
integer(hsize_t), allocatable, dimension(:) :: stride
integer(hsize_t), allocatable, dimension(:) :: block
character(len=80) :: filename
logical :: fexist
!logical, allocatable, dimension(:^D&) :: patchw

!HDF5 variables
integer(hid_t) :: file_id, attribute_id, dspaces_vars_blocks_nx_1,&
dsets_celldata1,memspaces_vars_blocks_nx_1
{#IFDEF STAGGERED
integer(hid_t) :: dspaces_vars_blocks_nx_2,dsets_celldata2,&
memspaces_vars_blocks_nx_2\}
integer(hid_t), dimension(2) :: dsets_celldata
!------------------------------------------------------------------------------------------------

write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".hdf5"

if(mype==0) then
  inquire(file=filename,exist=fexist)
  if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found!")
endif

call H5open_f(error)
call H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

call H5Aopen_f(file_id, "Nleafs", attribute_id, error)
dims(1) = 1
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, nleafs, dims, error)
call H5Aclose_f(attribute_id, error)
nleafs_active=nleafs

call H5Aopen_f(file_id, "LevMax", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, levmaxini, dims, error)
call H5Aclose_f(attribute_id, error)

call H5Aopen_f(file_id, "Dims", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, ndimini, dims, error)
call H5Aclose_f(attribute_id, error)

call H5Aopen_f(file_id, "VectorComponents", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, ndirini, dims, error)
call H5Aclose_f(attribute_id, error)

call H5Aopen_f(file_id, "nw", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, nwini, dims, error)
call H5Aclose_f(attribute_id, error)

{#IFDEF STAGGERED
call H5Aopen_f(file_id, "nws", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, nwsini, dims, error)
call H5Aclose_f(attribute_id, error)
\}

call H5Aopen_f(file_id, "NEqPar", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, neqparini, dims, error)
call H5Aclose_f(attribute_id, error)

call H5Aopen_f(file_id, "Iteration", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, it, dims, error)
call H5Aclose_f(attribute_id, error)

call H5Aopen_f(file_id, "Time", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_DOUBLE, t, dims, error)
call H5Aclose_f(attribute_id, error)

call H5Aopen_f(file_id, "HasGhostZones", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, HasGZ, dims, error)
call H5Aclose_f(attribute_id, error)

call H5Aopen_f(file_id, "GhostZones", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, GZ, dims, error)
call H5Aclose_f(attribute_id, error)


! check if settings are suitable for restart
if (levmaxini>mxnest) then
   if (mype==0) write(*,*) "number of levels in restart file = ",levmaxini
   if (mype==0) write(*,*) "mxnest = ",mxnest
   call mpistop("mxnest should be at least number of levels in restart file")
end if
if (ndimini/=ndim) then
   if (mype==0) write(*,*) "ndim in restart file = ",ndimini
   if (mype==0) write(*,*) "ndim = ",ndim
   call mpistop("reset ndim to ndim in restart file")
end if
if (ndirini/=ndir) then
   if (mype==0) write(*,*) "ndir in restart file = ",ndirini
   if (mype==0) write(*,*) "ndir = ",ndir
   call mpistop("reset ndir to ndir in restart file")
end if
if (nw/=nwini) then
   if (mype==0) write(*,*) "nw=",nw," and nw in restart file=",nwini
   call mpistop("currently, changing nw at restart is not allowed")
end if
if (GZ/=dixB) then
   if (mype==0) write(*,*) "number of ghostzones=",dixB," and number of ghostzones in restart file=",GZ
   call mpistop("currently, changing number of ghostzones at restart is not allowed with HDF5")
end if

!{#IFDEF STAGGERED
!if (nws/=nwsini) then
!   if (mype==0) write(*,*) "nws=",nws," and nws in restart file=",nwsini
!   call mpistop("currently, changing nws at restart is not allowed")
!end if}

{call H5Aopen_f(file_id, "ExtentX^D", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_INTEGER, nxini^D, dims, error)
call H5Aclose_f(attribute_id, error)\}

if (ixGhi^D/=nxini^D|.or.) then
   if (mype==0) write(*,*) "Error: reset resolution to ",nxini^D
   call mpistop("change with setamrvac")
end if

neqparini=min(neqparini,neqpar+nspecialpar)
dims(1)=neqparini
call H5Aopen_f(file_id, "EqPar", attribute_id, error)
call H5Aread_f(attribute_id, H5T_NATIVE_DOUBLE, eqpar, dims, error)
call H5Aclose_f(attribute_id, error)

call read_forest_hdf5(file_id)

allocate(stride(1:ndimini), stat = alloc_stat)
allocate(block(1:ndimini), stat = alloc_stat)
stride = (/{1^D&,},1,1/)
block = (/{1^D&,},1,1/)

call H5Dopen_f(file_id, "Cell-centered Variables", dsets_celldata1, error)

!Define memspace
if(HasGZ==1) then
  {dims4(^D)=nxini^D\}
else
  {dims4(^D)=nxini^D-2*dixB+1\}
end if
dims4(ndimini+1)=nwini
rank=ndimini+1
call H5Screate_simple_f(rank, dims4, memspaces_vars_blocks_nx_1, error)


{#IFDEF STAGGERED
call H5Dopen_f(file_id, "Staggered Variables", dsets_celldata2, error)

{dims4(^D)= 1+dims4(^D)\}
rank=ndimini+1
dims4(ndimini+1)=nwsini
call H5Screate_simple_f(rank, dims4, memspaces_vars_blocks_nx_2, error)
\}



do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  call alloc_node(igrid)
  call set_tmpGlobals(igrid)

  !Get hyperlab selection of dataset 
  call H5Dget_space_f(dsets_celldata1, dspaces_vars_blocks_nx_1, error)
  call H5Sget_simple_extent_dims_f(dspaces_vars_blocks_nx_1,dims_start,dims_count,error)
  {#IFDEF STAGGERED
  call H5Dget_space_f(dsets_celldata2, dspaces_vars_blocks_nx_2, error)\}

  dims_start(ndimini+2)=Morton_no-1
  dims_count(ndimini+2)=1
  {dims_start(^D) = 0\}
  if(HasGZ==1) then
    {dims_count(^D) = nxini^D\}
    {ixmin^D=1\}
    {ixmax^D=nxini^D\}
  else
    {dims_count(^D) = nxini^D-2*dixB+1\}
    {ixmin^D=dixB+1\}
    {ixmax^D=nxini^D-dixB+1\}
  end if
  dims_start(ndimini+1) = 0
  dims_count(ndimini+1) = nwini
  call H5Sselect_hyperslab_f(dspaces_vars_blocks_nx_1,H5S_SELECT_SET_F,&
                               dims_start,dims_count, error, stride, block)

  {dims4(^D) = dims_count(^D)\}
  dims4(ndimini+1) = nwini
  call H5Dread_f(dsets_celldata1, H5T_NATIVE_DOUBLE, pw(igrid)%w(ix^S,:), dims4,&
       error,memspaces_vars_blocks_nx_1,dspaces_vars_blocks_nx_1)
  
  call H5Sclose_f(dspaces_vars_blocks_nx_1,error)

  {#IFDEF STAGGERED
  {dims_count(^D) = 1+dims_count(^D)\}
  dims_count(ndimini+1) = nwsini
  call H5Sselect_hyperslab_f(dspaces_vars_blocks_nx_2,H5S_SELECT_SET_F,&
                               dims_start,dims_count, error, stride, block)
  dims4(ndimini+1) = nwsini
  {dims4(^D) = 1+dims_count(^D)\}
  {ixmax^D = ixmax^D+1\}
  call H5Dread_f(dsets_celldata2, H5T_NATIVE_DOUBLE, pws(igrid)%w(ix^S,:), dims4,&
       error, memspaces_vars_blocks_nx_2,dspaces_vars_blocks_nx_2)

  call H5Sclose_f(dspaces_vars_blocks_nx_2,error)
  \}

  !call conserve(ixG^LL,ixG^LL,pw(igrid)%w,px(igrid)%x,patchfalse)
 
end do

call H5Sclose_f(memspaces_vars_blocks_nx_1,error)
call H5Dclose_f(dsets_celldata1, error)
{#IFDEF STAGGERED 
call H5Sclose_f(memspaces_vars_blocks_nx_2,error)
call H5Dclose_f(dsets_celldata2, error)\}

call H5Fclose_f(file_id, error)
call H5close_f(error)

deallocate(stride, stat=alloc_stat)
deallocate(block, stat=alloc_stat)

end subroutine read_snapshot_hdf5
\}
!=============================================================================
