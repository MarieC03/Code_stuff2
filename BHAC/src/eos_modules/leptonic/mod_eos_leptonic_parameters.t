!================================================================================
!
!  mod_eos_leptonic_parameters.t
!
!  Parameter module for the 4D leptonic EOS table in BHAC.
!  Stores grid dimensions, axis arrays, variable indices, and global
!  bounds for the three-table architecture:
!
!    alltables_baryon  (nrho x ntemp x nye_baryon)   -- baryonic sector
!    alltables_ele     (nrho x ntemp x nyle)          -- electronic leptonic sector
!    alltables_muon    (nrho x ntemp x nymu)          -- muonic leptonic sector
!
!  The baryon table third axis is yp = ye + ymu (total proton fraction),
!  mapped to the original scollapse/compose ye axis so that existing
!  3D table readers require minimal modification.
!
!  Units throughout: code units (BHAC geometric units) identical to
!  those used in mod_eos_tabulated_parameters, i.e. the conversion
!  factors rho_gf, press_gf, eps_gf are applied at read time.
!
!  Copyright (C) 2024  Harry Ho-Yin Ng
!  Based on BHAC mod_eos_tabulated_parameters by Elias R. Most et al.
!
!================================================================================

module mod_eos_leptonic_parameters
  implicit none
  public

  ! -----------------------------------------------------------------------
  ! Table file names
  ! -----------------------------------------------------------------------
  character(len=256) :: leptonic_table_name  = 'leptonic_eos.h5'
  character(len=256) :: baryon_table_name    = 'baryon_eos.h5'

  !> Which baryon-table format to read:  'scollapse'  or  'compose'
  character(len=32)  :: baryon_table_type    = 'scollapse'

  ! -----------------------------------------------------------------------
  ! Numerical precision / root-finding settings
  ! -----------------------------------------------------------------------
  double precision, parameter :: lep_eos_precision = 1.0d-14
  integer,          parameter :: lep_eos_iter_max  = 50

  ! -----------------------------------------------------------------------
  ! Table dimensions
  ! -----------------------------------------------------------------------
  integer :: lep_nrho  = 0   !< number of density   points (shared by all tables)
  integer :: lep_ntemp = 0   !< number of temperature points (shared by all tables)
  integer :: lep_nye   = 0   !< Ye  axis of the *baryon* sub-table  (yp = ye+ymu)
  integer :: lep_nyle  = 0   !< Ye- axis of the *electronic* sub-table
  integer :: lep_nymu  = 0   !< Ymu axis of the *muonic* sub-table  (log scale)

  ! -----------------------------------------------------------------------
  ! Axis arrays  (allocated in the reader routines)
  ! -----------------------------------------------------------------------
  double precision, allocatable :: lep_logrho_table(:)   ! ln(rho / [code units])
  double precision, allocatable :: lep_logtemp_table(:)  ! ln(T   / [MeV])
  double precision, allocatable :: lep_ye_table(:)       ! yp = yle+ymu  (baryon table)
  double precision, allocatable :: lep_yle_table(:)      ! yle  (linear, electronic table)
  double precision, allocatable :: lep_logymu_table(:)   ! ln(ymu)        (muonic table)

  ! -----------------------------------------------------------------------
  ! Energy shift (as in stellarcollapse tables)
  ! -----------------------------------------------------------------------
  double precision :: lep_energy_shift = 0.0d0

  ! -----------------------------------------------------------------------
  ! Baryon mass
  ! -----------------------------------------------------------------------
  double precision :: lep_baryon_mass  = 1.660539040d-24   ! [g]  default = amu

  ! -----------------------------------------------------------------------
  ! Speed-of-sound flag (0 = Newtonian cs2 stored, 1 = relativistic h*cs2)
  ! -----------------------------------------------------------------------
  integer :: lep_have_rel_cs2 = 0

  ! -----------------------------------------------------------------------
  ! Variable index constants -- BARYON table  (layout: ft(nrho,ntemp,nye,nvars_baryon))
  ! -----------------------------------------------------------------------
  integer, parameter :: lep_nvars_baryon = 13
  integer, parameter :: i_lep_logpress   = 1   ! ln(P)
  integer, parameter :: i_lep_logenergy  = 2   ! ln(eps + e_shift)
  integer, parameter :: i_lep_entropy    = 3   ! s  [kb/baryon]
  integer, parameter :: i_lep_cs2        = 4   ! relativistic cs2 = cs2_nr / h
  integer, parameter :: i_lep_mu_e       = 5   ! mu_e  [MeV]   (baryonic table, for yp>yemax fallback)
  integer, parameter :: i_lep_mu_p       = 6   ! mu_p  [MeV]
  integer, parameter :: i_lep_mu_n       = 7   ! mu_n  [MeV]
  integer, parameter :: i_lep_xa         = 8   ! alpha-particle mass fraction
  integer, parameter :: i_lep_xh         = 9   ! heavy-nucleus mass fraction
  integer, parameter :: i_lep_xn         = 10  ! free neutron mass fraction
  integer, parameter :: i_lep_xp         = 11  ! free proton mass fraction
  integer, parameter :: i_lep_abar       = 12  ! mean mass number
  integer, parameter :: i_lep_zbar       = 13  ! mean proton number

  ! -----------------------------------------------------------------------
  ! Variable index constants -- ELECTRONIC table  (ft(nrho,ntemp,nyle,nvars_ele))
  ! -----------------------------------------------------------------------
  integer, parameter :: lep_nvars_ele    = 9
  integer, parameter :: i_lep_mu_ele     = 1   ! mu_e^leptonic  [MeV]  (includes rest mass)
  integer, parameter :: i_lep_yle_minus  = 2   ! Y_{e^-}
  integer, parameter :: i_lep_yle_plus   = 3   ! Y_{e^+}
  integer, parameter :: i_lep_press_eminus = 4 ! P_{e^-}
  integer, parameter :: i_lep_press_eplus  = 5 ! P_{e^+}
  integer, parameter :: i_lep_eps_eminus   = 6 ! eps_{e^-}
  integer, parameter :: i_lep_eps_eplus    = 7 ! eps_{e^+}
  integer, parameter :: i_lep_s_eminus     = 8 ! s_{e^-}
  integer, parameter :: i_lep_s_eplus      = 9 ! s_{e^+}

  ! -----------------------------------------------------------------------
  ! Variable index constants -- MUONIC table  (ft(nrho,ntemp,nymu,nvars_muon))
  ! -----------------------------------------------------------------------
  integer, parameter :: lep_nvars_muon   = 9
  integer, parameter :: i_lep_mu_mu      = 1   ! mu_mu  [MeV]  (includes rest mass)
  integer, parameter :: i_lep_ymu_minus  = 2   ! Y_{mu^-}
  integer, parameter :: i_lep_ymu_plus   = 3   ! Y_{mu^+}
  integer, parameter :: i_lep_press_muminus = 4 ! P_{mu^-}
  integer, parameter :: i_lep_press_muplus  = 5 ! P_{mu^+}
  integer, parameter :: i_lep_eps_muminus   = 6 ! eps_{mu^-}
  integer, parameter :: i_lep_eps_muplus    = 7 ! eps_{mu^+}
  integer, parameter :: i_lep_s_muminus     = 8 ! s_{mu^-}
  integer, parameter :: i_lep_s_muplus      = 9 ! s_{mu^+}

  ! -----------------------------------------------------------------------
  ! Allocated data arrays
  !   Layout chosen to match intep3d_many:  ft(nrho, ntemp, n3rd, nvars)
  ! -----------------------------------------------------------------------
  double precision, allocatable, save :: lep_tables_baryon(:,:,:,:)
  double precision, allocatable, save :: lep_tables_ele(:,:,:,:)
  double precision, allocatable, save :: lep_tables_muon(:,:,:,:)

  ! -----------------------------------------------------------------------
  ! Global EOS bounds (set by reader, used by con2prim and EOS wrappers)
  ! -----------------------------------------------------------------------
  double precision :: lep_eos_rhomin  = 0.0d0
  double precision :: lep_eos_rhomax  = 0.0d0
  double precision :: lep_eos_tempmin = 0.0d0
  double precision :: lep_eos_tempmax = 0.0d0
  double precision :: lep_eos_yemin   = 0.0d0   ! min yp  (baryon table 3rd axis)
  double precision :: lep_eos_yemax   = 0.0d0   ! max yp
  double precision :: lep_eos_ylemin  = 0.0d0   ! min yle (ele table 3rd axis)
  double precision :: lep_eos_ylemax  = 0.0d0   ! max yle
  double precision :: lep_eos_ymumin  = 0.0d0   ! min ymu (muon table 3rd axis, linear)
  double precision :: lep_eos_ymumax  = 0.0d0   ! max ymu
  double precision :: lep_eos_epsmin  = 0.0d0
  double precision :: lep_eos_epsmax  = 0.0d0
  double precision :: lep_eos_hmin    = 1.0d0   ! min specific enthalpy h = 1+eps+p/rho

  ! -----------------------------------------------------------------------
  ! Flags
  ! -----------------------------------------------------------------------
  logical :: lep_use_muons            = .true.   !< include muonic contribution
  logical :: lep_add_ele_contribution = .true.   !< include explicit electron leptonic table
  logical :: lep_fix_ymu_high_yp      = .true.   !< set ymu=ymumin when yle+ymu>yemax
  logical :: lep_extend_table_high    = .false.  !< extend beyond table at high rho/T

end module mod_eos_leptonic_parameters
