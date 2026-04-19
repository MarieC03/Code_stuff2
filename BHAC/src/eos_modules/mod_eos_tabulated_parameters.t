!> Module for eos
module mod_eos_tabulated_parameters

  implicit none
  public
  logical                                 :: fil_eps_range_bound = .false.
  character(len=128)                      :: eos_table_name = 'There is no eos table imported' 
  double precision                        :: eos_precision = 1.0d-14
  integer, parameter                      :: eos_iter_max = 20 ! max allowed iteration for the root finding

  
  character(len=128)                      :: table_type_input = 'There is no table type eos imported' 
  integer                                 :: table_type
  integer, parameter                      :: scollapse = 1
  integer, parameter                      :: compose   = 2


  integer                                 :: nrho,ntemp,nye

  integer                                 :: warn_from !warn from given reflevel

  double precision                        :: energy_shift = 0.0d0

  double precision                        :: t_max_hack = 240.0d0

  ! basics
  {#IFDEF M1
    {#IFDEF M1_RATES
    logical                               :: use_realistic_mp_table = .true.
    integer                               :: nvars = 25 ! if use realistic mp table, nvars = 25

    } 
    {#IFNDEF M1_RATES
    logical                                 :: use_realistic_mp_table = .false.
    integer                                 :: nvars = 6 ! if use realistic mp table, nvars = 25
    }
  }
  {#IFNDEF M1
  logical                                 :: use_realistic_mp_table = .false.  ! use realistic microphysics or not to change eos_mp_nvars to 25 (without quark) or 6
  integer                                 :: nvars = 6 ! if use realistic mp table, nvars = 25
  }
  integer                                 :: nthermo, ncomp, nav
  double precision, allocatable, save     :: eos_tables(:,:,:,:)
  integer, parameter                      :: i_logpress = 1
  integer, parameter                      :: i_logenergy = 2
  integer, parameter                      :: i_cs2 = 3
  integer, parameter                      :: i_mu_e = 4
  integer, parameter                      :: i_mu_p = 5
  integer, parameter                      :: i_mu_n = 6
  integer, parameter                      :: i_munu = 7
  integer, parameter                      :: i_entropy = 8
  integer, parameter                      :: i_dedT = 9
  integer, parameter                      :: i_dpdrhoe = 10
  integer, parameter                      :: i_dpderho = 11
  integer, parameter                      :: i_muhat = 12
  integer, parameter                      :: i_xa = 13
  integer, parameter                      :: i_xh = 14
  integer, parameter                      :: i_xn = 15
  integer, parameter                      :: i_xp = 16
  integer, parameter                      :: i_abar = 17
  integer, parameter                      :: i_zbar = 18
  integer, parameter                      :: i_gamma = 19
  ! for post-processing, tables can also have kano, hyperons, and other particles
  ! By modifying here and look for the indices of those particles, you can output their compositions
  integer, parameter                      :: i_quark_u = 20
  integer, parameter                      :: i_quark_d = 21
  integer, parameter                      :: i_quark_c = 22
  integer, parameter                      :: i_quark_s = 23
  integer, parameter                      :: i_quark_t = 24
  integer, parameter                      :: i_quark_b = 25

  double precision, allocatable :: logrho_table(:)
  double precision, allocatable :: logtemp_table(:)
  double precision, allocatable :: ye_table(:)

  ! integer that stores if the speed of
  ! sound in the table has already been
  ! divided by the specific enthalpy
  integer :: have_rel_cs2

!  logical :: Inversion_T_to_eps = .true.   
end module mod_eos_tabulated_parameters
