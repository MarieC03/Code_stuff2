module mod_m1_eas_param
  implicit none

  public

  {#IFDEF M1_TAU_OPT
  integer, parameter :: m1_num_eas = 6 
  }
  {#IFNDEF M1_TAU_OPT
  integer, parameter :: m1_num_eas = 5
  }
  integer, parameter :: Q_ems   = 1
  integer, parameter :: k_a   = 2
  integer, parameter :: k_s   = 3
  integer, parameter :: Q_ems_n = 4
  integer, parameter :: k_n   = 5
  integer, parameter :: tau_path   = 6


  integer, parameter :: fluid_vars = 3
  integer, parameter :: idx_rho = 1
  integer, parameter :: idx_T   = 2
  integer, parameter :: idx_Ye  = 3

  integer :: m1_i_nue = 1
  integer :: m1_i_nuebar = 2
  integer :: m1_i_nux = 3
  integer :: m1_i_mu = 4
  integer :: m1_i_mubar = 5
  integer :: m1_i_photon = 6

  double precision, allocatable :: m1_eas_ixD(:, :) !(m1_num_eas,^NS)  !(:,:,:) 
  ! not used:
  double precision, allocatable :: m1_eas({^D&:,}, :, :) !(ixI^L,m1_num_eas,^NS)  !(:,:,:)

end module mod_m1_eas_param
