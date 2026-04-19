module mod_m1_internal
  implicit none

  public
  
  !!integer, parameter :: m1_numvars_internal = 1+ ^NC + 1 + 1 + 1+ ^NC ! E,Fi,Gamma,J,zeta,Hi
  integer, parameter :: m1_numvars_internal = 1 + ^NC  ! only for E,Fi
  integer, parameter :: m1_energy_ = 1 
  {^C& integer, parameter :: m1_flux^C_ = m1_energy_+^C \}
  
  !integer, parameter :: m1_system_ndim = m1_numvars_internal
  ! old version:
  !integer, parameter :: m1_Gammadens_ = m1_flux^NC_ + 1
  !integer, parameter :: m1_Jdens_ = m1_ndens_ + 1
  !integer, parameter :: m1_zeta_ = m1_Jdens_ + 1
  !{^C& integer, parameter :: m1_Hlow_ = m1_Jdens_+^C \}
  
  
end module mod_m1_internal
