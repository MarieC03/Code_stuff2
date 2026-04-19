module mod_m1_tests

  {#IFDEF UNIT_TESTS 
  integer, parameter  :: m1_npress = ((^NC)**2-(^NC))/2 + (^NC)
  integer, parameter  :: m1_nvars = 1 + 1 + (^NC)
  !----------------------------------------------------
  !------ Unit test m1_update_closure_ixD -------------
  ! integer, parameter :: nw   = (3+(^NC))*(^NS) + (^NC)
  ! integer, parameter :: ndim = (^ND)
  !----------------------------------------------------
  !------ Unit test m1_update_closure -----------------
  !integer, parameter :: Nspecies1 = 1
  !integer, parameter :: nw = (m1_nvars+Nspecies1*m1_nvars) + Nspecies1 + 1+ (^NC) +1
  !integer, parameter :: ndim = (^ND)
  !integer, parameter :: ndir = (^NC)
  ! starting indicees for radiation variables in wprim()
  ! previous:
  !integer, parameter :: nrad1_ = 1    
  !integer, parameter :: erad1_ = nrad1_ + 1
  !integer, parameter :: frad1^C_ = erad1_ + ^C  
  !integer, parameter :: zeta1_ = m1_nvars + Nspecies1* m1_nvars   
  !integer, parameter :: u0_ = (zeta1_ + Nspecies1) +1
  !-------------------------------------------------------
  !------ Unit test m1_add_geometrical_source ------------
  !-------------------------------------------------------
  ! add the following to test m1_update_closure parameters
  !integer, parameter :: species2 = 2
  !integer, parameter :: nrad2_ = nrad1_ + species2 * m1_nvars 
  !integer, parameter :: erad2_ = erad1_ + species2 * m1_nvars 
  !integer, parameter :: frad2^C_ = frad1^C_ + species2 * m1_nvars 
  !integer, parameter :: species3 = 3
  !integer, parameter :: nrad3_ = nrad1_ + species3 * m1_nvars 
  !integer, parameter :: erad3_ = erad1_ + species3 * m1_nvars 
  !integer, parameter :: frad3^C_ = frad1^C_ + species3 * m1_nvars 
  !integer, parameter :: species4 = 4
  !integer, parameter :: nrad4_ = nrad1_ + species4 * m1_nvars 
  !integer, parameter :: erad4_ = erad1_ + species4 * m1_nvars 
  !integer, parameter :: frad4^C_ = frad1^C_ + species4 * m1_nvars  
  !-------------------------------------------------------
  !------ Unit test m1_get_fluxes ------------
  !----------------------------------------------------
  integer, parameter :: Nspecies1 = (^NS)
  integer, parameter :: ndim = (^ND)
  integer, parameter :: ndir = (^NC)
  ! explicit vars
  integer, parameter :: nrad1_ = 2 ! arbitrary number here   
  integer, parameter :: erad1_ = nrad1_ + 1
  {^C& integer, parameter :: frad1^C_ = erad1_ + ^C  \}
  !integer, parameter :: zeta1_ = nrad1_ + Nspecies1* m1_nvars   
  integer, parameter :: u0_ = (nrad1_-1) + (Nspecies1*m1_nvars) + 1 !(zeta1_ + Nspecies1)
  
  integer, parameter :: species2 = 2
  integer, parameter :: nrad2_ = nrad1_ + (species2-1) * m1_nvars 
  integer, parameter :: erad2_ = erad1_ + (species2-1) * m1_nvars 
  {^C& integer, parameter :: frad2^C_ = frad1^C_ + (species2-1) * m1_nvars \} 
  integer, parameter :: species3 = 3
  integer, parameter :: nrad3_ = nrad1_ + (species3-1) * m1_nvars 
  integer, parameter :: erad3_ = erad1_ + (species3-1) * m1_nvars 
  {^C& integer, parameter :: frad3^C_ = frad1^C_ + (species3-1) * m1_nvars \} 
  integer, parameter :: species4 = 4
  integer, parameter :: nrad4_ = nrad1_ + (species4-1) * m1_nvars 
  integer, parameter :: erad4_ = erad1_ + (species4-1) * m1_nvars 
  {^C& integer, parameter :: frad4^C_ = frad1^C_ + (species4-1) * m1_nvars \} 

  ! nrad,erad,frad*3 for each species,u0 
  !integer, parameter :: nwrad = (m1_nvars+Nspecies1*m1_nvars) + Nspecies1+ 1+ (^NC) +1 ! +1 for safety
  integer, parameter :: nwrad = (nrad1_-1) + (Nspecies1*m1_nvars) + 1+ (^NC) +1 !+Nspecies1 for zeta
  integer, parameter :: ncons = 2*nwrad     ! U, S
  integer, parameter :: nwflux = ncons+nwrad       ! U, S, F
  integer, parameter :: nw = nwflux ! for now

  !-------------------------
  ! -----implicit vars
  !---------------
  ! option 1 : adding an array with implicit vars outside of advect (& pass it to finite volume as optional)
  !---------------
  integer, parameter :: m1_nvars_impl = 5 !kappa_as, Q, Q_n,kappa_n
  integer, parameter :: kappa_a1_ = 1
  integer, parameter :: kappa_s1_ = kappa_a1_ + 1
  integer, parameter :: Q_er1_ = kappa_s1_ + 1
  integer, parameter :: Q_nr1_ = Q_er1_ + 1
  integer, parameter :: kappa_nr1_ = Q_nr1_ + 1
  
  integer, parameter :: nwradimpl = (kappa_a1_-1)+(Nspecies1*m1_nvars_impl)

  ! other species:
  integer, parameter :: kappa_a2_ = kappa_a1_ + (species2-1) * m1_nvars_impl
  integer, parameter :: kappa_s2_ = kappa_s1_ + (species2-1) * m1_nvars_impl
  integer, parameter :: Q_er2_ = Q_er1_ + (species2-1) * m1_nvars_impl
  integer, parameter :: Q_nr2_ = Q_nr1_ + (species2-1) * m1_nvars_impl
  integer, parameter :: kappa_nr2_ = kappa_nr1_ + (species2-1) * m1_nvars_impl

  integer, parameter :: kappa_a3_ = kappa_a1_ + (species3-1) * m1_nvars_impl
  integer, parameter :: kappa_s3_ = kappa_s1_ + (species3-1) * m1_nvars_impl
  integer, parameter :: Q_er3_ = Q_er1_ + (species3-1) * m1_nvars_impl
  integer, parameter :: Q_nr3_ = Q_nr1_ + (species3-1) * m1_nvars_impl
  integer, parameter :: kappa_nr3_ = kappa_nr1_ + (species3-1) * m1_nvars_impl

  integer, parameter :: kappa_a4_ = kappa_a1_ + (species4-1) * m1_nvars_impl
  integer, parameter :: kappa_s4_ = kappa_s1_ + (species4-1) * m1_nvars_impl
  integer, parameter :: Q_er4_ = Q_er1_ + (species4-1) * m1_nvars_impl
  integer, parameter :: Q_nr4_ = Q_nr1_ + (species4-1) * m1_nvars_impl
  integer, parameter :: kappa_nr4_ = kappa_nr1_ + (species4-1) * m1_nvars_impl
  !integer, parameter :: m1_nvars_impl = 5 !kappa_as, Q_eas
  integer, parameter :: m1_nvars_tot = m1_nvars + m1_nvars_impl
  !...
  !--------------------------------
  integer, parameter :: tau_ = 1
  {^C& integer, parameter :: s^C_ = tau_ + ^C \}
  integer, parameter :: ye_ = tau_ + ^NC + 1
  ! need T_eps_ and rho_ before we do wfluidprim cons to prim
  integer, parameter :: rho_ = ye_+1
  integer, parameter :: T_eps_ = rho_ + 1 

  } !end IFDEF UNIT_TESTS
  
end module mod_m1_tests

