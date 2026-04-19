module mod_m1_constants

  implicit none
  public

  ! unit conversions 
  double precision, parameter :: LENGTHGF = 6.77269222552442d-6   ! cm to CU
  double precision, parameter :: TIMEGF   = 2.03040204956746d5    ! s to CU
  double precision, parameter :: RHOGF    = 1.61887093132742d-18  ! g/cm^3 to CU
  double precision, parameter :: PRESSGF  = 1.80123683248503d-39  ! dyn / cm^2 to CU
  double precision, parameter :: EPSGF    = 1.11265005605362d-21  ! 1. / c_cgs**2
  double precision, parameter :: MEV_TO_ERG = 1.60217733d-6       ! MeV to erg
  double precision, parameter :: MEV_TO_ERG_KEN = 1.60217733d-6       ! MeV to erg
  double precision, parameter :: ERG_TO_MEV = 6.24150636d+5       ! erg to MeV

  ! inverses
  double precision, parameter :: INVRHOGF   = 6.1771447040563d+17
  double precision, parameter :: INVEPSGF   = 8.98755178736818d20
  double precision, parameter :: INVPRESSGF = 5.55174079257738d38

  ! physical constants
  double precision, parameter :: C2_CGS        = 8.9875517873681764d+20
  double precision, parameter :: CLITE_CM      = 2.99792458d+10
  double precision, parameter :: MSUN_GRAMM    = 1.98847E+33
  double precision, parameter :: MNUC_MEV      = 931.494061
  double precision, parameter :: MNUC_CGS      = 1.660539040d-24
  double precision, parameter :: MNUC_MSUN     = 8.351131764232556E-58  ! CU
  double precision, parameter :: MNUC_CU       = MNUC_CGS * PRESSGF*(LENGTHGF*LENGTHGF*LENGTHGF)*C2_CGS
  double precision, parameter :: CM3_TO_FM3    = 1.0d39
  double precision, parameter :: AVOGADRO      = 6.0221367d23
  double precision, parameter :: M_NEUTRON_MEV = 939.565379
  double precision, parameter :: SIGMA_0       = 1.76d-44
  double precision, parameter :: ALPHA         = 1.25d0 !1.23d0
  double precision, parameter :: QNP           = 1.293333
  double precision, parameter :: HC_MEVCM      = 1.23984172d-10
  double precision, parameter :: CV            = 0.5d0 + 2.0d0 * 0.23d0
  double precision, parameter :: CA            = 0.5d0
  double precision, parameter :: GAMMA_0       = 5.565d-2
  double precision, parameter :: PI            = 3.14159265358979323846
  double precision, parameter :: kB_MEV        = 8.61738568d-11
  
end module mod_m1_constants

