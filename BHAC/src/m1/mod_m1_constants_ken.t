module mod_m1_constants

  implicit none
  public

  ! Geometric Units: c = G = M_sun = 1
  ! These conversion factors transform physical CGS units to computational units (CU)
  
  ! ============================================================================
  ! FUNDAMENTAL UNIT CONVERSIONS (CGS to Geometric Units)
  ! ============================================================================
  
  ! Length: characteristic scale is GM_sun/c^2 ~ 1.477 km
  ! 1 cm in CGS = LENGTHGF in CU
  double precision, parameter :: LENGTHGF = 6.77269222552442d-6   ! cm to CU
  
  ! Time: characteristic scale is GM_sun/c^3 ~ 4.926 microseconds
  ! 1 second in CGS = TIMEGF in CU
  double precision, parameter :: TIMEGF   = 2.03040204956746d5    ! s to CU
  
  ! Density: characteristic scale is M_sun/(GM_sun/c^2)^3 ~ 6.177e17 g/cm^3
  ! 1 g/cm^3 in CGS = RHOGF in CU
  double precision, parameter :: RHOGF    = 1.61887093132742d-18  ! g/cm^3 to CU
  
  ! Pressure: characteristic scale is c^4/(G^2*M_sun) ~ 5.55e38 dyn/cm^2
  ! 1 dyn/cm^2 in CGS = PRESSGF in CU
  double precision, parameter :: PRESSGF  = 1.80123683248503d-39  ! dyn/cm^2 to CU
  
  ! Specific energy (energy per unit mass): 1/c^2 in CGS
  ! Converts from erg/g to dimensionless (since c=1 in geometric units)
  double precision, parameter :: EPSGF    = 1.11265005605362d-21  ! erg/g to CU (= 1/c^2)
  
  ! ============================================================================
  ! INVERSE CONVERSIONS (CU to CGS)
  ! ============================================================================
  
  double precision, parameter :: INVRHOGF   = 6.1771447040563d+17   ! CU to g/cm^3
  double precision, parameter :: INVEPSGF   = 8.98755178736818d20   ! CU to erg/g
  double precision, parameter :: INVPRESSGF = 5.55174079257738d38   ! CU to dyn/cm^2
  
  ! ============================================================================
  ! ENERGY UNIT CONVERSIONS
  ! ============================================================================
  
  ! Nuclear physics energy conversions
  double precision, parameter :: MEV_TO_ERG = 1.60217733d-6       ! MeV to erg
  double precision, parameter :: ERG_TO_MEV = 6.24150636d+5       ! erg to MeV
  
  ! ============================================================================
  ! FUNDAMENTAL PHYSICAL CONSTANTS
  ! ============================================================================
  
  ! Speed of light and related
  double precision, parameter :: C2_CGS        = 8.9875517873681764d+20  ! c^2 in CGS (cm^2/s^2)
  double precision, parameter :: CLITE_CM      = 2.99792458d+10          ! speed of light (cm/s)
  
  ! Nuclear mass constants
  double precision, parameter :: MNUC_MEV      = 931.494061              ! atomic mass unit (MeV/c^2)
  double precision, parameter :: MNUC_CGS      = 1.660539040d-24         ! atomic mass unit (g)
  double precision, parameter :: MNUC_MSUN      = 8.351131764232556E-58  ! CU
  double precision, parameter :: MNUC_CU       = MNUC_CGS * &
     PRESSGF*(LENGTHGF*LENGTHGF*LENGTHGF)*C2_CGS                         ! atomic mass unit in CU
  double precision, parameter :: M_NEUTRON_MEV = 939.565379              ! neutron mass (MeV/c^2)
  
  ! Avogadros number and volume conversion
  double precision, parameter :: AVOGADRO      = 6.0221367d23            ! particles per mole
  double precision, parameter :: CM3_TO_FM3    = 1.0d39                  ! cm^3 to fm^3 conversion
  
  ! ============================================================================
  ! WEAK INTERACTION PARAMETERS
  ! ============================================================================
  
  ! Neutrino cross section and weak interaction parameters
  double precision, parameter :: SIGMA_0       = 1.76d-44                ! reference cross section (cm^2)
  double precision, parameter :: ALPHA         = 1.25d0                  ! cross section power law index
  double precision, parameter :: QNP           = 1.293333                ! neutron-proton mass difference (MeV)
  double precision, parameter :: HC_MEVCM      = 1.23984172d-10          ! hbar*c (MeV*cm)
  double precision, parameter :: CV            = 0.5d0 + 2.0d0 * 0.23d0  ! vector coupling constant
  double precision, parameter :: CA            = 0.5d0                   ! axial coupling constant
  double precision, parameter :: GAMMA_0       = 5.565d-2                ! reference rate (MeV)
  
  ! ============================================================================
  ! MATHEMATICAL CONSTANTS
  ! ============================================================================
  
  double precision, parameter :: PI            = 3.14159265358979323846  ! pi
  double precision, parameter :: kB_MEV        = 8.61738568d-11          ! Boltzmann constant (MeV/K)



  ! ============================================================================
! UNITS IN THE CODE - KEY PHYSICAL QUANTITIES
! ============================================================================
!
! This code uses geometric units where c = G = M_sun = 1
! Below are the units for key physical quantities in computational units (CU)
!
! ------------------------------------------------------------------------
! NEUTRINO QUANTITIES
! ------------------------------------------------------------------------
!
! Energy Density (e_nu):
!   - Computational Units: [CU] = [M_sun / (GM_sun/c^2)^3]
!   - Physical Units: Use INVPRESSGF to convert to erg/cm^3
!   - Example: e_nu_cgs = e_nu_cu * INVPRESSGF * C2_CGS
!   - Note: In GR, energy density has same dimensions as pressure
!
! Number Density (n_nu):
!   - Computational Units: [CU] = [1 / (GM_sun/c^2)^3]
!   - Physical Units: Use INVRHOGF/MNUC_CGS to convert to particles/cm^3
!   - Example: n_nu_cgs = n_nu_cu * (INVRHOGF / MNUC_CGS)
!
! ------------------------------------------------------------------------
! THERMODYNAMIC QUANTITIES
! ------------------------------------------------------------------------
!
! Electron Fraction (Y_e):
!   - Computational Units: Dimensionless [CU] = dimensionless
!   - Physical Units: Dimensionless (number of electrons per baryon)
!   - Definition: Y_e = n_e / n_b = (protons) / (protons + neutrons)
!   - Range: 0 < Y_e < 1 (typically 0.1-0.5 in neutron star matter)
!   - Note: Y_e is invariant across unit systems
!
! Matter Temperature (T):
!   - Computational Units: [CU] = [MeV]
!   - Physical Units: MeV (same as neutrino temperature)
!   - To Kelvin: T_K = T_MeV / kB_MEV
!   - Typical range: 1-100 MeV in core-collapse supernovae
!
! Matter Density (rho):
!   - Computational Units: [CU] = [M_sun / (GM_sun/c^2)^3]
!   - Physical Units: Use INVRHOGF to convert to g/cm^3
!   - Example: rho_cgs = rho_cu * INVRHOGF
!   - Nuclear density: rho_nuc ~ 2.8e14 g/cm^3 ~ 0.16 fm^-3
!
! Specific Internal Energy (eps):
!   - Computational Units: [CU] = dimensionless (energy per mass with c=1)
!   - Physical Units: Use INVEPSGF to convert to erg/g
!   - Example: eps_cgs = eps_cu * INVEPSGF
!   - Note: In geometric units, eps = E/(mc^2) is dimensionless
!
! ------------------------------------------------------------------------
! CONVERSION EXAMPLES
! ------------------------------------------------------------------------
!
! Example 1: Convert neutrino energy density to physical units
!   e_nu_erg_per_cm3 = e_nu_cu * INVPRESSGF * C2_CGS
!
! Example 2: Convert matter density and Y_e to electron number density
!   n_e_per_cm3 = (rho_cu * INVRHOGF / MNUC_CGS) * Y_e
!
! Example 3: Convert temperature from MeV to Kelvin
!   T_kelvin = T_MeV / kB_MEV  ! kB_MEV = 8.617e-11 MeV/K
!
! Example 4: Compute baryon number density
!   n_b_per_cm3 = rho_cu * INVRHOGF / MNUC_CGS
!   n_b_per_fm3 = n_b_per_cm3 / CM3_TO_FM3
!
! ============================================================================

end module mod_m1_constants



