module mod_weakhub_parameters
  !use mod_lookup_table, only: LT4_t, LT3_t, LT2_t

  implicit none
  public
    !.............................................
    ! The precision of the real numbers used in the tables
    integer, parameter :: dp = kind(1.0d0)

    !> weakhub method: 1. multigroup + otf (default)
    !>                 2. multigroup + hybrid of otf + tabulated
    !>                 3. multigroup + full tabulated
    !> weakhub_grey method
    !>                 1. single group, integrated all opacities and kernels
    !>                 2. single group, analytical approx.
    integer                             :: Weakhub_method          = -1
    integer                             :: Weakhub_grey_method     = 1 !-1
    logical                             :: Create_weakhub_table    = .false.
    logical                             :: use_emissivity_table    = .false.
  type weakhub_t
    character(512)                      :: table_path              = ""
    character(512)                      :: energy_bin_spacing      = "none"    ! log or uni
    logical                             :: energy_bin_centroid     = .true.
    double precision                    :: eps_max           ! = rad_eps_max = rad_epsi(n_ebin)
    double precision                    :: eps_min           ! = rad_eps_min = rad_epsi(0)
    double precision, allocatable       :: energies(:)       ! size: 0:weakhub_table%n_groups+1
    double precision, allocatable       :: energies_i(:)     ! size: 0:weakhub_table%n_groups+1
    double precision, allocatable       :: dV_eps(:)         ! size: 1:weakhub_table%n_groups
    integer                             :: n_species         = -1
    integer                             :: n_groups          = -1
  
  end type weakhub_t

  type(weakhub_t) :: weakhub_table

  ! Table uses: 
    integer :: IVrho    = -1     ! in 4D Weakhub table, the number of points of rho 
    integer :: IVtemp   = -1     ! in 4D Weakhub table, the number of points of temp
    integer :: IVye     = -1     ! in 4D Weakhub table, the number of points of ye 
    integer :: IVymu    = -1     ! in 4D Weakhub table, the number of points of ymu
    integer :: IIIrho   = -1     ! in 3D Weakhub table, the number of points of rho 
    integer :: IIItemp  = -1     ! in 3D Weakhub table, the number of points of temp
    integer :: IIIyp    = -1     ! in 3D Weakhub table, the number of points of yp, Warning this is proton fraction!
    integer :: IIeta_e  = -1     ! in 2D Weakhub table, the number of points of eta_e
    integer :: IItemp   = -1     ! in 2D Weakhub table, the number of points of temp
    ! 4D table bounds
    double precision :: logrho_max_IV
    double precision :: logrho_min_IV
    double precision :: logtemp_max_IV
    double precision :: logtemp_min_IV
    double precision :: ye_max_IV
    double precision :: ye_min_IV
    double precision :: logymu_max_IV
    double precision :: logymu_min_IV

    ! 3D table bounds
    double precision :: logrho_max_III
    double precision :: logrho_min_III
    double precision :: logtemp_max_III
    double precision :: logtemp_min_III
    double precision :: yp_max_III
    double precision :: yp_min_III

    ! 2D table bounds
    double precision :: logeta_e_max_II    ! ln(eta_e_max)  for Kernels, mue included rest mass
    double precision :: logeta_e_min_II    ! ln(eta_e_min)  for Kernels
    double precision :: logtemp_max_II     ! ln(temp_max)   for Kernels
    double precision :: logtemp_min_II     ! ln(temp_min)   for Kernels   !as it > 0.01MeV need to limit the radius/nu_rho_min

    
    integer :: weakhub_table_type = -1 
    ! Does not exist only has 3D table
    integer :: table_allD = 1    ! Including 2D, 3D, 4D tables
    integer :: table_IVandII = 2    ! Including 2D, 4D only
    integer :: table_IV  = 3     ! only has 4D table
    integer :: table_II  = 4  ! only has 2D table
    integer :: table_III  = 5  ! only has 3D table
    integer :: table_grey  = 100 ! only has 4D opacity table and grey

    ! Opacity tables
    double precision :: Otable_rho_LB = 1.0d-15   ! Lower bound of rho for Otable
    double precision :: Otable_logtemp_LB = dlog(0.01d0) ! Lower bound of rho for Otable
    ! opacity as a function of (rho, T, ye, ymu), due to the memory problems, we can only do this fixme I need to add n_ebin
    ! They are as function of rho, T, ye, ymu becoz of medium modifications Un, mn_eff
    ! Virtual memory insufficient for assigning nrho*ntemp*nye*nymu, n_ebin as well! you need to use more memory when generate the code
    logical :: use_IVtable = .false.  ! use 4D table for opacity
    logical :: use_IIItable = .false. ! use 3D table for kernels
    logical :: use_IItable = .false.  ! use 2D table for kernels

    ! ----- Methods flags :
    logical, save                       :: use_emissivity_NNbrem          = .false.
    logical, save                       :: use_emissivity_epann           = .false.

    logical, save                       :: Sion_nuclei_scat               = .true.
    integer, save                       :: H_scat_method
    integer, parameter                  :: ILEAS_kappa_s_nuclei           = 1
    integer, parameter                  :: Burrow_kappa_s_nuclei          = 2
    logical, save                       :: gauss_quadrature_Iscat         = .true.
    logical, save                       :: Simple_coherent_Freedman_scat  = .false. 
    ! Corrections :
    logical, save                       :: weak_magnetism_approx          = .false.
    logical, save                       :: N_recoil_approx                = .false.
    logical, save                       :: NC_virial_and_RPA              = .false.
    logical, save                       :: Strangequark_coupling          = .false.
    logical, save                       :: Thompson_NNbrem_factor         = .true.
  
    ! Full kinematics beta processes
    logical, save                       :: beta_proc_full                 = .false.
    logical, save                       :: weak_magnetism_full            = .false.
    logical, save                       :: pseudoscalar                   = .false.
    logical, save                       :: nucleon_formfactor             = .false.
  
  !  logical, save                       :: table_inelastic_Processes      = .false.  ! use table format for IS
  !  logical, save                       :: table_pair_Processes           = .false.  ! use table format for pair processes
  

    ! monitor the requirements of corresponding types of kernels for different implicit solvers in weakhub (save computational cost)
    ! default true
    logical, save                       :: present_IS                     = .true.
    logical, save                       :: present_PP                     = .true.
    logical, save                       :: present_nunuPP                 = .true.

    ! Variables that physically related: 

    ! suppression factor for opacity/kernels for muon/anti-muon creation processes
    logical                             :: mu_suppression_factor = .false.
    ! activate the conservation law for nunu annihilation
    logical                             :: nunubar_pa_conserve = .True.
    ! density threshold for nunu annihilation to be non-zero (Density to assume nue nuebar in LTE with matter)
    double precision                    :: rho_nunupa = 10.0d0**(12.5d0)

    ! Weak correction factors
    double precision, allocatable       :: W_m_CC(:,:)
    double precision, allocatable       :: W_m_NC(:,:,:)
    double precision, allocatable       :: W_recoil_CC(:,:)
    double precision, allocatable       :: W_recoil_NC(:,:,:)

  !---------- Neutrino-related variables and their indices (can be used in leakage/weakhub) ---------!

  ! temerature in MeV
  integer  :: temp

  ! degeneracy parameters including rest mass, 1: electrons  2: p  3: n  4: mu  5: tau fixme: Leptons not nucleons, pls fix the names
  integer, allocatable  :: eta_nucleons(:)
  integer, parameter    :: N_tableparticles = 5

  ! n degeneracy - p dengeneracy - Qnp/temp
  integer  :: eta_hat

  ! nu_e degeneracy, from eta_e - eta_n + eta_p not necessarily in weak beta eqm;  1:n_spec
  integer, allocatable  :: eta_nu(:)

  ! lepton blocking terms for e- e+
  integer, allocatable  :: lepton_blocking(:)


  ! mass fractions;1. Xn   2. Xp   3.  Xa   4. XH   5. X2H   6. X3H   7. X3He
  integer, allocatable  :: mass_fraction(:)
  integer, parameter    :: n_Mfrac = 7

  integer  :: abar
  integer  :: zbar

  !nucleon final-state inhibition(blocking) factors for E/A
  integer  :: eta_pn
  integer  :: eta_np

  !nucleon final-state blocking factors for scattering 1. PP  2. NN
  integer, allocatable :: ksi_NN(:)

  ! degeneracy parameter in weak beta eqm; 1:n_spec
  integer, allocatable :: eta_nu_eq(:)

  ! baryonic number denisty w.r.t neutron mass
  integer, allocatable :: n_b

  ! Particle number density; 1. p   2. n   3. alp  4. Heavy nuclei 5. e-  6. 2H  7. 3H  8. 3He
  integer, allocatable :: n_N(:)
  integer, parameter   :: n_particle = 8

  ! Mean field Potential
  integer  :: U_n
  integer  :: U_p
  ! Effective mass of nucleons
  integer  :: mn_eff
  integer  :: mp_eff

  integer  :: S_A
  integer  :: S_V

  integer  :: omegap ! plasma frequency - the lowest limit for the value of plasmon mass
  integer  :: omega1
  integer  :: omegaA

  ! Specialized for leakage
  integer  :: lvol !cell volume in cm^3, indices <radial position>
  integer  :: lIarea !inverse cell area in cm^-2, indices <radial position>
  integer  :: lmass !cell mass in g, indices <radial position>

  integer, allocatable  :: lcoolingsource(:)   ! ndir are mom, ndir + 1 : tau, ndir +2 : Dye

  !number emission rate used in determining change in ye
  integer, allocatable  :: R_eff(:) !effective number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3
  integer, allocatable  :: R_loc(:) !local number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3
  integer, allocatable  :: R_diff(:) !diffusive number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3
  integer, allocatable  :: R_tot(:)

  integer, allocatable  :: Q_eff(:) !effective energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3
  !energy emission rate used in determining change in energy
  integer, allocatable  :: Q_loc(:) !local energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3
  integer, allocatable  :: Q_diff(:) !diffusive energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3
  integer, allocatable  :: Q_tot(:)

  integer, allocatable  :: kappa_a_bar_num(:)
  integer, allocatable  :: kappa_a_bar_en(:)
  integer, allocatable  :: kappa_s_bar_en(:)

  integer, allocatable :: bar_E_num(:)  ! integrated neu spectrum neu energy eq(22)  <neu species: j>
  integer, allocatable :: bar_E_en(:)  ! integrated neu spectrum neu energy eq(22)  <neu species: j>

  ! production timescale
  !integer, allocatable :: time_prod(:,:) !prod time scale    < neu species: j>
  integer, allocatable :: time_prod_num(:) !prod time scale    < neu species: j>
  integer, allocatable :: time_prod_en(:) !prod time scale    < neu species: j>
  ! diffusion timescale
  !integer, allocatable :: time_diff(:,:) !prod time scale    < neu species: j>
  integer, allocatable :: time_diff_num(:) !prod time scale    < neu species: j>
  integer, allocatable :: time_diff_en(:) !prod time scale    < neu species: j>



  integer, allocatable :: kappa_diff_tot_nue(:) !  < n_eps >
  integer, allocatable :: kappa_diff_tot_nua(:) !  < n_eps >
  integer, allocatable :: kappa_diff_tot_nux(:) !  < n_eps >


  integer, allocatable :: gammaeff_num(:)
  integer, allocatable :: gammaeff_en(:)

  integer, allocatable :: Q_abs_loc(:)
  integer, allocatable :: R_abs_loc(:)

  integer, allocatable :: Q_abs_eff(:)
  integer, allocatable :: R_abs_eff(:)

  !-------End of weakvar -------!


  ! for new variables for weakhub interactions kernels and eas
  integer, public    :: nweak = 0

  !> Maximum number of weakhub variables
  integer, parameter :: max_nvar_weak = 500

  !> weakhub variable names
  character(len=32)  :: weak_names(max_nvar_weak)

  contains

  !> Set weakhub var variable
  function var_set_weakhub_var(name_weak_var, ix) result(iw)
    character(len=*), intent(in)  :: name_weak_var !< weakhub var name
    integer, intent(in), optional :: ix        !< Optional index (to make var1, var2, ...)
    integer                       :: iw

    ! total number of weakhub variables
    nweak  = nweak + 1
    iw     = nweak

    if (.not. present(ix)) then
      weak_names(nweak) = name_weak_var
    else
      write(weak_names(nweak),"(A,I0)") name_weak_var, ix
    end if

  end function var_set_weakhub_var

end module mod_weakhub_parameters
