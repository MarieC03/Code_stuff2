 !****************************************************************************
 !>
 !> calculate_m1_eas_ixD
 !>
 !> Computes equilibrium emission/absorption/scattering rates for a single
 !> neutrino species at a single grid point, then applies the neutrino-
 !> temperature correction.
 !>
 !> The routine mirrors the two code paths in FIL/Margherita M1.hh calc_eas():
 !>
 !>   #ifdef M1_RATES  (compile flag M1_RATES = Weakhub tables are available)
 !>     Path A: Weakhub table rates + Kirchhoff emission + analytical plasma/brems
 !>   #else
 !>     Path B: Fully analytical rates (charged-current + scattering + plasma)
 !>   #endif
 !>
 !>  In both paths the rates are subsequently corrected by the neutrino
 !>  temperature T_nu following the same three correction modes as before.
 !>
 !> References
 !>   Ruffert et al. 1996, A&A 311, 532   (R96)
 !>   Rosswog & Liebendörfer 2003, MNRAS 342, 673
 !>   Bruenn 1985, ApJS 58, 771
 !****************************************************************************
  subroutine calculate_m1_eas_ixD(wrad,N,ix^D,Gamma,vel,Wlor,Jrad,ixI^L,x,eas, &
       speciesKSP, fluid_Prim, eta_rad, calc_eqRates_here, qdt)
    use mod_m1_internal
    use mod_m1_fermi
    use mod_m1_constants
    use mod_m1_eas_microphysical_gray
    use mod_eos_tabulated
    use mod_eos, only: small_rho, small_temp, big_ye, eos_yemin, eos_yemax
    use mod_Weakhub_reader
    use mod_m1_eas_param
    use, intrinsic :: ieee_arithmetic
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    {#IFNDEF UNIT_TESTS
    !include "amrvacdef.f"
    }
    implicit none

    ! ---- arguments --------------------------------------------------------
    integer,           intent(in)    :: ix^D, ixI^L, speciesKSP
    double precision,  intent(inout) :: wrad(1:m1_numvars_internal)
    double precision,  intent(inout) :: fluid_Prim(1:fluid_vars)
    double precision,  intent(inout) :: eas(1:m1_num_eas)
    double precision,  intent(in)    :: x(ixI^S, 1:ndim)
    double precision,  intent(in)    :: Gamma, Wlor, Jrad
    double precision,  intent(in)    :: vel(1:^NC)
    double precision,  intent(inout) :: N
    double precision,  intent(inout) :: eta_rad   !> neutrino degeneracy parameter
    logical,           intent(in)    :: calc_eqRates_here
    double precision,  intent(in)    :: qdt

    ! ---- local scalars ----------------------------------------------------
    double precision, parameter :: Qnp1 = 1.2935d0   ! MeV, n-p rest-mass diff.
    double precision, parameter :: me_mev = 0.510998910d0

    !> Primitive fluid variables
    double precision :: rho_fluid, T_fluid, ye_fluid, rho_cgs

    !> Chemical potentials and composition (from EOS)
    double precision :: eps, prs, ent, cs2, dedt, dpderho, dpdrhoe
    double precision :: xa, xh, xn, xp, abar1, zbar1
    double precision :: mu_e, mu_n, mu_p, muhat, munu

    !> Neutrino thermodynamics
    double precision :: chem_pot          ! neutrino chemical potential [MeV]
    double precision :: eta_e             ! electron degeneracy parameter
    double precision :: tau_opt           ! optical depth
    double precision :: T_nu              ! neutrino temperature [MeV, code units]
    double precision :: av_energy         ! average neutrino energy

    !> Blackbody functions  (code units)
    double precision :: Black_er, Black_nr

    !> Weakhub equilibrium rates (from table, before any correction)
    double precision :: eas_eq(1:m1_num_eas)

    !> Neutrino-nucleon blocking factors (Ruffert A12-A13 / Rosswog A9)
    double precision :: eta_np, eta_pn    ! cm^{-3}
    double precision :: nb                ! baryon number density  [cm^{-3}]
    double precision :: eta_hat           ! (mu_n - mu_p - Qnp1)/T

    !> Intermediate rate accumulators (all in CGS / MeV before conversion)
    ! Emissivities [MeV s^{-1} cm^{-3}]
    double precision :: Qems, Qems_n          ! per-process
    double precision :: Qems_sec, Qems_n_sec  ! secondary (plasma/brems/pair)
    ! Opacities  [cm^{-1}]
    double precision :: kappa_a_cc, kappa_n_cc   ! charged-current
    double precision :: kappa_s_nc               ! neutral-current scattering
    double precision :: kappa_a_K,  kappa_n_K    ! Kirchhoff-derived

    !> Blocking / Fermi-integral helpers
    double precision :: block_nue, block_nuebar
    double precision :: zeta_nue, zeta_nuebar     ! absorption prefactor ζ
    double precision :: nfac, efac                ! Fermi-ratio factors
    double precision :: efac_s                    ! scattering energy factor
    double precision :: ymp, ymn                  ! Pauli-blocked mass fractions
    double precision :: C_nucleus                 ! coherent nuclear scattering coeff.
    double precision :: fermi_fac, fermi_block    ! helpers for Fermi integrals
    double precision :: energy_av                 ! average energy helper

    !> Kirchhoff Blackbody in MeV [for internal detailed-balance only]
    double precision :: BB_e_mev, BB_n_mev        ! MeV s^{-1} cm^{-3}
    double precision :: g_factor                  ! species multiplicity

    !> T_nu correction bookkeeping
    double precision :: Tnu_ratio, tanh_factor

    !> Logical correction-mode flags (same semantics as the original code)
    logical :: eas_correction_standard
    logical :: eas_correction_Tnu
    logical :: eas_correction_tanh

    !> misc
    double precision :: UNITS_rho
    logical :: is_known_tau_opt

    ! ======================================================================
    !  Correction-mode selection  (same defaults as the original code)
    ! ======================================================================
    eas_correction_standard = .true.
    eas_correction_Tnu      = .false.
    eas_correction_tanh     = .false.

    ! ======================================================================
    !  Unit-conversion factors
    ! ======================================================================
    UNITS_rho = INVRHOGF

    ! ======================================================================
    !  Copy fluid primitives
    ! ======================================================================
    rho_fluid = fluid_Prim(idx_rho)
    ye_fluid  = fluid_Prim(idx_ye)
    T_fluid   = fluid_Prim(idx_T)

    ! ======================================================================
    !  EOS call: chemical potentials and composition
    ! ======================================================================
    {#IFNDEF UNIT_TESTS
    call tabulated_temp_get_all_one_grid(rho_fluid, T_fluid, ye_fluid, &
         eps, prs=prs, ent=ent, cs2=cs2, dedt=dedt,                   &
         dpderho=dpderho, dpdrhoe=dpdrhoe,                             &
         xa=xa, xh=xh, xn=xn, xp=xp,                                  &
         abar=abar1, zbar=zbar1,                                        &
         mu_e=mu_e, mu_n=mu_n, mu_p=mu_p, muhat=muhat, munu=munu)
    }

    ! ======================================================================
    !  Neutrino chemical potential and equilibrium degeneracy parameter
    ! ======================================================================
    select case(speciesKSP)
      case(m1_i_nue)
        chem_pot = mu_e + mu_p - mu_n - Qnp1
      case(m1_i_nuebar)
        chem_pot = -(mu_e + mu_p - mu_n - Qnp1)
      case default        ! nux, numu, numubar
        chem_pot = 0.0d0
    end select

    eta_rad = chem_pot / T_fluid
    eta_e   = mu_e     / T_fluid

    ! ======================================================================
    !  Optical-depth based damping of eta_rad (same as original code)
    ! ======================================================================
    is_known_tau_opt = .true.
    if (is_known_tau_opt) then
      rho_cgs  = rho_fluid * UNITS_rho
      tau_opt  = exp(DLOG10(10.d0) * (0.96d0 * (DLOG10(rho_cgs)/DLOG10(10.0d0) - 11.7d0)))
      {#IFDEF M1_TAU_OPT_USE
       tau_opt = eas(tau_path)
      }
      if (speciesKSP == m1_i_nue .or. speciesKSP == m1_i_nuebar) then
        eta_rad = eta_rad * (1.0d0 - exp(-tau_opt))
      end if
    end if

    ! ======================================================================
    !  Neutrino temperature and average energy
    ! ======================================================================
    call get_neutrino_temp_ixD(wrad, N, ix^D, Gamma, vel, Wlor, Jrad, &
         eta_rad, speciesKSP, fluid_Prim, av_energy, T_nu, T_fluid,   &
         timeCurrent, x(ix^D,:))

    ! ======================================================================
    !  Blackbody functions (code units, used for Q = BB * kappa_a)
    ! ======================================================================
    call blackbody_ixD(wrad, ix^D, speciesKSP, fluid_Prim, T_fluid, &
         eta_rad, Black_er, Black_nr)

    ! ======================================================================
    !  Initialise eas_eq to safe floor values
    ! ======================================================================
    eas_eq(Q_ems)   = 1.0d-30
    eas_eq(Q_ems_n) = 1.0d-30
    eas_eq(k_a)     = 1.0d-30
    eas_eq(k_s)     = 1.0d-30
    eas_eq(k_n)     = 1.0d-30

    ! ======================================================================
    !  Helper: nucleon-blocking factors η_np, η_pn  (Ruffert A12-A13)
    !  ρ < 2×10^{11}: Bruenn 1985 (3.1) non-degenerate limit.
    ! ======================================================================
    ! rho_cgs already set above (tau_opt block); reuse it here.
    nb       = rho_cgs * AVOGADRO         ! baryon number density [cm^{-3}]
    eta_hat  = (mu_n - mu_p - Qnp1) / T_fluid    ! η_hat = (μ_n - μ_p - Q)/T

    if (rho_cgs < 2.0d11) then
      ! Bruenn non-degenerate limit
      eta_pn = nb * xp
      eta_np = nb * xn
    else
      ! Degenerate limit (Ruffert A12-A13)
      if (abs(eta_hat) > 1.0d-10) then
        eta_pn = nb * (xn - xp) / (exp(eta_hat)  - 1.0d0)
        eta_np = nb * (xp - xn) / (exp(-eta_hat) - 1.0d0)
      else
        ! Limit as η_hat → 0
        eta_pn = nb * xp
        eta_np = nb * xn
      end if
    end if
    ! Rosswog (A9): positivity
    eta_pn = max(eta_pn, 0.0d0)
    eta_np = max(eta_np, 0.0d0)
    if (.not. ieee_is_finite(eta_pn)) eta_pn = 0.0d0
    if (.not. ieee_is_finite(eta_np)) eta_np = 0.0d0

    ! ======================================================================
    !  Species multiplicity factor for Blackbody in MeV
    ! ======================================================================
    if (speciesKSP == m1_i_nue .or. speciesKSP == m1_i_nuebar) then
      g_factor = 1.0d0
    else if (m1_use_muons) then
      g_factor = 1.0d0
    else
      g_factor = 4.0d0    ! nux = νμ + ν̄μ + ντ + ν̄τ in 3-species scheme
    end if

    ! Blackbody number/energy density in MeV (for internal Kirchhoff use)
    ! BB_n = g * 4π c / (ℏc)^3 * T^3 * F_2(η)   [MeV^{-1} s^{-1} cm^{-3}]  ... actually 1/s/cm^3
    ! BB_e = g * 4π c / (ℏc)^3 * T^4 * F_3(η)   [MeV s^{-1} cm^{-3}]
    BB_n_mev = g_factor * 4.0d0*PI * CLITE_CM / HC_MEVCM**3 * &
               T_fluid**3 * fermi_dirac(eta_rad, 2)
    BB_e_mev = g_factor * 4.0d0*PI * CLITE_CM / HC_MEVCM**3 * &
               T_fluid**4 * fermi_dirac(eta_rad, 3)
    BB_n_mev = max(BB_n_mev, 1.0d-60)
    BB_e_mev = max(BB_e_mev, 1.0d-60)

    ! ======================================================================
    !  Scattering opacity (neutral-current, same for both paths)
    !
    !  Ruffert A1, A8 (1996):
    !    κ_s = n_b * (Cs_n * Yn_eff + Cs_p * Yp_eff + C_A * Xh) * (T/me)^2 * F5/F3
    ! ======================================================================
    ! Coupling constants
    ! Cs_n = (1 + 5α^2)/24 · σ_0
    ! Cs_p = (4(Cv-1)^2 + 5α^2)/24 · σ_0   with Cv = 0.5 + 2*sin^2(θ_W)

    ! Pauli blocking for degenerate nucleons (Ruffert A8)
    ymn = xn / (1.0d0 + 2.0d0/3.0d0 * max(0.0d0, mu_n/T_fluid))
    ymp = xp / (1.0d0 + 2.0d0/3.0d0 * max(0.0d0, mu_p/T_fluid))
    if (.not. ieee_is_finite(ymn)) ymn = xn
    if (.not. ieee_is_finite(ymp)) ymp = xp

    ! Coherent nuclear scattering: C_A = σ_0/16 * Ā * (1 - Z̄/Ā)^2
    C_nucleus = 0.0d0
    if (abar1 > 0.0d0 .and. xh > 0.0d0) then
      C_nucleus = SIGMA_0 / 16.0d0 * abar1 * (1.0d0 - zbar1/abar1)**2
      if (.not. ieee_is_finite(C_nucleus)) C_nucleus = 0.0d0
    end if

    ! Energy factor for scattering: (T/me)^2 * F5(η)/F3(η)
    efac_s = (T_fluid/me_mev)**2 * &
             fermi_dirac(eta_rad, 5) / max(fermi_dirac(eta_rad, 3), 1.0d-40)

    kappa_s_nc = nb * efac_s * ( &
         (1.0d0 + 5.0d0*ALPHA**2) / 24.0d0 * SIGMA_0 * ymn + &
         (4.0d0*(CV-1.0d0)**2 + 5.0d0*ALPHA**2) / 24.0d0 * SIGMA_0 * ymp + &
         C_nucleus * xh )
    if (.not. ieee_is_finite(kappa_s_nc)) kappa_s_nc = 0.0d0
    kappa_s_nc = max(kappa_s_nc, 1.0d-60)   ! [cm^{-1}]

    ! ======================================================================
    !                   *** PATH SELECTION ***
    !
    !  #ifdef M1_RATES  →  Weakhub table path  (Path A)
    !  #else            →  Fully analytical     (Path B)
    ! ======================================================================

    {#IFDEF M1_RATES
    if(m1_rates_Weakhub) then
        !==========================================================================
        !  PATH A: Weakhub table rates
        !  Mirrors FIL calc_eas() when (M1_Weakhub::use_Weakhub_eas == true)
        !
        !  Step A1: read kappa_a, kappa_n, kappa_s from table via m1_get_eas
        !  Step A2: Kirchhoff emission Q = kappa_a * BB  (for νe, ν̄e, νx)
        !  Step A3: add analytical plasma+brems+pair for νx/νμ/ν̄μ
        !  Step A4: Kirchhoff κ_a contribution from Step A3
        !  Step A5: accumulate and apply T_nu correction
        !==========================================================================
    
        !> A1 – table lookup (already divided by LENGTHGF inside the pointer sub)
        {#IFNDEF M1_TAU_OPT
        call m1_get_eas(wrad, speciesKSP, ix^D, eas_eq, fluid_Prim)
        eas_eq(:) = eas_eq / LENGTHGF
        }
        {#IFDEF M1_TAU_OPT
        if (calc_eqRates_here) then
          call m1_get_eas(wrad, speciesKSP, ix^D, eas_eq, fluid_Prim)
          eas_eq(:) = eas_eq / LENGTHGF
        else
          eas_eq = eas   ! previously stored equilib. kappas from tau_path
        end if
        }
    
        !> Add analytical scattering on top of Weakhub scattering if desired,
        !> or just use the Weakhub value for k_s.
        !> (Weakhub already includes nucleon + nuclear scattering, so we replace
        !>  with the analytical value only when Weakhub does NOT provide k_s.)
        !> Current policy: use Weakhub k_s directly (same as FIL).
        ! eas_eq(k_s) is already set by the table call.
    
        !> A2 – Kirchhoff emission for ALL species from the table kappa_a/kappa_n
        !>  (FIL: calc_Kirchoff_emission(Q, R, kappa_a, kappa_n, true))
        !>  For νe, ν̄e:   Q = kappa_a * BB_e,  R = kappa_n * BB_n
        !>  For νx/νμ/ν̄μ: also applies (need_nux_emission = true)
        Qems   = 0.0d0
        Qems_n = 0.0d0
        ! Use BBs in code units (Black_er, Black_nr) → Q in code units:
        ! Q_code = kappa_a_code * Black_er
        ! but here we compute in MeV and convert at the end.
        ! In CGS: Q_mev = kappa_a_cgs * BB_e_mev
        ! kappa_a_code = kappa_a_cgs / LENGTHGF  →  kappa_a_cgs = kappa_a_code * LENGTHGF
        Qems_sec   = eas_eq(k_a) * LENGTHGF * BB_e_mev     ! [MeV s^-1 cm^-3]
        Qems_n_sec = eas_eq(k_n) * LENGTHGF * BB_n_mev     ! [s^-1 cm^-3]
        if (Qems_sec   > 0.0d0 .and. ieee_is_finite(Qems_sec))   Qems   = Qems_sec
        if (Qems_n_sec > 0.0d0 .and. ieee_is_finite(Qems_n_sec)) Qems_n = Qems_n_sec
    
        !> A3 – secondary emission from plasma, brems, pair for heavy leptons
        !>  (FIL: add_plasmon_decay_emission, add_brems_emission, add_pair_process_emission
        !>        accumulated into Q_tmp, R_tmp)
        if (speciesKSP == m1_i_nux .or. speciesKSP == m1_i_mu .or. &
            speciesKSP == m1_i_mubar) then
    
          Qems_sec   = 0.0d0
          Qems_n_sec = 0.0d0
    
          call m1_eas_analytic_plasmon(speciesKSP, eta_rad, eta_e, &
               ye_fluid, T_fluid, rho_fluid, eas, Qems, Qems_n)
          Qems_sec   = Qems_sec   + Qems
          Qems_n_sec = Qems_n_sec + Qems_n
    
          call m1_eas_analytic_brems(speciesKSP, eta_rad, ye_fluid, &
               T_fluid, rho_fluid, av_energy, T_nu, eas, Qems, Qems_n)
          Qems_sec   = Qems_sec   + Qems
          Qems_n_sec = Qems_n_sec + Qems_n
    
          call m1_eas_analytic_pair(speciesKSP, eta_rad, eta_e, &
               ye_fluid, T_fluid, rho_fluid, eas, Qems, Qems_n)
          Qems_sec   = Qems_sec   + Qems
          Qems_n_sec = Qems_n_sec + Qems_n
    
          !> A4 – Kirchhoff kappa from secondary emission  (FIL: add_Kirchoff_absorption_opacity)
          kappa_a_K = 0.0d0
          kappa_n_K = 0.0d0
          if (Qems_sec   > 0.0d0 .and. ieee_is_finite(Qems_sec))   &
               kappa_a_K = Qems_sec   / BB_e_mev
          if (Qems_n_sec > 0.0d0 .and. ieee_is_finite(Qems_n_sec)) &
               kappa_n_K = Qems_n_sec / BB_n_mev
    
          !> A5 – accumulate: add Kirchhoff kappas onto Weakhub kappas  [cm^{-1}]
          eas_eq(k_a) = eas_eq(k_a) + kappa_a_K
          eas_eq(k_n) = eas_eq(k_n) + kappa_n_K
    
          !> Recompute total emission including secondary processes
          Qems_sec = Qems_sec * MEV_TO_ERG_KEN * RHOGF * EPSGF / TIMEGF
          Qems_n_sec = Qems_n_sec * MNUC_CGS * RHOGF / TIMEGF
          Qems   = eas_eq(k_a) * LENGTHGF * BB_e_mev * MEV_TO_ERG_KEN * RHOGF * EPSGF / TIMEGF
          Qems_n = eas_eq(k_n) * LENGTHGF * BB_n_mev * MNUC_CGS * RHOGF / TIMEGF
        else
          ! νe, ν̄e: keep Kirchhoff emission from Weakhub kappas (Step A2)
          Qems   = Qems   * MEV_TO_ERG_KEN * RHOGF * EPSGF / TIMEGF
          Qems_n = Qems_n * MNUC_CGS * RHOGF / TIMEGF
        end if
    
        !> Load scattering from Weakhub (already in code units via eas_eq)
        !> Note: we use the analytical kappa_s only for floor safety.
        eas(k_s) = max(eas_eq(k_s), kappa_s_nc / LENGTHGF)
    
        !> Emission and absorption into eas (code units)
        eas(Q_ems)   = max(Qems,   1.0d-30)
        eas(Q_ems_n) = max(Qems_n, 1.0d-30)
        eas(k_a)     = max(eas_eq(k_a), 1.0d-30)
        eas(k_n)     = max(eas_eq(k_n), 1.0d-30)
    
     !=== end #ifdef M1_RATES (Path A) =====================================
    else if (m1_rates_analytical) then
        !==========================================================================
        !  PATH B: Fully analytical rates
        !  Mirrors FIL calc_eas() when (M1_Weakhub::use_Weakhub_eas == false)
        !
        !  Step B1: charged-current absorption kappa_a/kappa_n  (νe, ν̄e)
        !  Step B2: neutral-current scattering kappa_s          (all species)
        !  Step B3: secondary emission Q for νx/νμ/ν̄μ (brems + pair + plasmon)
        !  Step B4: Kirchhoff kappa_a for νx/νμ/ν̄μ from Step B3
        !  Step B5: Kirchhoff emission Q for νe, ν̄e from kappa_a (Step B1)
        !==========================================================================
    
        !> B1 – Charged-current absorption opacity  (FIL: add_charged_current_absorption_opacity)
        !>
        !>  ζ_nue     = η_np * σ_0/4 * (1+3α^2) / block_factor_n
        !>  ζ_nuebar  = η_pn * σ_0/4 * (1+3α^2) / block_factor_p
        !>
        !>  kappa_n[i] = ζ[i] * (T/me)^2 * F4(η_ν)/F2(η_ν)
        !>  kappa_a[i] = ζ[i] * (T/me)^2 * F5(η_ν)/F3(η_ν)
    
        kappa_a_cc = 0.0d0
        kappa_n_cc = 0.0d0
    
        select case(speciesKSP)
        case(m1_i_nue)
          ! νe + n → p + e⁻
          ! Lepton final-state blocking: 1 + exp( η_e - <ε_ν>/T )  with <ε_ν>/T = F5/F4(η_ν)
          block_nue = 1.0d0 + exp(eta_e - fermi_dirac(eta_rad, 5) / max(fermi_dirac(eta_rad, 4), 1.0d-40))
          zeta_nue  = 0.0d0
          if (eta_np > 0.0d0 .and. ieee_is_finite(block_nue) .and. block_nue > 0.0d0) then
            zeta_nue = eta_np * 0.25d0 * (1.0d0 + 3.0d0*ALPHA**2) * SIGMA_0 / block_nue
          end if
          nfac = (T_fluid/me_mev)**2 * fermi_dirac(eta_rad, 4) / max(fermi_dirac(eta_rad, 2), 1.0d-40)
          efac = (T_fluid/me_mev)**2 * fermi_dirac(eta_rad, 5) / max(fermi_dirac(eta_rad, 3), 1.0d-40)
          kappa_n_cc = max(zeta_nue * nfac, 0.0d0)
          kappa_a_cc = max(zeta_nue * efac, 0.0d0)
    
        case(m1_i_nuebar)
          ! ν̄e + p → n + e⁺
          block_nuebar = 1.0d0 + exp(-eta_e - fermi_dirac(eta_rad, 5) / max(fermi_dirac(eta_rad, 4), 1.0d-40))
          zeta_nuebar  = 0.0d0
          if (eta_pn > 0.0d0 .and. ieee_is_finite(block_nuebar) .and. block_nuebar > 0.0d0) then
            zeta_nuebar = eta_pn * 0.25d0 * (1.0d0 + 3.0d0*ALPHA**2) * SIGMA_0 / block_nuebar
          end if
          nfac = (T_fluid/me_mev)**2 * fermi_dirac(eta_rad, 4) / max(fermi_dirac(eta_rad, 2), 1.0d-40)
          efac = (T_fluid/me_mev)**2 * fermi_dirac(eta_rad, 5) / max(fermi_dirac(eta_rad, 3), 1.0d-40)
          kappa_n_cc = max(zeta_nuebar * nfac, 0.0d0)
          kappa_a_cc = max(zeta_nuebar * efac, 0.0d0)
    
        case default
          ! Heavy-lepton neutrinos have no charged-current in the analytical scheme
          kappa_a_cc = 0.0d0
          kappa_n_cc = 0.0d0
        end select
    
        if (.not. ieee_is_finite(kappa_a_cc)) kappa_a_cc = 0.0d0
        if (.not. ieee_is_finite(kappa_n_cc)) kappa_n_cc = 0.0d0
    
        ! Store (convert to code units): [cm^{-1}] / LENGTHGF → code units
        eas_eq(k_a) = max(kappa_a_cc / LENGTHGF, 1.0d-60)
        eas_eq(k_n) = max(kappa_n_cc / LENGTHGF, 1.0d-60)
        eas_eq(k_s) = max(kappa_s_nc / LENGTHGF, 1.0d-60)
    
        !> B3 – Secondary emission (νx, νμ, ν̄μ only)
        !>  (FIL: add_plasmon_decay_emission, add_brems_emission, add_pair_process_emission)
        Qems_sec   = 0.0d0
        Qems_n_sec = 0.0d0
    
        if (speciesKSP == m1_i_nux .or. speciesKSP == m1_i_mu .or. &
            speciesKSP == m1_i_mubar) then
    
          ! Plasmon decay
          call m1_eas_analytic_plasmon(speciesKSP, eta_rad, eta_e, &
               ye_fluid, T_fluid, rho_fluid, eas, Qems, Qems_n)
          Qems_sec   = Qems_sec   + Qems
          Qems_n_sec = Qems_n_sec + Qems_n
    
          ! NN Bremsstrahlung
          call m1_eas_analytic_brems(speciesKSP, eta_rad, ye_fluid, &
               T_fluid, rho_fluid, av_energy, T_nu, eas, Qems, Qems_n)
          Qems_sec   = Qems_sec   + Qems
          Qems_n_sec = Qems_n_sec + Qems_n
    
          ! Electron-positron pair annihilation
          call m1_eas_analytic_pair(speciesKSP, eta_rad, eta_e, &
               ye_fluid, T_fluid, rho_fluid, eas, Qems, Qems_n)
          Qems_sec   = Qems_sec   + Qems
          Qems_n_sec = Qems_n_sec + Qems_n
    
          Qems_sec   = max(Qems_sec,   1.0d-60)
          Qems_n_sec = max(Qems_n_sec, 1.0d-60)
    
          !> B4 – Kirchhoff absorption opacity from secondary emission
          !>  (FIL: add_Kirchoff_absorption_opacity)
          !>  kappa_a = Q_mev / BB_e_mev   [cm^{-1}]
          kappa_a_K = 0.0d0
          kappa_n_K = 0.0d0
          if (ieee_is_finite(Qems_sec)   .and. Qems_sec   > 0.0d0) &
               kappa_a_K = Qems_sec   / BB_e_mev
          if (ieee_is_finite(Qems_n_sec) .and. Qems_n_sec > 0.0d0) &
               kappa_n_K = Qems_n_sec / BB_n_mev
    
          ! Add Kirchhoff kappas for νx [cm^{-1}] converted to code units
          eas_eq(k_a) = eas_eq(k_a) + max(kappa_a_K / LENGTHGF, 0.0d0)
          eas_eq(k_n) = eas_eq(k_n) + max(kappa_n_K / LENGTHGF, 0.0d0)
    
          ! Emission in code units
          eas(Q_ems)   = max(Qems_sec   * MEV_TO_ERG_KEN * RHOGF * EPSGF / TIMEGF, 1.0d-30)
          eas(Q_ems_n) = max(Qems_n_sec * MNUC_CGS * RHOGF / TIMEGF, 1.0d-30)
    
    else
      !> B5 – Kirchhoff emission for νe, ν̄e from charged-current kappa_a
      !>  (FIL: calc_Kirchoff_emission – overwrites Q for NUE, NUE_BAR)
      !>  Q_mev = kappa_a_cgs * BB_e_mev
      !>  R     = kappa_n_cgs * BB_n_mev
      Qems   = kappa_a_cc * BB_e_mev    ! [MeV s^{-1} cm^{-3}]
      Qems_n = kappa_n_cc * BB_n_mev    ! [s^{-1} cm^{-3}]

      ! Safety checks (FIL: isfinite && > 0)
      if (.not. ieee_is_finite(Qems)   .or. Qems   <= 0.0d0) Qems   = 1.0d-60
      if (.not. ieee_is_finite(Qems_n) .or. Qems_n <= 0.0d0) Qems_n = 1.0d-60

      eas(Q_ems)   = max(Qems   * MEV_TO_ERG_KEN * RHOGF * EPSGF / TIMEGF, 1.0d-30)
      eas(Q_ems_n) = max(Qems_n * MNUC_CGS * RHOGF / TIMEGF, 1.0d-30)
    end if

    ! Copy equilibrium kappas to eas output
    eas(k_a) = max(eas_eq(k_a), 1.0d-30)
    eas(k_s) = max(eas_eq(k_s), 1.0d-30)
    eas(k_n) = max(eas_eq(k_n), 1.0d-30)

    end if   !=== end #ifndef M1_RATES (Path B) =====================================


    ! ======================================================================
    !  T_nu correction (identical logic for BOTH paths)
    !
    !  Applied only when timeCurrent >= m1_tset AND av_energy is non-trivial.
    !  Mirrors the original three correction modes exactly.
    ! ======================================================================
    if (timeCurrent >= m1_tset .and. av_energy >= m1_E_atmo*100.0d0) then

      Tnu_ratio = max(1.0d0, T_nu / T_fluid)

      if (eas_correction_standard) then
        ! ------------------------------------------------------------------
        ! Standard: scale all rates by (T_nu/T)^2 when T_nu > T_fluid
        ! Applied to νe, ν̄e  (all rates) and νx  (k_s only)
        ! ------------------------------------------------------------------
        if (speciesKSP == m1_i_nue .or. speciesKSP == m1_i_nuebar .or. speciesKSP == m1_i_mu .or. &
                 speciesKSP == m1_i_mubar) then
          eas(k_a)     = eas(k_a)     * Tnu_ratio**2
          eas(k_s)     = eas(k_s)     * Tnu_ratio**2
          eas(k_n)     = eas(k_n)     * Tnu_ratio**2
          eas(Q_ems)   = eas(Q_ems)   * Tnu_ratio**2
          eas(Q_ems_n) = eas(Q_ems_n) * Tnu_ratio**2
        else if (speciesKSP == m1_i_nux ) then
          eas(k_s) = eas(k_s) * Tnu_ratio**2
        end if

      else if (eas_correction_Tnu) then
        ! ------------------------------------------------------------------
        ! Ramp: uses (T_nu/T)^2 after 2*tset, linear (T_nu/T) before
        ! ------------------------------------------------------------------
        if (timeCurrent > 2.0d0 * m1_tset) then
          if (speciesKSP == m1_i_nue .or. speciesKSP == m1_i_nuebar .or. speciesKSP == m1_i_mu .or. &
                 speciesKSP == m1_i_mubar) then
            eas(k_a)     = eas(k_a)     * Tnu_ratio**2
            eas(k_s)     = eas(k_s)     * Tnu_ratio**2
            eas(k_n)     = eas(k_n)     * Tnu_ratio**2
          else if (speciesKSP == m1_i_nux ) then
            eas(k_s)     = eas(k_s)     * Tnu_ratio**2
          end if
        else
          if (speciesKSP == m1_i_nue .or. speciesKSP == m1_i_nuebar .or. speciesKSP == m1_i_mu .or. &
                 speciesKSP == m1_i_mubar) then
            eas(k_a)     = eas(k_a)     * Tnu_ratio
            eas(k_s)     = eas(k_s)     * Tnu_ratio
            eas(k_n)     = eas(k_n)     * Tnu_ratio
          else if (speciesKSP == m1_i_nux ) then
            eas(k_s)     = eas(k_s)     * Tnu_ratio
          end if
        end if

      else if (eas_correction_tanh) then
        ! ------------------------------------------------------------------
        ! tanh: smooth switch on/off around T_nu/T = 1
        ! ------------------------------------------------------------------
        tanh_factor = (1.0d0 + tanh((T_nu/T_fluid - 1.0d0) * 10.0d0)) * 0.5d0
        if (speciesKSP == m1_i_nue .or. speciesKSP == m1_i_nuebar .or. speciesKSP == m1_i_mu .or. &
                 speciesKSP == m1_i_mubar) then
          eas(k_a)     = eas(k_a)     * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
          eas(k_s)     = eas(k_s)     * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
          eas(k_n)     = eas(k_n)     * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
          eas(Q_ems)   = eas(Q_ems)   * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
          eas(Q_ems_n) = eas(Q_ems_n) * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
        else if (speciesKSP == m1_i_nux ) then
          eas(k_s)     = eas(k_s)     * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
          eas(k_a)     = eas(k_a)   ! unchanged for νx
          eas(k_n)     = eas(k_n)   ! unchanged for νx
        end if
      end if

    end if  ! T_nu correction

    ! ======================================================================
    !  Final floor to prevent zero/negative values reaching the solver
    ! ======================================================================
    eas(Q_ems)   = max(eas(Q_ems),   1.0d-30)
    eas(Q_ems_n) = max(eas(Q_ems_n), 1.0d-30)
    eas(k_a)     = max(eas(k_a),     1.0d-30)
    eas(k_s)     = max(eas(k_s),     1.0d-30)
    eas(k_n)     = max(eas(k_n),     1.0d-30)

  end subroutine calculate_m1_eas_ixD
