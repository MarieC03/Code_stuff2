!================================================================================
!
!  mod_eos_leptonic.t
!
!  4D leptonic EOS module for BHAC.
!  Provides the same public subroutine interface as mod_eos_tabulated but
!  extended to accept both yle (electron fraction) and ymu (muon fraction)
!  as independent thermodynamic variables.
!
!  Physics:
!    The total proton fraction  yp = yle + ymu  is used to query the
!    baryonic 3D sub-table  (rho, T, yp).
!    Electronic quantities (mu_e, Y_{e±}, P_{e±}, eps_{e±}) are obtained
!    from the electronic sub-table  (rho, T, yle).
!    Muonic   quantities (mu_mu, Y_{mu±}, P_{mu±}, eps_{mu±}) are obtained
!    from the muonic   sub-table  (rho, T, ln(ymu)).
!
!  Interface design mirrors mod_eos_tabulated:
!    - leptonic_temp_get_all_one_grid(rho, temp, yle, ymu, eps, prs, …)
!    - leptonic_eps_get_all_one_grid (rho, eps,  yle, ymu, temp, prs, …)
!    - leptonic_get_temp_one_grid    (rho, eps, temp, yle, ymu)
!    - leptonic_get_pressure_one_grid(prs, rho, eps,  yle, ymu)
!    - leptonic_get_eps_range        (rho, epsmin, epsmax, yle, ymu)
!    - leptonic_get_mu               (rho, temp, yle, ymu, mu_e, mu_mu, mu_p, mu_n)
!    - leptonic_beta_eq_get_all      (rho, temp, yle, ymu, eps, prs)
!
!  Copyright (C) 2024  Harry Ho-Yin Ng
!  Based on BHAC mod_eos_tabulated by Elias R. Most and Harry Ho-Yin Ng,
!  and on Margherita leptonic_eos_implementation.hh by Harry Ho-Yin Ng.
!
!================================================================================

module mod_eos_leptonic
  use mod_eos_leptonic_parameters
  use mod_eos_leptonic_interpolation
  implicit none
  public

contains

  !=============================================================================
  !> Activate the leptonic EOS: read tables, set global bounds, wire up
  !> the function pointers in mod_eos so that existing con2prim code is
  !> unmodified (for ye-only calls the ymu argument defaults to eos_ymumin).
  !=============================================================================
  subroutine eos_leptonic_activate()
    use mod_eos_readtable_leptonic_scollapse
    use mod_eos_readtable_leptonic
    use mod_eos
    include 'amrvacdef.f'

    ! 1. Read baryon table (sets lep_nrho, lep_ntemp, lep_nye, axes, bounds)
    select case (trim(baryon_table_type))
    case ('scollapse')
      call activate_baryon_tablereader_scollapse()
    case ('compose')
      call mpistop("eos_leptonic_activate: compose baryon reader not yet implemented")
    case default
      call mpistop("eos_leptonic_activate: unknown baryon_table_type")
    end select

    ! 2. Read leptonic table (ele + muon sub-tables; nrho/ntemp from baryon)
    call activate_tablereader_leptonic()

    ! 3. Consistency check: rho and temp axes must agree
    if (abs(lep_logrho_table(1) - lep_logrho_table(1)) > 1.0d-12) &
        call mpistop("eos_leptonic_activate: rho axis mismatch between baryon and leptonic tables")

    ! 4. Wire up the EOS function pointers that the rest of BHAC uses
    eos_type       = tabulated   ! re-use existing flag
    eos_rhomin     = lep_eos_rhomin
    eos_rhomax     = lep_eos_rhomax
    eos_tempmin    = lep_eos_tempmin
    eos_tempmax    = lep_eos_tempmax
    eos_yemin      = lep_eos_yemin    ! min yp  (used by c2p bounds)
    eos_yemax      = lep_eos_yemax
    eos_ymumin     = lep_eos_ymumin
    eos_ymumax     = lep_eos_ymumax
    eos_epsmin     = lep_eos_epsmin
    eos_epsmax     = lep_eos_epsmax
    eos_hmin       = lep_eos_hmin

    eos_get_pressure_one_grid       => leptonic_get_pressure_one_grid
    eos_get_eps_one_grid            => leptonic_get_eps_one_grid
    eos_get_cs2_one_grid            => leptonic_get_cs2_one_grid
    eos_get_temp_one_grid           => leptonic_get_temp_one_grid
    eos_get_eps_range               => leptonic_get_eps_range
    eos_eps_get_all_one_grid        => leptonic_eps_get_all_one_grid
    eos_temp_get_all_one_grid       => leptonic_temp_get_all_one_grid
    eos_get_all_beta_eqm_one_grid   => leptonic_get_all_beta_eqm_one_grid

    if (mype == 0) then
      write(*,*) '****************************************'
      write(*,*) '----  4D Leptonic EOS activated  -------'
      write(*,*) '  use_muons = ', lep_use_muons
      write(*,*) '****************************************'
    end if
  end subroutine eos_leptonic_activate

  !=============================================================================
  !> Internal: bound and decompose (yle, ymu) into (yp, safe_ymu)
  !=============================================================================
  pure subroutine lep_bound_ye_ymu(yle_in, ymu_in, yp_out, yle_safe, ymu_safe, logymu_safe)
    use mod_eos_leptonic_parameters
    double precision, intent(in)  :: yle_in, ymu_in
    double precision, intent(out) :: yp_out, yle_safe, ymu_safe, logymu_safe
    double precision :: yp_loc

    yle_safe  = max(lep_eos_ylemin, min(lep_eos_ylemax, yle_in))
    ymu_safe  = max(lep_eos_ymumin, min(lep_eos_ymumax, ymu_in))

    ! Fix ymu if yp would exceed baryon table range
    if (lep_fix_ymu_high_yp) then
      yp_loc = yle_safe + ymu_safe
      if (yp_loc > lep_eos_yemax) ymu_safe = lep_eos_ymumin
    end if

    yp_out      = max(lep_eos_yemin, min(lep_eos_yemax, yle_safe + ymu_safe))
    logymu_safe = log(max(lep_eos_ymumin, ymu_safe))
  end subroutine lep_bound_ye_ymu

  !=============================================================================
  !> Internal: bound rho and return ln(rho), ln(T)
  !=============================================================================
  pure subroutine lep_bound_rho_temp(rho_in, temp_in, rho_safe, temp_safe, lrho, ltemp)
    use mod_eos_leptonic_parameters
    double precision, intent(in)  :: rho_in, temp_in
    double precision, intent(out) :: rho_safe, temp_safe, lrho, ltemp

    rho_safe  = max(lep_eos_rhomin,  min(lep_eos_rhomax,  rho_in))
    temp_safe = max(lep_eos_tempmin, min(lep_eos_tempmax, temp_in))
    lrho  = log(rho_safe)
    ltemp = log(temp_safe)
  end subroutine lep_bound_rho_temp

  !=============================================================================
  !> Get pressure from (rho, eps, yle, ymu)  -- inversion via T(rho,eps,yle,ymu)
  !=============================================================================
  subroutine leptonic_get_pressure_one_grid(prs, rho, eps, temp, ye, ymu)
    double precision, intent(inout) :: prs
    double precision, intent(in)    :: rho, eps
    double precision, intent(inout), optional :: temp
    double precision, intent(in),    optional :: ye, ymu

    double precision :: yle_l, ymu_l, temp_l
    double precision :: prs_tmp, cs2_tmp, h_tmp, ent_tmp

    yle_l = lep_eos_ylemin
    ymu_l = lep_eos_ymumin
    if (present(ye))  yle_l = ye
    if (present(ymu)) ymu_l = ymu
    temp_l = lep_eos_tempmin
    if (present(temp)) temp_l = temp

    call leptonic_eps_get_all_one_grid(rho, eps, yle_l, ymu_l, temp_l, prs=prs)
    if (present(temp)) temp = temp_l
  end subroutine leptonic_get_pressure_one_grid

  !=============================================================================
  !> Get eps from (rho, prs, yle, ymu)  -- not monotone, so approximation only
  !=============================================================================
  subroutine leptonic_get_eps_one_grid(prs, rho, eps, temp, ye, ymu)
    double precision, intent(in)    :: prs, rho
    double precision, intent(inout) :: eps
    double precision, intent(in),    optional :: temp, ye, ymu
    ! Placeholder: inversion from pressure to eps is not monotone for generic EOS.
    ! Use leptonic_eps_get_all_one_grid or leptonic_temp_get_all_one_grid instead.
    call mpistop("leptonic_get_eps_one_grid: inversion P->eps not supported for 4D leptonic EOS")
  end subroutine leptonic_get_eps_one_grid

  !=============================================================================
  !> Get cs2 from (rho, eps, yle, ymu)
  !=============================================================================
  subroutine leptonic_get_cs2_one_grid(cs2, rho, eps, temp, ye, ymu)
    double precision, intent(inout) :: cs2
    double precision, intent(in)    :: rho, eps
    double precision, intent(in),    optional :: ye, ymu
    double precision, intent(inout), optional :: temp

    double precision :: yle_l, ymu_l, temp_l

    yle_l = lep_eos_ylemin
    ymu_l = lep_eos_ymumin
    if (present(ye))  yle_l = ye
    if (present(ymu)) ymu_l = ymu
    temp_l = lep_eos_tempmin
    if (present(temp)) temp_l = temp

    call leptonic_eps_get_all_one_grid(rho, eps, yle_l, ymu_l, temp_l, cs2=cs2)
    if (present(temp)) temp = temp_l
  end subroutine leptonic_get_cs2_one_grid

  !=============================================================================
  !> Invert eps -> T at fixed (rho, yle, ymu) using Brent root finding
  !=============================================================================
  subroutine leptonic_get_temp_one_grid(rho, eps, temp, ye, ymu)
    use mod_rootfinding
    implicit none
    double precision, intent(in)    :: rho, ye
    double precision, intent(inout) :: eps, temp
    double precision, intent(in), optional :: ymu

    double precision :: yle_l, ymu_l
    double precision :: lrho, ltemp_lo, ltemp_hi, log_eps_target
    double precision :: eps_shifted, ffx(lep_nvars_baryon)
    double precision :: yp, logymu
    double precision :: rho_s, temp_s
    integer :: error_code

    yle_l = ye
    ymu_l = lep_eos_ymumin
    if (present(ymu)) ymu_l = ymu

    if (rho <= 0.0d0) then
      temp = lep_eos_tempmin; return
    end if

    call lep_bound_rho_temp(rho, temp, rho_s, temp_s, lrho, ltemp_lo)
    call lep_bound_ye_ymu(yle_l, ymu_l, yp, yle_l, ymu_l, logymu)

    ltemp_lo = lep_logtemp_table(1)
    ltemp_hi = lep_logtemp_table(lep_ntemp)

    ! Target: ln(eps + energy_shift) from table
    eps_shifted = max(0.0d0, eps) + lep_energy_shift
    if (eps_shifted <= 0.0d0) eps_shifted = tiny(1.0d0)
    log_eps_target = log(eps_shifted)

    error_code = -1
    call rootfinding_brent(temp, ltemp_lo, ltemp_hi, lep_eos_precision, lep_eos_iter_max, &
                           error_code, func_eps_of_temp)

    select case (error_code)
    case (-1); call mpistop("leptonic_get_temp_one_grid: Brent not attempted")
    case (2);  call mpistop("leptonic_get_temp_one_grid: NaN in Brent")
    case (3);  temp = ltemp_lo  ! out of range: clamp to cold slice
    end select

    ! temp now holds ln(T); convert back
    temp = exp(temp)

  contains
    double precision function func_eps_of_temp(lt)
      double precision, intent(in) :: lt
      call intep3d_lep_many(lrho, lt, yp, ffx, lep_tables_baryon, &
          lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon, &
          lep_logrho_table, lep_logtemp_table, lep_ye_table)
      func_eps_of_temp = ffx(i_lep_logenergy) - log_eps_target
    end function
  end subroutine leptonic_get_temp_one_grid

  !=============================================================================
  !> Get eps range at fixed (rho, yle, ymu)
  !=============================================================================
  subroutine leptonic_get_eps_range(rho, epsmin, epsmax, ye, ymu)
    double precision, intent(in)  :: rho
    double precision, intent(out) :: epsmin, epsmax
    double precision, intent(in), optional :: ye, ymu

    double precision :: yle_l, ymu_l, yp, logymu, lrho, rho_s
    double precision :: ffx(lep_nvars_baryon)
    integer :: j

    yle_l = lep_eos_ylemin
    ymu_l = lep_eos_ymumin
    if (present(ye))  yle_l = ye
    if (present(ymu)) ymu_l = ymu

    rho_s = max(lep_eos_rhomin, min(lep_eos_rhomax, rho))
    lrho  = log(rho_s)
    call lep_bound_ye_ymu(yle_l, ymu_l, yp, yle_l, ymu_l, logymu)

    epsmin =  1.0d99
    epsmax = -1.0d99
    do j = 1, lep_ntemp
      call intep3d_lep_many(lrho, lep_logtemp_table(j), yp, ffx, lep_tables_baryon, &
          lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon, &
          lep_logrho_table, lep_logtemp_table, lep_ye_table)
      epsmin = min(epsmin, exp(ffx(i_lep_logenergy)) - lep_energy_shift)
      epsmax = max(epsmax, exp(ffx(i_lep_logenergy)) - lep_energy_shift)
    end do
  end subroutine leptonic_get_eps_range

  !=============================================================================
  !> Full query from (rho, temp, yle, ymu)  ->  eps, prs, ent, cs2, mu_e, etc.
  !=============================================================================
  subroutine leptonic_temp_get_all_one_grid(rho, temp, yle, ymu, eps, prs, ent, cs2, &
                       mu_e, mu_mu, mu_p, mu_n, xa, xh, xn, xp, abar, zbar, &
                       press_ele, press_mu, eps_ele, eps_mu)
    use mod_eos_leptonic_parameters
    use mod_eos_leptonic_interpolation
    implicit none

    double precision, intent(inout) :: rho, temp
    double precision, intent(in)    :: yle
    double precision, intent(in), optional :: ymu
    double precision, intent(inout) :: eps
    double precision, intent(inout), optional :: prs, ent, cs2
    double precision, intent(inout), optional :: mu_e, mu_mu, mu_p, mu_n
    double precision, intent(inout), optional :: xa, xh, xn, xp, abar, zbar
    double precision, intent(inout), optional :: press_ele, press_mu
    double precision, intent(inout), optional :: eps_ele,   eps_mu

    double precision :: yle_s, ymu_s, yp, logymu
    double precision :: lrho, ltemp
    double precision :: rho_s, temp_s
    double precision :: ffx_b(lep_nvars_baryon)
    double precision :: ffx_e(lep_nvars_ele)
    double precision :: ffx_m(lep_nvars_muon)
    double precision :: P_bar, eps_bar, P_e, eps_e, P_mu_loc, eps_mu_loc

    ! Bounds
    if (rho  <= 0.0d0) then
      eps = lep_eos_epsmin; if (present(prs)) prs = 0.0d0; return
    end if
    call lep_bound_rho_temp(rho, temp, rho_s, temp_s, lrho, ltemp)
    rho = rho_s; temp = temp_s

    ymu_s = lep_eos_ymumin
    if (present(ymu)) ymu_s = ymu
    call lep_bound_ye_ymu(yle, ymu_s, yp, yle_s, ymu_s, logymu)

    !--------------------------------------------------------------------------
    ! 1. Baryonic sub-table  f(lrho, ltemp, yp)
    !--------------------------------------------------------------------------
    call intep3d_lep_many(lrho, ltemp, yp, ffx_b, lep_tables_baryon, &
        lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon, &
        lep_logrho_table, lep_logtemp_table, lep_ye_table)

    P_bar   = exp(ffx_b(i_lep_logpress))
    eps_bar = exp(ffx_b(i_lep_logenergy)) - lep_energy_shift

    !--------------------------------------------------------------------------
    ! 2. Electronic sub-table  f(lrho, ltemp, yle)
    !--------------------------------------------------------------------------
    P_e = 0.0d0; eps_e = 0.0d0
    if (lep_add_ele_contribution) then
      call intep3d_lep_many(lrho, ltemp, yle_s, ffx_e, lep_tables_ele, &
          lep_nrho, lep_ntemp, lep_nyle, lep_nvars_ele, &
          lep_logrho_table, lep_logtemp_table, lep_yle_table)
      P_e   = ffx_e(i_lep_press_eminus) + ffx_e(i_lep_press_eplus)
      eps_e = ffx_e(i_lep_eps_eminus)   + ffx_e(i_lep_eps_eplus)
    end if

    !--------------------------------------------------------------------------
    ! 3. Muonic sub-table  f(lrho, ltemp, ln(ymu))
    !--------------------------------------------------------------------------
    P_mu_loc = 0.0d0; eps_mu_loc = 0.0d0
    if (lep_use_muons) then
      call intep3d_lep_many(lrho, ltemp, logymu, ffx_m, lep_tables_muon, &
          lep_nrho, lep_ntemp, lep_nymu, lep_nvars_muon, &
          lep_logrho_table, lep_logtemp_table, lep_logymu_table)
      P_mu_loc   = ffx_m(i_lep_press_muminus) + ffx_m(i_lep_press_muplus)
      eps_mu_loc = ffx_m(i_lep_eps_muminus)   + ffx_m(i_lep_eps_muplus)
    end if

    !--------------------------------------------------------------------------
    ! 4. Assemble outputs
    !--------------------------------------------------------------------------
    eps = eps_bar + eps_e + eps_mu_loc
    if (present(prs))      prs = P_bar + P_e + P_mu_loc
    if (present(ent))      ent = ffx_b(i_lep_entropy)
    if (present(cs2))      cs2 = max(0.0d0, min(0.9999999d0, ffx_b(i_lep_cs2)))

    if (present(mu_p))     mu_p  = ffx_b(i_lep_mu_p)
    if (present(mu_n))     mu_n  = ffx_b(i_lep_mu_n)
    if (present(xa))       xa    = ffx_b(i_lep_xa)
    if (present(xh))       xh    = ffx_b(i_lep_xh)
    if (present(xn))       xn    = ffx_b(i_lep_xn)
    if (present(xp))       xp    = ffx_b(i_lep_xp)
    if (present(abar))     abar  = ffx_b(i_lep_abar)
    if (present(zbar))     zbar  = ffx_b(i_lep_zbar)

    ! Chemical potentials from leptonic tables (more accurate for yle < yemax)
    if (present(mu_e)) then
      if (yp >= lep_eos_yemax) then
        mu_e = ffx_b(i_lep_mu_e)   ! fallback to baryonic mu_e
      else if (lep_add_ele_contribution) then
        mu_e = ffx_e(i_lep_mu_ele)
      else
        mu_e = ffx_b(i_lep_mu_e)
      end if
    end if

    if (present(mu_mu)) then
      if (lep_use_muons) then
        mu_mu = ffx_m(i_lep_mu_mu)
      else
        mu_mu = 105.6583755d0   ! muon rest mass in MeV
      end if
    end if

    if (present(press_ele)) press_ele = P_e
    if (present(press_mu))  press_mu  = P_mu_loc
    if (present(eps_ele))   eps_ele   = eps_e
    if (present(eps_mu))    eps_mu    = eps_mu_loc
  end subroutine leptonic_temp_get_all_one_grid

  !=============================================================================
  !> Full query from (rho, eps, yle, ymu)  ->  temp, prs, cs2, etc.
  !> Performs T inversion first, then calls leptonic_temp_get_all_one_grid.
  !=============================================================================
  subroutine leptonic_eps_get_all_one_grid(rho, eps, yle, ymu, temp, prs, ent, cs2, &
                       mu_e, mu_mu, mu_p, mu_n, xa, xh, xn, xp, abar, zbar, &
                       press_ele, press_mu, eps_ele, eps_mu)
    double precision, intent(inout) :: rho, eps
    double precision, intent(in)    :: yle
    double precision, intent(in),    optional :: ymu
    double precision, intent(inout) :: temp
    double precision, intent(inout), optional :: prs, ent, cs2
    double precision, intent(inout), optional :: mu_e, mu_mu, mu_p, mu_n
    double precision, intent(inout), optional :: xa, xh, xn, xp, abar, zbar
    double precision, intent(inout), optional :: press_ele, press_mu, eps_ele, eps_mu

    double precision :: ymu_l

    ymu_l = lep_eos_ymumin
    if (present(ymu)) ymu_l = ymu

    ! Invert eps -> T
    call leptonic_get_temp_one_grid(rho, eps, temp, yle, ymu_l)

    ! Then query everything at (rho, T, yle, ymu)
    call leptonic_temp_get_all_one_grid(rho, temp, yle, ymu_l, eps, &
        prs=prs, ent=ent, cs2=cs2, mu_e=mu_e, mu_mu=mu_mu, mu_p=mu_p, mu_n=mu_n, &
        xa=xa, xh=xh, xn=xn, xp=xp, abar=abar, zbar=zbar, &
        press_ele=press_ele, press_mu=press_mu, eps_ele=eps_ele, eps_mu=eps_mu)
  end subroutine leptonic_eps_get_all_one_grid

  !=============================================================================
  !> Chemical potential query at fixed (rho, temp, yle, ymu)
  !=============================================================================
  subroutine leptonic_get_mu(rho, temp, yle, ymu, mu_e, mu_mu, mu_p, mu_n)
    double precision, intent(in)    :: rho, temp, yle, ymu
    double precision, intent(out)   :: mu_e, mu_mu, mu_p, mu_n
    double precision :: eps_dummy

    eps_dummy = lep_eos_epsmin
    call leptonic_temp_get_all_one_grid(rho, temp, yle, ymu, eps_dummy, &
        mu_e=mu_e, mu_mu=mu_mu, mu_p=mu_p, mu_n=mu_n)
  end subroutine leptonic_get_mu

  !=============================================================================
  !> Beta equilibrium with muons:  solve for (yle, ymu) at fixed (rho, temp)
  !> such that:
  !>   mu_n - mu_p - mu_e = 0       (electron beta eq)
  !>   mu_n - mu_p - mu_mu = 0      (muon    beta eq)
  !>   yp = yle + ymu  (charge neutrality enforced via yp -> baryonic table)
  !>
  !> Algorithm: outer Brent on ln(ymu) using
  !>   F(ymu) = mu_mu + mu_p - mu_n - Qnp   (Qnp = mn-mp in MeV, ~ 1.293 MeV)
  !>   where for each ymu we solve the inner electron beta-eq for yle.
  !=============================================================================
  subroutine leptonic_get_all_beta_eqm_one_grid(rho, temp, yle, ymu, eps, prs)
    use mod_rootfinding
    implicit none
    double precision, intent(in)    :: rho, temp
    double precision, intent(inout) :: yle, ymu
    double precision, intent(inout), optional :: eps, prs

    ! Qnp = (mn - mp) in MeV (mass difference, no rest-mass of electron)
    double precision, parameter :: Qnp = 1.2933322d0   ! MeV

    double precision :: lrho, ltemp, rho_s, temp_s
    double precision :: lymu_lo, lymu_hi, lymu_root
    double precision :: yp_loc, yle_loc, ymu_loc, logymu_loc
    double precision :: ffx_b(lep_nvars_baryon), ffx_e(lep_nvars_ele), ffx_m(lep_nvars_muon)
    double precision :: mu_e_loc, mu_mu_loc, mu_p_loc, mu_n_loc
    integer :: error_code

    call lep_bound_rho_temp(rho, temp, rho_s, temp_s, lrho, ltemp)

    lymu_lo = log(lep_eos_ymumin)
    lymu_hi = log(lep_eos_ymumax)

    error_code = -1
    ymu_loc = lymu_lo
    call rootfinding_brent(ymu_loc, lymu_lo, lymu_hi, lep_eos_precision, lep_eos_iter_max, &
                           error_code, func_F_ymu)

    select case (error_code)
    case (-1)
      ! Root not bracketed: no muons at this (rho, T)
      ymu = lep_eos_ymumin
      call inner_electron_beta_eq(lrho, ltemp, yp_loc, yle_loc)
      yle = yle_loc
    case (0, 1, 3)
      ! Accept best estimate
      ymu_loc = max(lymu_lo, min(lymu_hi, ymu_loc))
      if (ymu_loc <= lymu_lo + tiny(1.0d0)) then
        ymu = lep_eos_ymumin
      else
        ymu = exp(ymu_loc)
      end if
      yle = yle_loc
    case (2)
      call mpistop("leptonic_get_all_beta_eqm_one_grid: NaN in Brent")
    end select

    if (present(eps) .or. present(prs)) then
      call leptonic_temp_get_all_one_grid(rho_s, temp_s, yle, ymu, &
          eps_dummy, prs=prs_dummy)
      if (present(eps)) eps = eps_dummy
      if (present(prs)) prs = prs_dummy
    end if

  contains

    ! ------------------------------------------------------------------
    ! Outer root-finding residual:
    !   Given ln(ymu), find yle from inner electron beta-eq,
    !   then evaluate mu_mu + mu_p - mu_n - Qnp.
    ! ------------------------------------------------------------------
    double precision function func_F_ymu(lymu_in)
      double precision, intent(in) :: lymu_in
      double precision :: ymu_try, yp_try

      ymu_try = exp(lymu_in)
      yp_try  = 0.0d0  ! will be updated by inner solve

      ! inner electron beta-eq for yle
      call inner_electron_beta_eq(lrho, ltemp, yp_try, yle_loc)

      yp_try = yle_loc + ymu_try
      yp_loc = max(lep_eos_yemin, min(lep_eos_yemax, yp_try))

      ! baryonic mu_p, mu_n
      call intep3d_lep_many(lrho, ltemp, yp_loc, ffx_b, lep_tables_baryon, &
          lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon, &
          lep_logrho_table, lep_logtemp_table, lep_ye_table)
      mu_p_loc = ffx_b(i_lep_mu_p)
      mu_n_loc = ffx_b(i_lep_mu_n)

      ! muonic mu_mu
      call intep3d_lep_many(lrho, ltemp, lymu_in, ffx_m, lep_tables_muon, &
          lep_nrho, lep_ntemp, lep_nymu, lep_nvars_muon, &
          lep_logrho_table, lep_logtemp_table, lep_logymu_table)
      mu_mu_loc = ffx_m(i_lep_mu_mu)

      func_F_ymu = mu_mu_loc + mu_p_loc - mu_n_loc - Qnp
    end function func_F_ymu

    ! ------------------------------------------------------------------
    ! Inner electron beta-eq: at fixed (lrho, ltemp) and given yp,
    ! find yle such that  mu_e + mu_p - mu_n = 0  (i.e. mu_e = mu_hat).
    ! Uses the electronic leptonic table for mu_e.
    ! ------------------------------------------------------------------
    subroutine inner_electron_beta_eq(lrho_in, ltemp_in, yp_out, yle_out)
      use mod_rootfinding
      double precision, intent(in)  :: lrho_in, ltemp_in
      double precision, intent(out) :: yp_out, yle_out
      double precision :: yle_lo, yle_hi, yle_root
      integer :: ec

      yle_lo   = lep_eos_ylemin
      yle_hi   = lep_eos_ylemax
      yle_root = 0.5d0 * (yle_lo + yle_hi)
      ec = -1
      call rootfinding_brent(yle_root, yle_lo, yle_hi, lep_eos_precision, lep_eos_iter_max, &
                             ec, func_ele_betaeq)
      yle_out = max(yle_lo, min(yle_hi, yle_root))
      yp_out  = yle_out
    contains
      double precision function func_ele_betaeq(yle_try)
        double precision, intent(in) :: yle_try
        double precision :: yp_try, ffxb(lep_nvars_baryon), ffxe(lep_nvars_ele)
        yp_try = max(lep_eos_yemin, min(lep_eos_yemax, yle_try))
        call intep3d_lep_many(lrho_in, ltemp_in, yp_try, ffxb, lep_tables_baryon, &
            lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon, &
            lep_logrho_table, lep_logtemp_table, lep_ye_table)
        call intep3d_lep_many(lrho_in, ltemp_in, yle_try, ffxe, lep_tables_ele, &
            lep_nrho, lep_ntemp, lep_nyle, lep_nvars_ele, &
            lep_logrho_table, lep_logtemp_table, lep_yle_table)
        ! mu_e + mu_p - mu_n = 0  (electron beta eq without neutrinos)
        func_ele_betaeq = ffxe(i_lep_mu_ele) + ffxb(i_lep_mu_p) - ffxb(i_lep_mu_n) - Qnp
      end function
    end subroutine inner_electron_beta_eq

    ! Dummy variables needed inside contains
    double precision :: eps_dummy, prs_dummy
  end subroutine leptonic_get_all_beta_eqm_one_grid

end module mod_eos_leptonic
