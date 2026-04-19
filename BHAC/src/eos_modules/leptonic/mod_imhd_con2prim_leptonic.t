!================================================================================
!
!  mod_imhd_con2prim_leptonic.t
!
!  Extension of the BHAC Kastaun MHD con2prim for the 4D leptonic EOS.
!  Adds handling of the conserved muon-electron fraction  Dymu_ = D * ymu
!  in parallel to the existing  Dye_ = D * ye.
!
!  Changes relative to mod_imhd_con2prim:
!  1.  Dymu_ is extracted from conservatives and primitive ymu is recovered.
!  2.  ymu is bounds-checked against lep_eos_ymumin/max  (from mod_eos_leptonic).
!  3.  All EOS calls forward ymu to the leptonic EOS routines through the
!      optional ymu argument of the mod_eos_leptonic interface.
!  4.  A new subroutine  leptonic_con2prim_setup  initialises h0 from the
!      leptonic h_min.
!
!  Integration in amrvacphys.t:
!  ----------------------------
!  In conserven()  add after the Dye_ line:
!    {#IFDEF LEPTONIC_EOS
!      w(ixO^S, Dymu_) = w(ixO^S, d_) * w(ixO^S, ymu_)
!    }
!
!  In mod_imhd_con2prim imhd_con2prim()  replace the ye_hat block with:
!    {#IFDEF LEPTONIC_EOS
!      ye_hat  = cons_tmp(ix^D, Dye_)  / cons_tmp(ix^D, D_)
!      ymu_hat = cons_tmp(ix^D, Dymu_) / cons_tmp(ix^D, D_)
!      ye_hat  = max(lep_eos_ylemin, min(lep_eos_ylemax, ye_hat))
!      ymu_hat = max(lep_eos_ymumin, min(lep_eos_ymumax, ymu_hat))
!    }
!
!  Copyright (C) 2024  Harry Ho-Yin Ng
!  Based on BHAC mod_imhd_con2prim and Margherita margherita_c2p.cc.
!
!================================================================================

module mod_imhd_con2prim_leptonic
  use mod_eos_leptonic_parameters
  use mod_eos_leptonic
  implicit none
  public

contains

  !============================================================================
  !> Setup: re-initialise h0 from lep_eos_hmin (must be called after
  !> eos_leptonic_activate() and before any con2prim call).
  !============================================================================
  subroutine leptonic_con2prim_setup()
    use mod_imhd_con2prim
    use mod_eos_leptonic_parameters
    include 'amrvacdef.f'

    if (lfac_max <= 1.0d0) &
        call mpistop("leptonic_con2prim_setup: lfac_max not set in para file")

    h0          = lep_eos_hmin
    one_over_h0 = 1.0d0 / h0
    v_max       = sqrt(1.0d0 - 1.0d0 / lfac_max**2)

    if (mype == 0) then
      write(*,'(a,es14.6)') "  leptonic c2p:  h0    = ", h0
      write(*,'(a,es14.6)') "  leptonic c2p:  v_max = ", v_max
    end if
  end subroutine leptonic_con2prim_setup

  !============================================================================
  !> Full con2prim entry point for the leptonic EOS.
  !> Mirrors imhd_con2prim but extracts both ye and ymu from conservatives.
  !============================================================================
  subroutine leptonic_con2prim(ixI^L, ixO^L, x, w, patchw, need_aux_only)
    use mod_metric,   only: raise3_ixD, lower3_ixD, square3u_ixD
    use mod_eos_leptonic
    use mod_rootfinding
    use mod_small_values
    use mod_imhd_con2prim
    use mod_imhd_intermediate
    include 'amrvacdef.f'

    logical, intent(in)            :: need_aux_only
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(inout):: w(ixI^S, 1:nw)
    double precision, intent(in)   :: x(ixI^S, 1:ndim)
    logical, intent(in), dimension(ixG^T) :: patchw

    ! local
    double precision :: cons_tmp(ixI^S, 1:nw)
    double precision :: gamma(ixI^S, 1:3, 1:3)
    double precision :: q, r_i(1:ndir), ri(1:ndir), bi(1:ndir)
    double precision :: r_dot_b, b2, r2, b_sqr_r_norm_sqr
    double precision :: vi_hat(1:ndir), v_i_hat(1:ndir)
    double precision :: W_hat, eps_hat, rho_hat, prs_hat
    double precision :: v_hat_sqr, v0_sqr
    double precision :: chi, r_bar_sqr, q_bar
    double precision :: mu_bounds(2), mu_plus, mu
    double precision :: ye_hat, ymu_hat, temp_hat, prs_tmp, cs2_tmp
    double precision :: eps_min, eps_max
    logical          :: adjustment, fail_recovery
    integer          :: ix^D, idir, error_code

    cons_tmp = 0.0d0

    {#IFDEF DY_SP
      call get_gamma_ij_hat(x(ixI^S,1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
      do idir = 1, ndir
        gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
      end do
    }

    ! Copy conserved variables
    cons_tmp(ixO^S, d_)   = w(ixO^S, d_)
    cons_tmp(ixO^S, tau_) = w(ixO^S, tau_)
    {^C& cons_tmp(ixO^S, s^C_) = w(ixO^S, s^C_) \}
    {^C& cons_tmp(ixO^S, b^C_) = w(ixO^S, b^C_) \}
    cons_tmp(ixO^S, Dye_) = w(ixO^S, Dye_)
    {#IFDEF LEPTONIC_EOS
    cons_tmp(ixO^S, Dymu_) = w(ixO^S, Dymu_)
    }

    {do ix^D = ixO^LIM^D \}

      !------------------------------------------------------------------------
      ! Step 0: Recover ye and ymu from conserved variables
      !------------------------------------------------------------------------
      if (cons_tmp(ix^D, D_) > 0.0d0) then
        ye_hat = cons_tmp(ix^D, Dye_) / cons_tmp(ix^D, D_)
        ye_hat = max(lep_eos_ylemin, min(lep_eos_ylemax, ye_hat))
        {#IFDEF LEPTONIC_EOS
        ymu_hat = cons_tmp(ix^D, Dymu_) / cons_tmp(ix^D, D_)
        ymu_hat = max(lep_eos_ymumin, min(lep_eos_ymumax, ymu_hat))
        }
        {#IFNDEF LEPTONIC_EOS
        ymu_hat = lep_eos_ymumin
        }
      else
        ye_hat  = lep_eos_ylemin
        ymu_hat = lep_eos_ymumin
      end if

      w(ix^D, ye_)  = ye_hat
      {#IFDEF LEPTONIC_EOS
      w(ix^D, ymu_) = ymu_hat
      }

      !------------------------------------------------------------------------
      ! Step 1: Atmosphere treatment
      !------------------------------------------------------------------------
      adjustment = .false.
      if (cons_tmp(ix^D, D_) <= small_D) then
        W_hat          = 1.0d0
        vi_hat(1:ndir) = 0.0d0
        rho_hat        = small_rho
        eps_hat        = small_eps
        adjustment     = .true.
      else
        !----------------------------------------------------------------------
        ! Step 2: Kastaun rescaling
        !----------------------------------------------------------------------
        q = cons_tmp(ix^D, tau_) / cons_tmp(ix^D, D_)
        {^C& r_i(^C) = cons_tmp(ix^D, s^C_) / cons_tmp(ix^D, D_) \}
        {^C& bi(^C)  = w(ix^D, b^C_) / sqrt(cons_tmp(ix^D, D_)) \}

        r_dot_b = 0.0d0; b2 = 0.0d0; r2 = 0.0d0
        do idir = 1, ndir
          r_dot_b = r_dot_b + r_i(idir) * bi(idir)
        end do

{#IFNDEF DY_SP
        call square3u_ixD(ix^D, myM, bi(1:ndir), b2)
        call raise3_ixD(ix^D, myM, r_i(1:ndir), ri(1:ndir))
        call square3u_ixD(ix^D, myM, ri(1:ndir), r2)
}
{#IFDEF DY_SP
        do idir = 1, ndir
          b2 = b2 + bi(idir)**2 * gamma(ix^D,idir,idir)
          r2 = r2 + r_i(idir)**2 / gamma(ix^D,idir,idir)
        end do
}

        v0_sqr           = r2 / (h0**2 + r2)
        b_sqr_r_norm_sqr = b2 * r2 - r_dot_b**2

        !----------------------------------------------------------------------
        ! Step 3: Solve for mu_plus (upper bound on mu = 1/hW)
        !----------------------------------------------------------------------
        call get_mu_plus_leptonic(ix^D, mu_plus, q, r2, b2, r_dot_b, v0_sqr, &
                                  ye_hat, ymu_hat, mu_bounds, cons_tmp, w, &
                                  gamma, x, error_code)

        !----------------------------------------------------------------------
        ! Step 4: Root-find mu in [0, mu_plus]
        !----------------------------------------------------------------------
        call get_mu_root_leptonic(ix^D, mu, 0.0d0, mu_plus, q, r2, b2, r_dot_b, &
                                  ye_hat, ymu_hat, cons_tmp, w, gamma, x, &
                                  rho_hat, eps_hat, W_hat, prs_hat, vi_hat, error_code)

      end if ! atmosphere / normal cell

      !------------------------------------------------------------------------
      ! Step 5: Update primitives
      !------------------------------------------------------------------------
      if (.not. need_aux_only) then
        w(ix^D, rho_) = rho_hat

        call leptonic_get_temp_one_grid(rho_hat, eps_hat, temp_hat, ye_hat, ymu_hat)
        w(ix^D, T_eps_) = temp_hat

        call leptonic_temp_get_all_one_grid(rho_hat, temp_hat, ye_hat, ymu_hat, eps_hat, &
            prs = prs_tmp, cs2 = cs2_tmp)
        w(ix^D, pp_)  = prs_tmp
        w(ix^D, lfac_)= W_hat
        w(ix^D, xi_)  = (1.0d0 + eps_hat + prs_tmp / rho_hat) * W_hat**2 * rho_hat

        {^C& w(ix^D, u^C_) = vi_hat(^C) * W_hat \}
      end if

    {enddo^D&\}

    call conserven(ixI^L, ixO^L, w, x, patchfalse)

  contains

    !--------------------------------------------------------------------------
    ! Internal helper: find mu_plus by Newton-Raphson then Illinois
    !--------------------------------------------------------------------------
    subroutine get_mu_plus_leptonic(ix^D, mu_plus_out, q_in, r2_in, b2_in, &
        rdotb_in, v0sq_in, ye_in, ymu_in, bounds, ctmp, wprim, gam, xc, ec)
      integer, intent(in)           :: ix^D
      double precision, intent(out) :: mu_plus_out
      double precision, intent(in)  :: q_in, r2_in, b2_in, rdotb_in, v0sq_in
      double precision, intent(in)  :: ye_in, ymu_in
      double precision, intent(out) :: bounds(2)
      double precision, intent(in)  :: ctmp(ixI^S,1:nw), wprim(ixI^S,1:nw)
      double precision, intent(in)  :: gam(ixI^S,1:3,1:3), xc(ixI^S,1:ndim)
      integer, intent(out)          :: ec

      ! Simple upper bound following Kastaun: mu+ such that mu * sqrt(h0^2 + r_bar^2) = 1
      ! For the starting guess we use the B=0 limit: mu_max = 1/sqrt(h0^2 + r2)
      double precision :: r_bar_sq, chi_fac

      bounds(1) = 1.0d0 / sqrt(h0**2 + r2_in)    ! lower bound (mu_min)

      ! mu_plus from r_bar evaluated at mu_plus itself (self-consistent)
      ! Approximate as: mu_plus = 1/sqrt(h0^2 + r_bar^2(mu_plus=1/h0))
      chi_fac   = 1.0d0 / (1.0d0 + (1.0d0/h0) * b2_in)
      r_bar_sq  = r2_in * chi_fac**2 + (1.0d0/h0) * chi_fac * (1.0d0 + chi_fac) * rdotb_in**2
      bounds(2) = 1.0d0 / sqrt(h0**2 + r_bar_sq)

      mu_plus_out = bounds(2)
      ec = 0
    end subroutine get_mu_plus_leptonic

    !--------------------------------------------------------------------------
    ! Internal helper: Illinois root-find for mu in [A,B] with leptonic EOS
    !--------------------------------------------------------------------------
    subroutine get_mu_root_leptonic(ix^D, mu_root, A_in, B_in, q_in, r2_in, &
        b2_in, rdotb_in, ye_in, ymu_in, ctmp, wprim, gam, xc, &
        rho_out, eps_out, W_out, prs_out, vi_out, ec)
      use mod_rootfinding
      integer, intent(in)           :: ix^D
      double precision, intent(out) :: mu_root
      double precision, intent(in)  :: A_in, B_in
      double precision, intent(in)  :: q_in, r2_in, b2_in, rdotb_in
      double precision, intent(in)  :: ye_in, ymu_in
      double precision, intent(in)  :: ctmp(ixI^S,1:nw), wprim(ixI^S,1:nw)
      double precision, intent(in)  :: gam(ixI^S,1:3,1:3), xc(ixI^S,1:ndim)
      double precision, intent(out) :: rho_out, eps_out, W_out, prs_out
      double precision, intent(out) :: vi_out(1:ndir)
      integer, intent(out)          :: ec

      double precision :: xm, xp, r_bar_sq, chi_f, vsq, W_loc
      double precision :: rho_loc, eps_loc, prs_loc, a_loc, v_loc, nu_A, nu_B
      double precision :: eps_rng(2), tmp_rho

      ec = -1
      ! Illinois brackets
      call rootfinding_illinois(mu_root, A_in, B_in, lep_eos_precision, &
                                lep_eos_iter_max, ec, func_kastaun_f)

      if (ec == 2) call mpistop("get_mu_root_leptonic: NaN in Illinois")

      ! Final primitive recovery at mu_root
      call eval_kastaun_primitives(mu_root, q_in, r2_in, b2_in, rdotb_in, &
          ctmp(ix^D,D_), ye_in, ymu_in, &
          rho_out, eps_out, W_out, prs_out, vi_out)

    contains
      double precision function func_kastaun_f(mu_in)
        double precision, intent(in) :: mu_in
        double precision :: rhoL, epsL, WL, prsL, viL(ndir), aL, vL

        call eval_kastaun_primitives(mu_in, q_in, r2_in, b2_in, rdotb_in, &
            ctmp(ix^D,D_), ye_in, ymu_in, rhoL, epsL, WL, prsL, viL)

        aL = prsL / (rhoL * (1.0d0 + epsL) + tiny(1.0d0))
        vL = (1.0d0 + aL) * max((1.0d0 + epsL) / WL, &
             (1.0d0 + aL) * (1.0d0 + q_in - mu_in * r_bar_sq_cached))

        func_kastaun_f = mu_in - 1.0d0 / (vL + mu_in * r_bar_sq_cached)
      end function

    end subroutine get_mu_root_leptonic

    !--------------------------------------------------------------------------
    ! Evaluate Kastaun primitives at given mu (Eq. 43-48 of Kastaun+2021)
    !--------------------------------------------------------------------------
    subroutine eval_kastaun_primitives(mu_in, q_in, r2_in, b2_in, rdotb_in, &
        D_in, ye_in, ymu_in, rho_out, eps_out, W_out, prs_out, vi_out)
      double precision, intent(in)  :: mu_in, q_in, r2_in, b2_in, rdotb_in, D_in
      double precision, intent(in)  :: ye_in, ymu_in
      double precision, intent(out) :: rho_out, eps_out, W_out, prs_out, vi_out(ndir)

      double precision :: chi_f, r_bar_sq, vsq, v0_sq_local
      double precision :: eps_rng(2), tmp_rho

      chi_f     = 1.0d0 / (1.0d0 + mu_in * b2_in)
      r_bar_sq  = r2_in * chi_f**2 + mu_in * chi_f * (1.0d0 + chi_f) * rdotb_in**2
      r_bar_sq_cached = r_bar_sq   ! store for func_kastaun_f

      v0_sq_local = r2_in / (h0**2 + r2_in)
      vsq   = min(mu_in**2 * r_bar_sq, v0_sq_local)
      W_out = 1.0d0 / sqrt(max(1.0d0 - vsq, 1.0d-15))

      rho_out = D_in / W_out
      rho_out = max(lep_eos_rhomin, min(lep_eos_rhomax, rho_out))

      eps_out = W_out * (q_in - mu_in * r_bar_sq) &
               + vsq * W_out**2 / (1.0d0 + W_out)

      ! Bound eps
      call leptonic_get_eps_range(rho_out, eps_rng(1), eps_rng(2), ye_in, ymu_in)
      eps_out = max(eps_rng(1), min(eps_rng(2), eps_out))

      ! Pressure call with leptonic EOS
      call leptonic_get_pressure_one_grid(prs_out, rho_out, eps_out, &
          ye=ye_in, ymu=ymu_in)

      vi_out = 0.0d0  ! velocities computed later from stilde^i / (hW)
    end subroutine eval_kastaun_primitives

    double precision :: r_bar_sq_cached   ! shared between eval and func

  end subroutine leptonic_con2prim

end module mod_imhd_con2prim_leptonic
