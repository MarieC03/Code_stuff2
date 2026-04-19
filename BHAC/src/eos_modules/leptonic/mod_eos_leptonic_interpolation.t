!================================================================================
!
!  mod_eos_leptonic_interpolation.t
!
!  Interpolation helper for the 4D leptonic EOS in BHAC.
!
!  Provides:
!    intep3d_lep_many   – vectorised trilinear interpolation over nvars
!    intep3d_lep        – single-variable trilinear interpolation
!
!  The routines are identical in algorithm to mod_eos_interpolation but
!  operate on the leptonic table arrays and support a log-scale third axis
!  (needed for the muonic table, where ymu is stored as ln(ymu)).
!
!  Copyright (C) 2024  Harry Ho-Yin Ng
!  Based on BHAC mod_eos_interpolation.
!
!================================================================================

module mod_eos_leptonic_interpolation
  implicit none
  public

contains

  !----------------------------------------------------------------------------
  !> Trilinear interpolation of nvars quantities simultaneously.
  !> All three coordinate axes are given as uniform-spacing arrays (in the
  !> sense used by intep3d_many in BHAC: spacing determined from first two
  !> elements, so axes need NOT be exactly uniform but spacing is approximated
  !> as constant from element [1] to [2]).
  !>
  !> ft  : table array of shape (nx, ny, nz, nvars)
  !> x   : ln(rho)   [code units]
  !> y   : ln(T)     [MeV]
  !> z   : third coordinate (yle linear / ln(ymu) for muon table / yp linear)
  !----------------------------------------------------------------------------
  subroutine intep3d_lep_many(x, y, z, f, ft, nx, ny, nz, nvars, xt, yt, zt)
    implicit none

    integer,          intent(in)  :: nx, ny, nz, nvars
    double precision, intent(in)  :: x, y, z
    double precision, intent(in)  :: xt(nx), yt(ny), zt(nz)
    double precision, intent(in)  :: ft(nx, ny, nz, nvars)
    double precision, intent(out) :: f(nvars)

    double precision :: dxi, dyi, dzi, xx, yy, zz
    integer          :: ix, iy, iz

    ! ---- locate cell index (equidistant spacing assumed) -------------------
    ix = max(1, min(int((x - xt(1)) / (xt(2) - xt(1))) + 1, nx-1))
    iy = max(1, min(int((y - yt(1)) / (yt(2) - yt(1))) + 1, ny-1))
    iz = max(1, min(int((z - zt(1)) / (zt(2) - zt(1))) + 1, nz-1))

    ! ---- fractional coordinates --------------------------------------------
    dxi = 1.0d0 / (xt(ix+1) - xt(ix))
    dyi = 1.0d0 / (yt(iy+1) - yt(iy))
    dzi = 1.0d0 / (zt(iz+1) - zt(iz))

    xx  = (x - xt(ix)) * dxi
    yy  = (y - yt(iy)) * dyi
    zz  = (z - zt(iz)) * dzi

    ! ---- trilinear interpolation -------------------------------------------
    f(1:nvars) = &
        ft(ix  , iy  , iz  , 1:nvars) * (1.0d0-xx) * (1.0d0-yy) * (1.0d0-zz) &
      + ft(ix+1, iy  , iz  , 1:nvars) * xx          * (1.0d0-yy) * (1.0d0-zz) &
      + ft(ix  , iy+1, iz  , 1:nvars) * (1.0d0-xx) * yy          * (1.0d0-zz) &
      + ft(ix  , iy  , iz+1, 1:nvars) * (1.0d0-xx) * (1.0d0-yy) * zz          &
      + ft(ix+1, iy+1, iz  , 1:nvars) * xx          * yy          * (1.0d0-zz) &
      + ft(ix+1, iy  , iz+1, 1:nvars) * xx          * (1.0d0-yy) * zz          &
      + ft(ix  , iy+1, iz+1, 1:nvars) * (1.0d0-xx) * yy          * zz          &
      + ft(ix+1, iy+1, iz+1, 1:nvars) * xx          * yy          * zz

  end subroutine intep3d_lep_many

  !----------------------------------------------------------------------------
  !> Single-variable trilinear interpolation (convenience wrapper)
  !----------------------------------------------------------------------------
  subroutine intep3d_lep(x, y, z, f, ft, nx, ny, nz, xt, yt, zt)
    implicit none

    integer,          intent(in)  :: nx, ny, nz
    double precision, intent(in)  :: x, y, z
    double precision, intent(in)  :: xt(nx), yt(ny), zt(nz)
    double precision, intent(in)  :: ft(nx, ny, nz)
    double precision, intent(out) :: f

    double precision :: ftmp(nx, ny, nz, 1), fout(1)

    ftmp(:,:,:,1) = ft(:,:,:)
    call intep3d_lep_many(x, y, z, fout, ftmp, nx, ny, nz, 1, xt, yt, zt)
    f = fout(1)
  end subroutine intep3d_lep

end module mod_eos_leptonic_interpolation
