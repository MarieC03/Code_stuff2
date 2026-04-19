!> Moyt(2:ny-1)dule for eos
module mod_eos_interpolation

  implicit none
  public

  contains

  subroutine intep3d(x, y, z, f, &
                  ft, &
                  nx, ny, nz, &
                  xt, yt, zt, &
                  ix_out, iy_out, iz_out, &
                  d1, d2, d3)
     implicit none
     integer, intent(in)   ::  nx, ny, nz
     double precision, intent(in)                   ::  x, y, z
     double precision, intent(out)                  ::  f
     integer, intent(out), optional                 ::  ix_out, iy_out, iz_out
     double precision, intent(out), optional        ::  d1, d2, d3
     double precision, intent(in)                   ::  xt(nx), yt(ny), zt(nz)
     double precision, intent(in)                   ::  ft(nx, ny, nz)

     integer                                        ::  ix, iy, iz
     double precision                               ::  dxi, dyi, dzi
     double precision                               ::  xx, yy, zz

!------  determine spacing parameters of (equidistant!!!) table
     ix = max(0, min(int((x - xt(1)) / (xt(2) - xt(1))) + 1, nx-1))
     iy = max(0, min(int((y - yt(1)) / (yt(2) - yt(1))) + 1, ny-1))
     iz = max(0, min(int((z - zt(1)) / (zt(2) - zt(1))) + 1, nz-1))

     if(present(ix_out)) ix_out = ix
     if(present(iy_out)) iy_out = iy
     if(present(iz_out)) iz_out = iz

     ! determine spacing parameters of (equidistant!!!) table

     dxi   = 1.d0 / ( xt(ix+1)-xt(ix) )
     dyi   = 1.d0 / ( yt(iy+1)-yt(iy) )
     dzi   = 1.d0 / ( zt(iz+1)-zt(iz) )

     xx = ( x - xt(ix) ) * dxi
     yy = ( y - yt(iy) ) * dyi
     zz = ( z - zt(iz) ) * dzi

     f = ft(ix  , iy  , iz  ) * (1.d0 - xx) * (1.d0 - yy) * (1.d0 - zz) &
       + ft(ix+1, iy  , iz  ) * xx          * (1.d0 - yy) * (1.d0 - zz) &
       + ft(ix  , iy+1, iz  ) * (1.d0 - xx) * yy          * (1.d0 - zz) &
       + ft(ix  , iy  , iz+1) * (1.d0 - xx) * (1.d0 - yy) * zz          &
       + ft(ix+1, iy+1, iz  ) * xx          * yy          * (1.d0 - zz) &
       + ft(ix+1, iy  , iz+1) * xx          * (1.d0 - yy) * zz          &
       + ft(ix  , iy+1, iz+1) * (1.d0 - xx) * yy          * zz          &
       + ft(ix+1, iy+1, iz+1) * xx          * yy          * zz

     if(present(d1)) d1 = &
         ft(ix  , iy  , iz  ) * ( - dxi   ) * (1.d0 - yy) * (1.d0 - zz) &
       + ft(ix+1, iy  , iz  ) * dxi         * (1.d0 - yy) * (1.d0 - zz) &
       + ft(ix  , iy+1, iz  ) * ( - dxi   ) * yy          * (1.d0 - zz) &
       + ft(ix  , iy  , iz+1) * ( - dxi   ) * (1.d0 - yy) * zz          &
       + ft(ix+1, iy+1, iz  ) * dxi         * yy          * (1.d0 - zz) &
       + ft(ix+1, iy  , iz+1) * dxi         * (1.d0 - yy) * zz          &
       + ft(ix  , iy+1, iz+1) * ( - dxi   ) * yy          * zz          &
       + ft(ix+1, iy+1, iz+1) * dxi         * yy          * zz

     if(present(d2)) d2 = &
         ft(ix  , iy  , iz  ) * (1.d0 - xx) * ( - dyi   ) * (1.d0 - zz) &
       + ft(ix+1, iy  , iz  ) * xx          * ( - dyi   ) * (1.d0 - zz) &
       + ft(ix  , iy+1, iz  ) * (1.d0 - xx) * dyi         * (1.d0 - zz) &
       + ft(ix  , iy  , iz+1) * (1.d0 - xx) * ( - dyi   ) * zz          &
       + ft(ix+1, iy+1, iz  ) * xx          * dyi         * (1.d0 - zz) &
       + ft(ix+1, iy  , iz+1) * xx          * ( - dyi   ) * zz          &
       + ft(ix  , iy+1, iz+1) * (1.d0 - xx) * dyi         * zz          &
       + ft(ix+1, iy+1, iz+1) * xx          * dyi         * zz

     if(present(d3)) d3 = &
         ft(ix  , iy  , iz  ) * (1.d0 - xx) * (1.d0 - yy) * ( - dzi   ) &
       + ft(ix+1, iy  , iz  ) * xx          * (1.d0 - yy) * ( - dzi   ) &
       + ft(ix  , iy+1, iz  ) * (1.d0 - xx) * yy          * ( - dzi   ) &
       + ft(ix  , iy  , iz+1) * (1.d0 - xx) * (1.d0 - yy) * dzi         &
       + ft(ix+1, iy+1, iz  ) * xx          * yy          * ( - dzi   ) &
       + ft(ix+1, iy  , iz+1) * xx          * (1.d0 - yy) * dzi         &
       + ft(ix  , iy+1, iz+1) * (1.d0 - xx) * yy          * dzi         &
       + ft(ix+1, iy+1, iz+1) * xx          * yy          * zz

  end subroutine

  subroutine intep3d_many( x, y, z, f, ft, nx, ny, nz, nvars, xt, yt, zt)

      implicit none

!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  8-point Lagrange linear interpolation formula
!
!     x        input vector of first  variable
!     y        input vector of second variable
!     z        input vector of third  variable
!
!     f        output vector of interpolated function values
!
!     kt       vector length of input and output vectors
!
!     ft       3d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     nz       z-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!     zt       vector of z-coordinates of table
!
!---------------------------------------------------------------------

     integer, intent(in)            :: nx,ny,nz,nvars
     double precision, intent(in)   :: x, y, z, xt(nx), yt(ny), zt(nz), ft(nx,ny,nz,nvars)
     double precision, intent(out)  :: f(nvars)

     double precision dxi,dyi,dzi, xx, yy, zz
     integer ix,iy,iz

     ! determine spacing parameters of (equidistant!!!) table

     ix = max(0, min(int((x - xt(1)) / (xt(2) - xt(1))) + 1, nx-1))
     iy = max(0, min(int((y - yt(1)) / (yt(2) - yt(1))) + 1, ny-1))
     iz = max(0, min(int((z - zt(1)) / (zt(2) - zt(1))) + 1, nz-1))

     dxi   = 1.d0 / ( xt(ix+1)-xt(ix) )
     dyi   = 1.d0 / ( yt(iy+1)-yt(iy) )
     dzi   = 1.d0 / ( zt(iz+1)-zt(iz) )

     xx  = (x - xt(ix)) * dxi
     yy  = (y - yt(iy)) * dyi
     zz  = (z - zt(iz)) * dzi

     f(1:nvars) = &
         ft(ix  , iy  , iz  , 1:nvars) * (1.d0 - xx) * (1.d0 - yy) * (1.d0 - zz) &
       + ft(ix+1, iy  , iz  , 1:nvars) * xx          * (1.d0 - yy) * (1.d0 - zz) &
       + ft(ix  , iy+1, iz  , 1:nvars) * (1.d0 - xx) * yy          * (1.d0 - zz) &
       + ft(ix  , iy  , iz+1, 1:nvars) * (1.d0 - xx) * (1.d0 - yy) * zz          &
       + ft(ix+1, iy+1, iz  , 1:nvars) * xx          * yy          * (1.d0 - zz) &
       + ft(ix+1, iy  , iz+1, 1:nvars) * xx          * (1.d0 - yy) * zz          &
       + ft(ix  , iy+1, iz+1, 1:nvars) * (1.d0 - xx) * yy          * zz          &
       + ft(ix+1, iy+1, iz+1, 1:nvars) * xx          * yy          * zz

  end subroutine intep3d_many





end module mod_eos_interpolation

