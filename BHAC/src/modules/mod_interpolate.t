!================================================================================
!
!    BHAC (The Black Hole Accretion Code) solves the equations of
!    general relativistic magnetohydrodynamics and other hyperbolic systems
!    in curved spacetimes.
!
!    Copyright (C) 2019 Oliver Porth, Hector Olivares, Yosuke Mizuno, Ziri Younsi,
!    Luciano Rezzolla, Elias Most, Bart Ripperda and Fabio Bacchini
!
!    This file is part of BHAC.
!
!    BHAC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BHAC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BHAC.  If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

!=======================================================================
module mod_interpolate
  !
  ! Contains routines to find points and to interpolate to them.
  !
  implicit none

  !=======================================================================
contains
  !=======================================================================
{#IFNDEF UNIT_TESTS
  logical function point_in_domain(x)

    include 'amrvacdef.f'
    double precision, dimension(ndim), intent(in)  :: x
    integer                                        :: idim
    !----------------------------------------------------------------------

    point_in_domain = .true.

    do idim=1,ndim
       select case(idim)
          {case (^D)
          if (x(^D) .lt. xprobmin^D) then
             point_in_domain = .false.
             exit
          end if
          if (x(^D) .ge. xprobmax^D) then
             point_in_domain = .false.
             exit
          end if
          \}
       end select
    end do

  end function point_in_domain
  !====================================================================
  subroutine find_point_ipe(x,igrid_point,ipe_point)

    use mod_forest, only: tree_node_ptr, tree_root
    use mod_slice, only: get_igslice
    include 'amrvacdef.f'

    double precision, dimension(ndim), intent(in)   :: x
    integer, intent(out)                            :: igrid_point, ipe_point

    integer, dimension(ndir,nlevelshi)              :: ig
    integer                                         :: idim, ic(ndim)
    type(tree_node_ptr)                             :: branch
    !--------------------------------------------------------------------

    ! first check if the point is in the domain
    if (.not. point_in_domain(x)) then
       igrid_point = -1
       ipe_point   = -1
       return
    end if

    ! get the index on each level
    do idim = 1, ndim
       call get_igslice(idim,x(idim),ig(idim,:))
    end do

    ! traverse the tree until leaf is found
    branch=tree_root({ig(^D,1)})
    do while (.not.branch%node%leaf)
       {ic(^D)=ig(^D,branch%node%level+1) - 2 * branch%node%ig^D +2\}
       branch%node => branch%node%child({ic(^D)})%node
    end do

    igrid_point = branch%node%igrid
    ipe_point   = branch%node%ipe

  end subroutine find_point_ipe
  !===========================================================================
  subroutine interpolate_var(igrid,ixI^L,gf,x,xloc,gfloc)

    include 'amrvacdef.f'
    integer, intent(in)                   :: igrid,ixI^L
    double precision, intent(in)          :: gf(ixI^S)
    double precision, intent(in)          :: x(ixI^S,1:ndim)
    double precision, intent(in)          :: xloc(1:ndim)
    double precision, intent(out)         :: gfloc
    integer                               :: ic^D, ic1^D, ic2^D, idir
    double precision                      :: xd^D
    {#IFDEF D2
    double precision                      :: c00, c10
    }
    {#IFDEF D3
    double precision                      :: c0, c1, c00, c10, c01, c11
    }
    character(len=1024)                   :: line
    !---------------------------------------------------------------------------

    ! flat interpolation:
    {ic^D = int((xloc(^D)-rnode(rpxmin^D_,igrid))/rnode(rpdx^D_,igrid)) + 1 + dixB \}
    !gfloc = gf(ic^D)



    if ({ic^D.lt.ixImin^D .or. ic^D.gt.ixImax^D|.or.}) then
       line = ''
       write(line,"(a)") 'Trying to flat-interpolate from out of grid!'
       write(line,"(a,a,^NDes14.6)") trim(line),' position: ',xloc(1:ndim)
       write(line,"(a,a,i4.3)") trim(line),' index: ', ic^D
       {^D&
       if (ic^D.lt.ixImin^D .or. ic^D.gt.ixImax^D) then
          write(line,"(a,a,i3.2)") trim(line),' Bounds exceeded in direction: ',^D
       else
          write(line,"(a,a,i3.2)") trim(line),' Bounds OK in direction: ',^D
       end if
       write(line,"(a,a,2es14.6)") trim(line),' grid range: ', x(ixImin^DD,^D), x(ixImax^DD,^D)
       \}
       call mpistop(line)
    end if
    
    ! linear interpolation:
    {
    if (x({ic^DD},^D) .lt. xloc(^D)) then
       ic1^D = ic^D
    else
       ic1^D = ic^D -1
    end if
    ic2^D = ic1^D + 1
    \}
    
    {^D&
    if (ic1^D.lt.ixImin^D .or. ic2^D.gt.ixImax^D) then
       line = ''
       write(line,"(a)") 'Trying to interpolate from out of grid!'
       write(line,"(a,a,i3.2)") trim(line),' direction: ',^D
       write(line,"(a,a,^NDes14.6)") trim(line),' position: ',xloc(1:ndim)
       write(line,"(a,a,2i4.3)") trim(line),' indices: ', ic1^D,ic2^D
       write(line,"(a,a,2es14.6)") trim(line),' grid range: ', x(ixImin^DD,^D), x(ixImax^DD,^D)
       call mpistop(line)
    end if
    \}
    
    
    {#IFDEF D1
    xd1 = (xloc(1)-x(ic11,1)) / (x(ic21,1) - x(ic11,1))
    gfloc  = gf(ic11) * (1.0d0 - xd1) + gf(ic21) * xd1
    }
    {#IFDEF D2
    xd1 = (xloc(1)-x(ic11,ic12,1)) / (x(ic21,ic12,1) - x(ic11,ic12,1))
    xd2 = (xloc(2)-x(ic11,ic12,2)) / (x(ic11,ic22,2) - x(ic11,ic12,2))
    c00 = gf(ic11,ic12) * (1.0d0 - xd1) + gf(ic21,ic12) * xd1
    c10 = gf(ic11,ic22) * (1.0d0 - xd1) + gf(ic21,ic22) * xd1
    gfloc  = c00 * (1.0d0 - xd2) + c10 * xd2
    }
    {#IFDEF D3
    xd1 = (xloc(1)-x(ic11,ic12,ic13,1)) / (x(ic21,ic12,ic13,1) - x(ic11,ic12,ic13,1))
    xd2 = (xloc(2)-x(ic11,ic12,ic13,2)) / (x(ic11,ic22,ic13,2) - x(ic11,ic12,ic13,2))
    xd3 = (xloc(3)-x(ic11,ic12,ic13,3)) / (x(ic11,ic12,ic23,3) - x(ic11,ic12,ic13,3))

    c00 = gf(ic11,ic12,ic13) * (1.0d0 - xd1) + gf(ic21,ic12,ic13) * xd1
    c10 = gf(ic11,ic22,ic13) * (1.0d0 - xd1) + gf(ic21,ic22,ic13) * xd1
    c01 = gf(ic11,ic12,ic23) * (1.0d0 - xd1) + gf(ic21,ic12,ic23) * xd1
    c11 = gf(ic11,ic22,ic23) * (1.0d0 - xd1) + gf(ic21,ic22,ic23) * xd1

    c0  = c00 * (1.0d0 - xd2) + c10 * xd2
    c1  = c01 * (1.0d0 - xd2) + c11 * xd2

    gfloc = c0 * (1.0d0 - xd3) + c1 * xd3
    }

  end subroutine interpolate_var
  !=======================================================================

  subroutine lagrange_interpolation(xx,yy,x,y,inverse)
    ! Interpolate xx=[x1,x2,...xn]
    !             yy=[f(x1),f(x2),...,f(xn)]
    ! Using Lagrange polynomial P(x)
    ! Returns y = P(x)
    implicit none
    double precision,dimension(:),intent(in) :: xx, yy
    double precision,intent(in)              :: x
    double precision,intent(out)             :: y
    logical,intent(in),optional              :: inverse

    !local variables:
    logical :: inv
    integer :: j,k,n,m, i_start, i_end, step
    double precision :: p

    inv=.false.
    if (present(inverse)) inv=inverse

    !check number of points:
    n = size(xx); m = size(yy)
    if ( (n/=m).or.(n<2) ) error stop &
        'Error: vectors must be the same size.'

    if (inv) then
      i_start=n; i_end=1; step = -1
    else
      i_start=1; i_end=n; step = 1
    end if

    !sum each of the Pj(x) terms:
    y = 0.0d0
    do j=i_start,i_end,step
        !compute Pj(x):
        p = yy(j)
        do k=i_start,i_end,step
            if (k/=j) p = p * (x-xx(k)) / (xx(j)-xx(k))
        end do
        y = y + p
    end do
   end subroutine lagrange_interpolation


!  subroutine lagrange_interpolation(xx,yy,x,y)
!      !
!      ! Interpolate xx=[x1,x2,...xn]
!      !             yy=[f(x1),f(x2),...,f(xn)]
!      ! Using Lagrange polynomial P(x)
!      ! Returns y = P(x)
!      !
!      implicit none
!
!      double precision,dimension(:),intent(in) :: xx
!      double precision,dimension(:),intent(in) :: yy
!      double precision,intent(in)              :: x
!      double precision,intent(out)             :: y
!
!      !local variables:
!      integer :: j,k,n,m
!      double precision :: p
!
!      !check number of points:
!      n = size(xx)
!      m = size(yy)
!      if ( (n/=m).or.(n<2) ) error stop &
!          'Error: vectors must be the same size.'
!
!      !sum each of the Pj(x) terms:
!      y = 0.0d0
!      do j=1,n
!          !compute Pj(x):
!          p = yy(j)
!          do k=1,n
!              if (k/=j) p = p * (x-xx(k)) / (xx(j)-xx(k))
!          end do
!          y = y + p
!      end do
!
!   end subroutine lagrange_interpolation

!fixme: fix for more generic, not only for cfc
   subroutine metric_interpolation(ixI^L,ixO^L,idims,w,x,wi,xi)
     include 'amrvacdef.f'
   
     integer, intent(in) :: ixI^L, ixO^L, idims
     double precision, dimension(ixI^S,1:nw), intent(in) :: w
     double precision, dimension(ixI^S,1:nw), intent(out) :: wi
     double precision, dimension(ixI^S,1:ndim), intent(in) :: x, xi
   
     ! local vars
     integer :: ix^D
     integer :: iw
     integer :: n_lo, n_hi
     double precision, allocatable :: xx(:), yy(:)
  
     {#IFNDEF DY_SP 
        call mpistop('You do not use Dynamical spacetime, &
                      why you need metric_interpolation during evolution')
     }
   
     if ( mod(fv_n_interp,2) == 0 ) then
        n_lo = fv_n_interp / 2 - 1
        n_hi = fv_n_interp / 2
     else
        n_lo = (fv_n_interp+1) / 2 - 1
        n_hi = (fv_n_interp-1) / 2
     endif
     allocate(xx(fv_n_interp))
     allocate(yy(fv_n_interp))
   
     do iw = nmetric_lo, nmetric_hi
   
        {^IFONED
        do ix1 = ixOmin1, ixOmax1
           xx = x(ix1-n_lo:ix1+n_hi, 1)
           yy = w(ix1-n_lo:ix1+n_hi, iw)
           call lagrange_interpolation(xx,yy,&
                   !xi(ix1,1),wi(ix1, iw))
                   xi(ix1,1),wi(ix1, iw),(xi(ix^D, idims)<0.0d0))
        end do
        }
        {^IFTWOD
        select case (idims)
        case (1)
           do ix2 = ixOmin2, ixOmax2
              do ix1 = ixOmin1, ixOmax1
                 xx = x(ix1-n_lo:ix1+n_hi, ix2, 1)
                 yy = w(ix1-n_lo:ix1+n_hi, ix2, iw)
                 call lagrange_interpolation(xx,yy,&
                         !xi(ix1, ix2, 1), wi(ix1, ix2, iw))
                         xi(ix1, ix2, 1), wi(ix1, ix2, iw),(xi(ix^D, idims)<0.0d0))
              end do
           end do
        case (2)
           do ix1 = ixOmin1, ixOmax1
              do ix2 = ixOmin2, ixOmax2
                 xx = x(ix1, ix2-n_lo:ix2+n_hi, 2)
                 yy = w(ix1, ix2-n_lo:ix2+n_hi, iw)
                 call lagrange_interpolation(xx,yy,&
                         !xi(ix1, ix2, 2), wi(ix1, ix2, iw))
                         xi(ix1, ix2, 2), wi(ix1, ix2, iw),(xi(ix^D, idims)<0.0d0))
              end do
           end do
        case default
           call mpistop(" idims can only be 1 or 2 in 2D")
        end select
        }
        {^IFTHREED
        select case (idims)
        case (1)
           do ix3 = ixOmin3, ixOmax3
              do ix2 = ixOmin2, ixOmax2
                 do ix1 = ixOmin1, ixOmax1
                    xx = x(ix1-n_lo:ix1+n_hi, ix2, ix3, 1)
                    yy = w(ix1-n_lo:ix1+n_hi, ix2, ix3, iw)
                    call lagrange_interpolation(xx,yy,&
                            !xi(ix^D, 1), wi(ix^D, iw))
                            xi(ix^D, 1), wi(ix^D, iw),(xi(ix^D, idims)<0.0d0))
                 end do
              end do
           end do
        case (2)
           do ix3 = ixOmin3, ixOmax3
              do ix1 = ixOmin1, ixOmax1
                 do ix2 = ixOmin2, ixOmax2
                    xx = x(ix1, ix2-n_lo:ix2+n_hi, ix3, 2)
                    yy = w(ix1, ix2-n_lo:ix2+n_hi, ix3, iw)
                    call lagrange_interpolation(xx,yy,&
                            !xi(ix^D, 2), wi(ix^D, iw))
                            xi(ix^D, 2), wi(ix^D, iw),(xi(ix^D, idims)<0.0d0))
                 end do
              end do
           end do
        case (3)
           do ix1 = ixOmin1, ixOmax1
              do ix2 = ixOmin2, ixOmax2
                 do ix3 = ixOmin3, ixOmax3
                    xx = x(ix1, ix2, ix3-n_lo:ix3+n_hi, 3)
                    yy = w(ix1, ix2, ix3-n_lo:ix3+n_hi, iw)
                    call lagrange_interpolation(xx,yy,&
                            !xi(ix^D, 3), wi(ix^D, iw))
                            xi(ix^D, 3), wi(ix^D, iw),(xi(ix^D, idims)<0.0d0))
                 end do
              end do
           end do
        case default
           call mpistop(" idims can only be 1, 2 or 3 in 3D")
        end select
        }
   
     end do
   
     deallocate(xx)
     deallocate(yy)
   end subroutine metric_interpolation
 }
   !==================================================================

   subroutine TrilinearInterpolation(x, y, z, x1, x2, y1, y2, z1, z2, f111, f211, f121, f221, f112, f212, f122, f222, result)
      ! Inputs:
      ! x, y, z   : The coordinates of the point where we want to interpolate
      ! x1, x2    : The x-coordinates of the grid points
      ! y1, y2    : The y-coordinates of the grid points
      ! z1, z2    : The z-coordinates of the grid points
      ! f111, f211, f121, f221, f112, f212, f122, f222 : Function values at the 8 corners of the cube
      ! Output:
      ! result    : Interpolated value at the point (x, y, z)     
      double precision, intent(in)  :: x, y, z      ! Coordinates of the interpolation point
      double precision, intent(in)  :: x1, x2       ! x-coordinates of the grid points
      double precision, intent(in)  :: y1, y2       ! y-coordinates of the grid points
      double precision, intent(in)  :: z1, z2       ! z-coordinates of the grid points
      double precision, intent(in)  :: f111, f211   ! Function values at the corners of the grid
      double precision, intent(in)  :: f121, f221
      double precision, intent(in)  :: f112, f212
      double precision, intent(in)  :: f122, f222
      double precision, intent(out)  :: result      ! Interpolated result
  
      double precision :: xd, yd, zd
      double precision :: c00, c10, c01, c11, c0, c1
  
      ! Calculate the normalized distances between the points
      xd = (x - x1) / (x2 - x1)
      yd = (y - y1) / (y2 - y1)
      zd = (z - z1) / (z2 - z1)
  
      ! Interpolate along the x-axis for the 4 z-constant planes
      c00 = f111 * (1.0 - xd) + f211 * xd
      c10 = f121 * (1.0 - xd) + f221 * xd
      c01 = f112 * (1.0 - xd) + f212 * xd
      c11 = f122 * (1.0 - xd) + f222 * xd
  
      ! Interpolate along the y-axis
      c0 = c00 * (1.0 - yd) + c10 * yd
      c1 = c01 * (1.0 - yd) + c11 * yd
  
      ! Interpolate along the z-axis
      result = c0 * (1.0 - zd) + c1 * zd
  
    end subroutine TrilinearInterpolation

   subroutine QuadrilinearInterpolation(x, y, z, q, x1, x2, y1, y2, z1, z2, q1, q2, &
        f1111, f2111, f1211, f2211, f1121, f2121, f1221, f2221, &
        f1112, f2112, f1212, f2212, f1122, f2122, f1222, f2222, result)
      double precision, intent(in)  :: x, y, z, q
      double precision, intent(in)  :: x1, x2, y1, y2, z1, z2, q1, q2
      double precision, intent(in)  :: f1111, f2111, f1211, f2211
      double precision, intent(in)  :: f1121, f2121, f1221, f2221
      double precision, intent(in)  :: f1112, f2112, f1212, f2212
      double precision, intent(in)  :: f1122, f2122, f1222, f2222
      double precision, intent(out) :: result

      double precision :: qd
      double precision :: result_q1, result_q2

      call TrilinearInterpolation(x, y, z, x1, x2, y1, y2, z1, z2, &
           f1111, f2111, f1211, f2211, f1121, f2121, f1221, f2221, result_q1)
      call TrilinearInterpolation(x, y, z, x1, x2, y1, y2, z1, z2, &
           f1112, f2112, f1212, f2212, f1122, f2122, f1222, f2222, result_q2)

      qd = (q - q1) / (q2 - q1)
      result = result_q1 * (1.0d0 - qd) + result_q2 * qd

   end subroutine QuadrilinearInterpolation

end module mod_interpolate
!=======================================================================
