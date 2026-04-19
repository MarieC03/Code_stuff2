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

!=============================================================================
module mod_rtsafe2D

  implicit none

  integer, save                         :: MAXIT=1000
  !$OMP THREADPRIVATE(MAXIT)
  
contains
  !-----------------------------------------------------------------------------

  !=============================================================================
  subroutine rt1D(func,x,absacc,ierror)

    double precision, intent(in)            :: absacc ! accuracy
    double precision, intent(inout)         :: x ! initial guess and final solution
    integer, intent(out)                    :: ierror

    interface

       subroutine func(s,f,df)
         double precision, intent(in)      :: s
         double precision, intent(out)     :: f, df
       end subroutine func

    end interface

    ! .. local ..
    integer                                :: it
    double precision                       :: f, df
    double precision                       :: dx, err
    !-----------------------------------------------------------------------------
    ierror = 1

    do it=1,MAXIT

       call func(x,f,df)
       
       if (df .eq. 0.0d0) then
          ierror = 2
          return
       else
          dx = - f/df
          x = x + dx
          if (abs(dx) .lt. absacc) then
             ierror = 0
             return
          end if
       end if
    end do
    
  end subroutine rt1D
  !=============================================================================
  subroutine rtsafe2D(f1,f2,xl,xh,yl,yh,absacc,xroot,yroot,niter,ierror,validate)

    double precision, intent(in)             :: xl, xh, yl, yh, absacc
    double precision, intent(inout)          :: xroot, yroot
    integer, intent(out)                     :: niter, ierror

    interface

       subroutine f1(x,y,f,dfdx,dfdy,return_derivatives)
         double precision, intent(in)        :: x, y
         double precision, intent(out)       :: f, dfdx, dfdy
         logical, intent(in)                 :: return_derivatives
       end subroutine f1

       subroutine f2(x,y,f,dfdx,dfdy,return_derivatives)
         double precision, intent(in)        :: x, y
         double precision, intent(out)       :: f, dfdx, dfdy
         logical, intent(in)                 :: return_derivatives
       end subroutine f2

       subroutine validate(xl,xh,yl,yh,xroot,yroot,xroot_old,yroot_old)
         double precision, intent(in)        :: xl, xh, yl, yh, xroot_old, yroot_old
         double precision, intent(inout)     :: xroot, yroot
       end subroutine validate
       
    end interface

    optional                                 :: validate
    
    ! .. local ..
    double precision                         :: myf1, myf2, df1dx, df1dy, df2dx, df2dy
    double precision                         :: dxold, dyold, dx, dy, detJ
    double precision                         :: tmpx, tmpy, xroot_old, yroot_old, errx
    integer                                  :: it, iextra
    integer, parameter                       :: nextra = 1
    !-----------------------------------------------------------------------------

    ierror = 0
    niter  = 0
    iextra = 0

    ! ------------------------------
    ! Iterate:
    ! ------------------------------
    do it = 1, MAXIT

       xroot_old = xroot; yroot_old = yroot
       call f1(xroot,yroot,myf1,df1dx,df1dy,.true.)
       call f2(xroot,yroot,myf2,df2dx,df2dy,.true.)

       
       ! Get the next step:
       detJ = df1dx*df2dy - df1dy*df2dx
       dx = (df1dy*myf2 - df2dy*myf1) / detJ
       dy = (df2dx*myf1 - df1dx*myf2) / detJ

       ! ------------------------------
       ! x-Direction:
       ! ------------------------------
       tmpx = xroot
       xroot = xroot + dx


       ! ------------------------------
       ! y-Direction:
       ! ------------------------------
       tmpy = yroot
       yroot = yroot + dy

       
       ! ------------------------------
       ! Check if we have reached the accuracy goal:
       ! only in the first variable
       ! ------------------------------
       if (yroot .eq. 0.0d0) then
          errx = abs(dy)
       else
          errx = abs(dy/yroot)
       end if
       if ((errx .lt. absacc) &
            .or. (tmpx .eq. xroot .and. tmpy .eq. yroot) ) then
          iextra = iextra + 1
       end if
       
       ! ------------------------------
       ! Correct within ranges (this is problem dependent!):
       ! ------------------------------
       if (present(validate)) then
          call validate(xl,xh,yl,yh,xroot,yroot,xroot_old,yroot_old)
       end if       

       niter = niter + 1

       ! ------------------------------
       ! Get out after the extra-steps:
       ! ------------------------------
       if (iextra .gt. nextra) exit

    end do

    if (niter .ge. maxit) &
         ierror = 3 ! Exceeded maximum number of iterations

  end subroutine rtsafe2D
  !=============================================================================
  
end module mod_rtsafe2D
!=============================================================================
