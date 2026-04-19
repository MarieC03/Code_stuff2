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
module mod_ngmres

  use mod_lu

! This module packages the nonlinear GMRES for solving systems of nonlinear equations
! Inspired by C.T. Kelley's Matlab implementation of nsolgm.m
! Fabio Bacchini, October 2018

contains
!=============================================================================
subroutine nsolgm(x, sol, n, ierr, ferr, abs_tol, rel_tol, etamax, maxit, lmaxit, residual)

  ! Nonlinear solver - Inexact Newton iteration
  ! Fabio Bacchini, October 2018

  implicit none

  ! Input variables
  integer,intent(in)                        :: n ! Dimension of the nonlinear system
  double precision,intent(in),dimension(n)  :: x ! Initial guess vector
  double precision,intent(out),dimension(n) :: sol ! Solution vector
  double precision,intent(in)               :: abs_tol, rel_tol, etamax ! Tolerances
  double precision,intent(out)              :: ferr ! Final error
  integer,intent(in)                        :: maxit, lmaxit ! Outer and inner maximum iterations
  integer,intent(inout)                     :: ierr ! Final error code
  ! Local variables
  double precision,dimension(n)             :: f0, step, xk
  double precision                          :: r, normf, normfold
  double precision                          :: stop_tol, eta, etaold, etanew, gg=0.9d0
  integer                                   :: ii
  external residual

  ! Initialise solver
  xk = x
  call residual(xk, f0, n) ! Initial residual
  normf = norm2(f0)/dsqrt(dble(n))
  eta = etamax

  ! Main loop
  stop_tol = abs_tol + rel_tol*norm2(f0)
  r = 1.d0
  ii = 0
  
  ! Outer nonlinear iteration    
  do while((normf .gt. stop_tol) .and. (ii .lt. maxit))

    ii = ii + 1
    normfold = normf

    ! Call linear GMRES for inner iteration: determine step length
    call fdgmres(f0, xk, step, n, eta, lmaxit, residual)
    
    ! Update solution
    xk = xk + step

    ! Recalculate errors
    call residual(xk, f0, n)
    normf = norm2(f0)/dsqrt(dble(n))
    r = normf/normfold

    ! Adjust eta
    if (etamax .gt. 0.d0) then
      etaold = eta
      etanew = gg*r**2
      if (gg*etaold**2 .gt. .1d0) then
        etanew = max(etanew,gg*etaold**2)
      end if
      eta = max(min(etanew,etamax),.5d0*stop_tol/normf)
    end if

  end do ! End of outer nonlinear iteration

  ! Finalise
  sol = xk
  
  ferr = normf
  ierr = 0
  if (ferr .gt. stop_tol) ierr = 1

end subroutine nsolgm
!=============================================================================
subroutine fdgmres(f0, xc, step, n, eta, maxit, residual)

  ! GMRES solver - Inner linear iteration in a nonlinear solver
  ! Fabio Bacchini, October 2018

  implicit none

  ! Input variables
  integer,intent(in)                          :: n ! Dimension of the linear system
  double precision,intent(in),dimension(n)    :: xc, f0 ! Current solution and residual vectors
  double precision,intent(out),dimension(n)   :: step ! Solution of the linear system
  double precision,intent(in)                 :: eta ! Tolerance
  integer,intent(in)                          :: maxit ! Maximum iterations
  ! Local variables
  double precision,dimension(maxit,maxit)     :: A, h
  double precision,dimension(n,maxit)         :: v
  double precision,dimension(maxit+1)         :: c, s, g
  double precision,dimension(maxit)           :: y
  double precision,dimension(n)               :: x, rhs
  double precision                            :: rho, normv1, normv2, hr, nu, err_tol
  integer,dimension(maxit)                    :: P
  integer                                     :: ii, jj, d, ierr
  external residual

  ! Initialise solver
  rhs = -f0 ! RHS of the linear system
  rho = norm2(rhs) ! Initial residual norm
  err_tol = eta*rho

  ! Initialise vectors
  h = 0.d0
  v = 0.d0
  c = 0.d0
  s = 0.d0
  g = 0.d0
  g(1) = rho
  x = 0.d0

  ! Check for early convergence
  if (rho .lt. err_tol) then
    step = 0.d0
    return
  end if

  v(:,1) = rhs/rho
  ii = 0

  ! Main loop
  do while((rho .gt. err_tol) .and. (ii .lt. maxit)) ! Start of (inner) linear iteration
    ii = ii + 1

    ! Directional derivative
    call dirder(xc, v(:,ii), f0, n, v(:,ii+1), residual)
    normv1 = norm2(v(:,ii+1))

    ! Modified Gram-Schmidt orthogonalisation
    do jj=1,ii
      h(jj,ii) = dot_product(v(:,jj),v(:,ii+1))
      v(:,ii+1) = v(:,ii+1) - h(jj,ii)*v(:,jj)
    end do
    h(ii+1,ii) = norm2(v(:,ii+1))
    normv2 = h(ii+1,ii)

    ! Reorthogonalisation (if necessary)
    if (normv1 + normv2*1.d-3 .eq. normv1) then ! Probably needs to be changed to normv1 < some tol
      do jj=1,ii
        hr = dot_product(v(:,jj),v(:,ii+1))
        h(jj,ii) = h(jj,ii) + hr
        v(:,ii+1) = v(:,ii+1) - hr*v(:,jj)
      end do
      h(ii+1,ii) = norm2(v(:,ii+1))
    end if

    if (h(ii+1,ii) .ne. 0.d0) then
      v(:,ii+1) = v(:,ii+1)/h(ii+1,ii)
!    else
!      write(*,*) "Happy breakdown needed on v"
    end if

    ! Givens rotation
    if (ii .gt. 1) call givapp(c(1:ii-1), s(1:ii-1), h(1:ii,ii), ii-1)
    nu = norm2(h(ii:ii+1,ii))

    if (nu .ne. 0.d0) then
      c(ii) = h(ii,ii)/nu
      s(ii) = -h(ii+1,ii)/nu
      h(ii,ii) = c(ii)*h(ii,ii) - s(ii)*h(ii+1,ii)
      h(ii+1,ii) = 0.d0
      call givapp(c(ii), s(ii), g(ii:ii+1), 1)
!    else
!      write(*,*) "Happy breakdown needed on nu"
    end if

    ! Update error
    rho = abs(g(ii+1))
!if (rho .ne. rho) then
!write(*,*) "NaN rho in inner it"
!call mpistop()
!end if
!write(*,*) "Inner it: it, rho, err_tol", ii, rho, err_tol

  end do ! End of (inner) linear iteration

  ! Finalise: solve the linear problem
  A(1:ii,1:ii) = h(1:ii,1:ii)
!    write(*,*) "A",A(1:ii,1:ii)
!    write(*,*) "g",g(1:ii)
    
  call LUDCMP_d(A(1:ii,1:ii),ii,P(1:ii),d,ierr) ! LU decomposition of h
  ! Do the rest only if A is nonsingular
  if (ierr .ne. 1) then
    y(1:ii) = g(1:ii)
    call LUBKSB_d(A(1:ii,1:ii),ii,P(1:ii),y(1:ii)) ! Solve hy = g
  else
!    write(*,*) "Singular matrix while inverting!"
    return
  end if

  ! Update solution
  step = x + matmul(v(1:n,1:ii),y(1:ii))

end subroutine fdgmres
!=============================================================================
subroutine dirder(x, v, f0, n, vnew, residual)

  ! Directional derivative - approximates a Jacobian matrix
  ! Fabio Bacchini, October 2018

  implicit none

  ! Input variables
  integer,intent(in)                        :: n ! Number of variables
  double precision,intent(in),dimension(n)  :: v, f0, x
  double precision,intent(out),dimension(n) :: vnew
  ! Local variables
  double precision,dimension(n)             :: delta, f1
  double precision                          :: eps
  external residual

  ! Initial sanity check
  if (norm2(v) .eq. 0.d0) then
    vnew = 0.d0
    return
  end if
  
  eps = 1.d-7!2.d0*sqrt(epsilon(1.d0))
  ! Scale increment -- SEEMS TO CAUSE TROUBLE. CAN BE AVOIDED ASSUMING norm2(x)=norm2(v)
  eps = eps/norm2(v)
  if (norm2(x) .gt. 0.d0) eps = eps*norm2(x)
!write(*,*) norm2(x),norm2(v)

  ! Compute derivatives
  delta = x + eps*v
  call residual(delta, f1, n)
  vnew = (f1 - f0) / eps

end subroutine dirder
!=============================================================================
subroutine givapp(c, s, v, n)

  ! Performs n Givens rotations
  ! Fabio Bacchini, October 2018

  implicit none

  ! Input variables
  integer,intent(in)                            :: n ! Number of variables
  double precision,intent(inout),dimension(n+1) :: v
  double precision,intent(in),dimension(n)      :: c, s
  ! Local variables
  double precision                              :: w1, w2
  integer                                       :: ii

  do ii=1,n
    w1 = c(ii)*v(ii) - s(ii)*v(ii+1)
    w2 = s(ii)*v(ii) + c(ii)*v(ii+1)
    v(ii:ii+1) = (/ w1, w2 /)
  end do

end subroutine givapp
!=============================================================================
!subroutine residual(x, res, n)
!
!!  implicit none
!
!  ! Input variables
!  integer,intent(in)                        :: n ! Number of variables
!  double precision,intent(in),dimension(n)  :: x
!  double precision,intent(out),dimension(n) :: res
!
!  res = x
!
!end subroutine residual
!=============================================================================
end module mod_ngmres
!=============================================================================
