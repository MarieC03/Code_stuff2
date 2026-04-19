!weno_eps_machine = 1.0e-52  to be same as FIL
module mod_weno
  implicit none
  private

  public :: WENO5limiter
  public :: WENO5limitervar
  public :: WENOZPlimitervar
  public :: WENOZ5limitervar
contains

  subroutine WENO5limiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,var,rec)
    include 'amrvacdef.f'
  
    integer, intent(in)             :: ixI^L, iL^L, idims, var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,3), d_array(3)
    double precision                :: beta(ixI^S,3), beta_coeff(2)
    double precision                :: tau(ixI^S), tmp(ixI^S)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,3), alpha_sum(ixI^S), flux(ixI^S)
    integer                         :: i
    double precision, parameter     :: weno_eps_machine = 1.0d-52
    !double precision, parameter     :: weno_eps_machine = 1.0d-18
    double precision                :: lambda
    double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    lambda = dxdim**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    
    if (rec) then
      !   reconstruction variation
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      !   interpolation variation
      d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
      u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
      u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
      u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    endif
    
    !> left side
    f_array(iL^S,1) = u1_coeff(1) * w(iLmm^S) + u1_coeff(2) * w(iLm^S) + u1_coeff(3) * w(iL^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLm^S)  + u2_coeff(2) * w(iL^S)  + u2_coeff(3) * w(iLp^S)
    f_array(iL^S,3) = u3_coeff(1) * w(iL^S)   + u3_coeff(2) * w(iLp^S) + u3_coeff(3) * w(iLpp^S)  
  
    beta(iL^S,1) = beta_coeff(1) * (w(iLmm^S) + w(iL^S) - 2.0d0*w(iLm^S))**2 &
         + beta_coeff(2) * (w(iLmm^S) - 4.0d0 * w(iLm^S) + 3.0d0*w(iL^S))**2
    beta(iL^S,2) = beta_coeff(1) * (w(iLm^S) + w(iLp^S) - 2.0d0 * w(iL^S))**2 &
         + beta_coeff(2) * (w(iLm^S) - w(iLp^S))**2
    beta(iL^S,3) = beta_coeff(1) * (w(iL^S) + w(iLpp^S) - 2.0d0 * w(iLp^S))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S) - 4.0d0 * w(iLp^S) + w(iLpp^S))**2
 
    alpha_sum(iL^S) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
         alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(2)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                      (beta(iL^S,i) + weno_eps_machine))**2)
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(3)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
      do i = 1,3
        tmp(iL^S) = (tau(iL^S) + weno_eps_machine) / (beta(iL^S,i) + weno_eps_machine)
        alpha_array(iL^S,i) = d_array(i) * (1.0d0 + tmp(iL^S)**2 + lambda/tmp(iL^S))
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    end select

    flux(iL^S) = 0.0d0
    do i = 1,3
      flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
  
    !> left value at right interface
    wLC(iL^S) = flux(iL^S)
  
    !> right side
    f_array(iL^S,1) = u1_coeff(1) * w(iLppp^S) + u1_coeff(2) * w(iLpp^S) + u1_coeff(3) * w(iLp^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLpp^S)  + u2_coeff(2) * w(iLp^S)  + u2_coeff(3) * w(iL^S)
    f_array(iL^S,3) = u3_coeff(1) * w(iLp^S)   + u3_coeff(2) * w(iL^S)   + u3_coeff(3) * w(iLm^S)  
  
    beta(iL^S,1) = beta_coeff(1) * (w(iLppp^S) + w(iLp^S) - 2.0d0*w(iLpp^S))**2 &
         + beta_coeff(2) * (w(iLppp^S) - 4.0d0 * w(iLpp^S) + 3.0d0*w(iLp^S))**2
    beta(iL^S,2) = beta_coeff(1) * (w(iLpp^S) + w(iL^S) - 2.0d0 * w(iLp^S))**2 &
         + beta_coeff(2) * (w(iLpp^S) - w(iL^S))**2
    beta(iL^S,3) = beta_coeff(1) * (w(iLp^S) + w(iLm^S) - 2.0d0 * w(iL^S))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S) - 4.0d0 * w(iL^S) + w(iLm^S))**2
  
    alpha_sum(iL^S) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(2) 
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                      (beta(iL^S,i) + weno_eps_machine))**2)
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(3)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
      do i = 1,3
        tmp(iL^S) = (tau(iL^S) + weno_eps_machine) / (beta(iL^S,i) + weno_eps_machine)
        alpha_array(iL^S,i) = d_array(i) * (1.0d0 + tmp(iL^S)**2 + lambda/tmp(iL^S))
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    end select
    flux(iL^S) = 0.0d0
    do i = 1,3
      flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
  
    !> right value at right interface
    wRC(iL^S) = flux(iL^S)

  end subroutine WENO5limiter

  subroutine minmod(ixI^L,ixO^L,a,b,minm)

    include 'amrvacdef.f'

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S)
    double precision, intent(out):: minm(ixI^S)

    minm(ixO^S) = (sign(one,a(ixO^S))+sign(one,b(ixO^S)))/2.0d0 * &
         min(abs(a(ixO^S)),abs(b(ixO^S)))

  end subroutine minmod

  subroutine median(ixI^L,ixO^L,a,b,c,med)

    include 'amrvacdef.f'

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S), c(ixI^S)
    double precision, intent(out):: med(ixI^S)

    double precision             :: tmp1(ixI^S),tmp2(ixI^S)

    tmp1(ixO^S) = b(ixO^S) - a(ixO^S); tmp2(ixO^S) = c(ixO^S) - a(ixO^S)

    med(ixO^S) = a(ixO^S) + (sign(one,tmp1(ixO^S))+sign(one,tmp2(ixO^S)))/2.0d0 * &
         min(abs(tmp1(ixO^S)),abs(tmp2(ixO^S)))

  end subroutine median


subroutine WENO5limitervar(ixI^L,iL^L,idims,w,wLC,wRC)

  include 'amrvacdef.f'

  integer, intent(in)             :: ixI^L, iL^L, idims
  double precision, intent(in)    :: w(ixI^S)

  double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S)
  ! .. local ..
  integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
  double precision                :: f_array(ixI^S,3), d_array(3), beta(ixI^S,3), beta_coeff(2), alpha_array(ixI^S,3), tau_5(ixI^S), alpha_sum(ixI^S), usum(ixI^S), flux(ixI^S)
  integer                         :: i
  double precision, parameter     :: weno_eps_machine = 1.0d-52
  !----------------------------------------------------------------------------



  ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.

  iLm^L=iL^L-kr(idims,^D);
  iLmm^L=iLm^L-kr(idims,^D);
  iLp^L=iL^L+kr(idims,^D);
  iLpp^L=iLp^L+kr(idims,^D);



  f_array(iL^S,1) = 3.0d0/8.0d0 * w(iLmm^S) - 10.0d0/8.0d0 * w(iLm^S) + 15.0d0/8.0d0 * w(iL^S) !first stencil
  f_array(iL^S,2) = -1.0d0/8.0d0 * w(iLm^S) + 6.0d0/8.0d0 * w(iL^S) + 3.0d0/8.0d0 * w(iLp^S)
  f_array(iL^S,3) = 3.0d0/8.0d0 * w(iL^S) + 6.0d0/8.0d0 * w(iLp^S) - 1.0d0/8.0d0 * w(iLpp^S)

  d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)


  beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)

  beta(iL^S,1) = beta_coeff(1) * (w(iLmm^S) + w(iL^S) - 2.0d0*w(iLm^S))**2 &
       + beta_coeff(2) * (w(iLmm^S) - 4.0d0 * w(iLm^S) + 3.0d0*w(iL^S))**2

  beta(iL^S,2) = beta_coeff(1) * (w(iLm^S) + w(iLp^S) - 2.0d0 * w(iL^S))**2 &
       + beta_coeff(2) * (w(iLm^S) - w(iLp^S))**2

  beta(iL^S,3) = beta_coeff(1) * (w(iL^S) + w(iLpp^S) - 2.0d0 * w(iLp^S))**2 &
       + beta_coeff(2) * (3.0d0 * w(iL^S) - 4.0d0 * w(iLp^S) + w(iLpp^S))**2


  alpha_sum(iL^S) = 0.0d0
  do i = 1,3
     alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
     alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
  end do


  flux(iL^S) = 0.0d0
  do i = 1,3
     flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
  end do

  !left value at right interface
  wLC(iL^S) = flux(iL^S)

  !now other side

  iLppp^L=iLpp^L+kr(idims,^D);

  f_array(iL^S,1) = 3.0d0/8.0d0 * w(iLppp^S) - 10.0d0/8.0d0 * w(iLpp^S) + 15.0d0/8.0d0 * w(iLp^S) !first stencil
  f_array(iL^S,2) = -1.0d0/8.0d0 * w(iLpp^S) + 6.0d0/8.0d0 * w(iLp^S) + 3.0d0/8.0d0 * w(iL^S)
  f_array(iL^S,3) = 3.0d0/8.0d0 * w(iLp^S) + 6.0d0/8.0d0 * w(iL^S) - 1.0d0/8.0d0 * w(iLm^S)

  beta(iL^S,1) = beta_coeff(1) * (w(iLppp^S) + w(iLp^S) - 2.0d0*w(iLpp^S))**2 &
       + beta_coeff(2) * (w(iLppp^S) - 4.0d0 * w(iLpp^S) + 3.0d0*w(iLp^S))**2

  beta(iL^S,2) = beta_coeff(1) * (w(iLpp^S) + w(iL^S) - 2.0d0 * w(iLp^S))**2 &
       + beta_coeff(2) * (w(iLpp^S) - w(iL^S))**2

  beta(iL^S,3) = beta_coeff(1) * (w(iLp^S) + w(iLm^S) - 2.0d0 * w(iL^S))**2 &
       + beta_coeff(2) * (3.0d0 * w(iLp^S) - 4.0d0 * w(iL^S) + w(iLm^S))**2

  alpha_sum(iL^S) = 0.0d0
  do i = 1,3
     alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
     alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
  end do

  flux(iL^S) = 0.0d0
  do i = 1,3
     flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
  end do

  !right value at right interface
  wRC(iL^S) = flux(iL^S)

end subroutine WENO5limitervar


subroutine WENOZPlimitervar(ixI^L,iL^L,idims,dxdim,w,wLC,wRC)
  ! WENOZ+ Limiter, see Acker, Borges and Costa (2016)
  ! Needs at least three ghost cells.  Set dixB=3.

  include 'amrvacdef.f'

  integer, intent(in)             :: ixI^L, iL^L, idims
  double precision, intent(in)    :: dxdim
  double precision, intent(in)    :: w(ixI^S)
  double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S)
  ! .. local ..
  integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
  double precision                :: f_array(ixI^S,3), d_array(3), beta(ixI^S,3)
  double precision                :: beta_coeff(2), alpha_array(ixI^S,3), tau(ixI^S)
  double precision                :: alpha_sum(ixI^S), tmp(ixI^S), flux(ixI^S)
  double precision                :: lambda
  integer                         :: i, iw
  double precision, parameter     :: weno_eps_machine = 1.0d-52
  !double precision, parameter     :: weno_eps_machine = 1.0d-42
  double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0
  !----------------------------------------------------------------------------

  ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.

  iLm^L=iL^L-kr(idims,^D);
  iLmm^L=iLm^L-kr(idims,^D);
  iLp^L=iL^L+kr(idims,^D);
  iLpp^L=iLp^L+kr(idims,^D);
  iLppp^L=iLpp^L+kr(idims,^D);

  d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
  beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)

  lambda = dxdim**weno_dx_exp

  ! ==================================================
  ! Left side
  ! ==================================================

  f_array(iL^S,1) = 3.0d0/8.0d0 * w(iLmm^S) - 10.0d0/8.0d0 * w(iLm^S) + 15.0d0/8.0d0 * w(iL^S) !first stencil
  f_array(iL^S,2) = -1.0d0/8.0d0 * w(iLm^S) + 6.0d0/8.0d0 * w(iL^S) + 3.0d0/8.0d0 * w(iLp^S)
  f_array(iL^S,3) = 3.0d0/8.0d0 * w(iL^S) + 6.0d0/8.0d0 * w(iLp^S) - 1.0d0/8.0d0 * w(iLpp^S)


  beta(iL^S,1) = beta_coeff(1) * (w(iLmm^S) + w(iL^S) - 2.0d0*w(iLm^S))**2 &
       + beta_coeff(2) * (w(iLmm^S) - 4.0d0 * w(iLm^S) + 3.0d0*w(iL^S))**2
  beta(iL^S,2) = beta_coeff(1) * (w(iLm^S) + w(iLp^S) - 2.0d0 * w(iL^S))**2 &
       + beta_coeff(2) * (w(iLm^S) - w(iLp^S))**2
  beta(iL^S,3) = beta_coeff(1) * (w(iL^S) + w(iLpp^S) - 2.0d0 * w(iLp^S))**2 &
       + beta_coeff(2) * (3.0d0 * w(iL^S) - 4.0d0 * w(iLp^S) + w(iLpp^S))**2

  tau(iL^S) = abs( beta(iL^S,1) - beta(iL^S,3) )

  alpha_sum(iL^S) = 0.0d0
  do i = 1,3
     tmp(iL^S) = (tau(iL^S) + weno_eps_machine ) / ( beta(iL^S,i) + weno_eps_machine )
     alpha_array(iL^S,i) = d_array(i) * (1.0d0 + tmp(iL^S)**2 + lambda/tmp(iL^S) )
     alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
  end do


  flux(iL^S) = 0.0d0
  do i = 1,3
     flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
  end do

  ! left value at right interface
  wLC(iL^S) = flux(iL^S)

  ! ==================================================
  ! Right side
  ! ==================================================

  f_array(iL^S,1) = 3.0d0/8.0d0 * w(iLppp^S) - 10.0d0/8.0d0 * w(iLpp^S) + 15.0d0/8.0d0 * w(iLp^S) !first stencil
  f_array(iL^S,2) = -1.0d0/8.0d0 * w(iLpp^S) + 6.0d0/8.0d0 * w(iLp^S) + 3.0d0/8.0d0 * w(iL^S)
  f_array(iL^S,3) = 3.0d0/8.0d0 * w(iLp^S) + 6.0d0/8.0d0 * w(iL^S) - 1.0d0/8.0d0 * w(iLm^S)

  beta(iL^S,1) = beta_coeff(1) * (w(iLppp^S) + w(iLp^S) - 2.0d0*w(iLpp^S))**2 &
       + beta_coeff(2) * (w(iLppp^S) - 4.0d0 * w(iLpp^S) + 3.0d0*w(iLp^S))**2
  beta(iL^S,2) = beta_coeff(1) * (w(iLpp^S) + w(iL^S) - 2.0d0 * w(iLp^S))**2 &
       + beta_coeff(2) * (w(iLpp^S) - w(iL^S))**2
  beta(iL^S,3) = beta_coeff(1) * (w(iLp^S) + w(iLm^S) - 2.0d0 * w(iL^S))**2 &
       + beta_coeff(2) * (3.0d0 * w(iLp^S) - 4.0d0 * w(iL^S) + w(iLm^S))**2


  tau(iL^S) = abs( beta(iL^S,1) - beta(iL^S,3) )

  alpha_sum(iL^S) = 0.0d0
  do i = 1,3
     tmp(iL^S) = (tau(iL^S) + weno_eps_machine ) / ( beta(iL^S,i) + weno_eps_machine )
     alpha_array(iL^S,i) = d_array(i) * (1.0d0 + tmp(iL^S)**2 + lambda/tmp(iL^S) )
     alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
  end do

  flux(iL^S) = 0.0d0
  do i = 1,3
     flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
  end do

  ! right value at right interface
  wRC(iL^S) = flux(iL^S)

end subroutine WENOZPlimitervar

subroutine WENOZ5limitervar(ixI^L,iL^L,idims,dxdim,w,wLC,wRC)
  ! WENOZ Limiter, see Acker, Borges+ (2018)
  ! Needs at least three ghost cells.  Set dixB=3.

  include 'amrvacdef.f'

  integer, intent(in)             :: ixI^L, iL^L, idims
  double precision, intent(in)    :: dxdim
  double precision, intent(in)    :: w(ixI^S)
  double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S)
  ! .. local ..
  integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
  double precision                :: f_array(ixI^S,3), d_array(3), beta(ixI^S,3)
  double precision                :: beta_coeff(2), alpha_array(ixI^S,3), tau(ixI^S)
  double precision                :: alpha_sum(ixI^S), tmp(ixI^S), flux(ixI^S)
  double precision                :: lambda
  integer                         :: i, iw
  double precision, parameter     :: weno_eps_machine = 1.0d-52
  !double precision, parameter     :: weno_eps_machine = 1.0d-42
  double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0
  !----------------------------------------------------------------------------

  ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.

  iLm^L=iL^L-kr(idims,^D);
  iLmm^L=iLm^L-kr(idims,^D);
  iLp^L=iL^L+kr(idims,^D);
  iLpp^L=iLp^L+kr(idims,^D);
  iLppp^L=iLpp^L+kr(idims,^D);

  d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
  beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)

  lambda = dxdim**weno_dx_exp

  ! ==================================================
  ! Left side
  ! ==================================================

  f_array(iL^S,1) = 3.0d0/8.0d0 * w(iLmm^S) - 10.0d0/8.0d0 * w(iLm^S) + 15.0d0/8.0d0 * w(iL^S) !first stencil
  f_array(iL^S,2) = -1.0d0/8.0d0 * w(iLm^S) + 6.0d0/8.0d0 * w(iL^S) + 3.0d0/8.0d0 * w(iLp^S)
  f_array(iL^S,3) = 3.0d0/8.0d0 * w(iL^S) + 6.0d0/8.0d0 * w(iLp^S) - 1.0d0/8.0d0 * w(iLpp^S)


  beta(iL^S,1) = beta_coeff(1) * (w(iLmm^S) + w(iL^S) - 2.0d0*w(iLm^S))**2 &
       + beta_coeff(2) * (w(iLmm^S) - 4.0d0 * w(iLm^S) + 3.0d0*w(iL^S))**2
  beta(iL^S,2) = beta_coeff(1) * (w(iLm^S) + w(iLp^S) - 2.0d0 * w(iL^S))**2 &
       + beta_coeff(2) * (w(iLm^S) - w(iLp^S))**2
  beta(iL^S,3) = beta_coeff(1) * (w(iL^S) + w(iLpp^S) - 2.0d0 * w(iLp^S))**2 &
       + beta_coeff(2) * (3.0d0 * w(iL^S) - 4.0d0 * w(iLp^S) + w(iLpp^S))**2

  tau(iL^S) = abs( beta(iL^S,1) - beta(iL^S,3) )

  alpha_sum(iL^S) = 0.0d0
  do i = 1,3
    alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                  (beta(iL^S,i) + weno_eps_machine))**2)
    alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
  end do

  flux(iL^S) = 0.0d0
  do i = 1,3
     flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
  end do

  ! left value at right interface
  wLC(iL^S) = flux(iL^S)

  ! ==================================================
  ! Right side
  ! ==================================================

  f_array(iL^S,1) = 3.0d0/8.0d0 * w(iLppp^S) - 10.0d0/8.0d0 * w(iLpp^S) + 15.0d0/8.0d0 * w(iLp^S) !first stencil
  f_array(iL^S,2) = -1.0d0/8.0d0 * w(iLpp^S) + 6.0d0/8.0d0 * w(iLp^S) + 3.0d0/8.0d0 * w(iL^S)
  f_array(iL^S,3) = 3.0d0/8.0d0 * w(iLp^S) + 6.0d0/8.0d0 * w(iL^S) - 1.0d0/8.0d0 * w(iLm^S)

  beta(iL^S,1) = beta_coeff(1) * (w(iLppp^S) + w(iLp^S) - 2.0d0*w(iLpp^S))**2 &
       + beta_coeff(2) * (w(iLppp^S) - 4.0d0 * w(iLpp^S) + 3.0d0*w(iLp^S))**2
  beta(iL^S,2) = beta_coeff(1) * (w(iLpp^S) + w(iL^S) - 2.0d0 * w(iLp^S))**2 &
       + beta_coeff(2) * (w(iLpp^S) - w(iL^S))**2
  beta(iL^S,3) = beta_coeff(1) * (w(iLp^S) + w(iLm^S) - 2.0d0 * w(iL^S))**2 &
       + beta_coeff(2) * (3.0d0 * w(iLp^S) - 4.0d0 * w(iL^S) + w(iLm^S))**2


  tau(iL^S) = abs( beta(iL^S,1) - beta(iL^S,3) )

  alpha_sum(iL^S) = 0.0d0
  do i = 1,3
    alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                  (beta(iL^S,i) + weno_eps_machine))**2)
    alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
  end do


  flux(iL^S) = 0.0d0
  do i = 1,3
     flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
  end do

  ! right value at right interface
  wRC(iL^S) = flux(iL^S)

end subroutine WENOZ5limitervar

end module mod_weno
