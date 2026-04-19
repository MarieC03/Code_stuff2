module mod_rootfinding
  implicit none
contains

  recursive subroutine rootfinding_brent(z, zmin, zmax, tolerance, iter_max, return_code, func)
    !> root
    double precision, intent(out):: z
    !> the range of the root
    double precision, intent(in) :: zmin, zmax
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximum iterations
    integer, intent(in)          :: iter_max
    !> return code:
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = the root is not bracketed;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
    end interface

    integer :: it_root
    double precision :: a, b, c, d, e
    double precision :: fa, fb, fc
    double precision :: p, q, r, s
    double precision :: xm, tol
    double precision, parameter :: eps = epsilon(zmin) ! machine precision

    return_code = -1

    a  = zmin
    b  = zmax
    fa = func(a)
    fb = func(b)

    ! check if the bounds are already the root
    if ( dabs(fa) == 0.0d0 ) then
       z = a
       return_code = 0
       return
    else if ( dabs(fb) == 0.0d0 ) then
       z = b
       return_code = 0
       return
    else if ( fa * fb > 0.0d0 ) then
       ! the root is not bracketed
       return_code = 3
       return
       !call mpistop("error in brent: the root is not bracketed!")
    end if

    c  = b
    fc = fb

    do it_root = 1, iter_max

       if ( fb * fc > 0.0d0 ) then
          ! Rename a, b, c and adjust bounding interval d
          c  = a
          fc = fa
          d  = b - a
          e  = d
       end if

       if ( dabs(fc) < dabs(fb) ) then
          a  = b
          b  = c
          c  = a
          fa = fb
          fb = fc
          fc = fa
       end if

       ! Convergence check.
       tol = 2.0d0 * eps * dabs(b) + 0.5d0 * tolerance
       xm  = 0.5d0 * ( c - b )
       if ( dabs(xm) <= tol .or. fb == 0.0d0 ) then
          z = b
          return_code = 0
          return
       end if

       if ( dabs(e) >= tol .and. dabs(fa) > dabs(fb) ) then
          ! Attempt inverse quadratic interpolation
          s = fb / fa
          if (a == c) then
             p = 2.0d0 * xm * s
             q = 1.0d0 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * (2.0d0 * xm * q * (q-r) - (b-a) * (r - 1.0d0))
             q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
          end if
          if (p > 0.0d0) q = -q  ! Check whether in bounds.
          p = dabs(p)
          if ( 2.0d0 * p < min(3.0d0 * xm * q - dabs(tol*q), dabs(e*q)) ) then
             ! Accept interpolation
             e = d
             d = p / q
          else
             ! Interpolation failed, use bisection
             d = xm
             e = d
          end if
       else
          ! Bounds decreasing too slowly, use bisection
          d = xm
          e = d
       end if
       a  = b  ! move last best guess to a
       fa = fb
       b  = b + merge(d, sign(tol,xm), dabs(d) > tol) ! evaluate new trail root
       fb = func(b)
    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_brent

  recursive subroutine rootfinding_illinois_old(z, zmin, zmax, tolerance, iter_max, return_code, func)
    !> root
    double precision, intent(out):: z
    !> the range of the root
    double precision, intent(in) :: zmin, zmax
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximun iterations
    integer, intent(in)          :: iter_max
    !> return code:
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = the root is not bracketed;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
    end interface

    integer :: it_root
    integer :: side
    double precision :: zm, zp
    double precision :: f, fm, fp
    double precision :: z_new

    return_code = -1

    zm = zmin
    zp = zmax
    fm = func(zm)
    fp = func(zp)


    ! check if the bounds are already the root
    if ( dabs(fm) <= tolerance ) then
       z = zm
       return_code = 0
       return
    else if ( dabs(fp) <= tolerance ) then
       z = zp
       return_code = 0
       return
    else if ( fm * fp > 0.0d0 ) then
       ! check if the root is bracketed
       !write(*,*) "Warning in illinois: the root is not bracketed!"
       !write(*,*) "root range = ", zm, zp
       !write(*,*) "f(root) range = ", fm, fp

       if ( dabs(fp) > dabs(fm) ) then
          z = zm
          f = fm
       else
          z = zp
          f = fp
       end if

       z_new = tolerance
       do it_root=1,10
          z_new = z_new * 1.0d1
          if ( dabs(f) <= z_new) then
!             return_code = 3
             return_code = 0
             return
          end if
       end do
       return_code = 3
       return
       !call mpistop("error in illinois: the root is not bracketed!")
    end if

    do it_root = 1, iter_max

       z = (fm * zp - fp * zm) / (fm - fp)
       f = func(z)

       if (dabs(zp-zm) <= (tolerance * dabs(zp+zm)) .or. &
          (dabs(f) <= tolerance)  ) then
          return_code = 0
          return
       end if

       if ( (f * fp) > 0.0d0 ) then
          !f and fp have same sign, copy z to zp
          zp = z
          fp = f
          if (side == 1) fm = fm * 0.5d0
          side = 1
       else if ( (f * fm) > 0.0d0 ) then
          !f and fm have same sign, copy z to zm
          zm = z
          fm = f
          if (side == 2) fp = fp * 0.5d0
          side = 2
       else !it looks like zero
          return_code = 0
          return
       end if

       if ( isnan(z) ) then
          return_code = 2
          return
       end if

    end do
    ! if the root is not found
    return_code = 1
    !call mpistop("fail to find the root")
  end subroutine rootfinding_illinois_old

!  Reference:
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
! from FIL
  recursive subroutine rootfinding_zero_brent(z, a, b, t, func)
    !> root
    double precision, intent(out):: z
    !> the range of the root [A, B]
    double precision, intent(in) :: a, b
    !> the tolerance
    double precision, intent(in) :: t 
    !> maximun iterations
    !integer, intent(in)          :: iter_max
    !> return code:
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = the root is not bracketed;
    !integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
    end interface
   
    double precision :: c, d, e, fa, fb, fc, m, p, q, r, s, sa, sb, tol
    double precision, parameter :: macheps = 2.22045d-16 ! in c++ it is 2.22045e-16

    sa = a
    sb = b
    fa = func(sa)
    fb = func(sb)

    c  = sa
    fc = fa
    e  = sb - sa
    d  = e
   

    do !infite loop

       if ( abs(fc) < abs(fb) ) then
          sa = sb
          sb = c
          c  = sa
          fa = fb
          fb = fc
          fc = fa
       endif

       tol = 2.0d0 * macheps * abs(sb) + t
       m = 0.5d0 * (c - sb)
   
       if ((abs(m) .le. tol) .or. (fb == 0.0d0)) exit

       if ((abs(e) .lt. tol) .or. (abs(fa) .le. abs(fb)) ) then 
          e = m
          d = e
       else
          s = fb/fa
          if (sa == c) then
              p = 2.0d0 * m * s
              q = 1.0d0 - s
          else
              q = fa/fc
              r = fb/fc
              p = s * (2.0d0 * m * q * (q-r) - (sb-sa) * (r-1.0d0)) 
              q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
          endif
 
          if (p > 0.0d0) then
              q = -q
          else
              p = -p
          endif

          s = e
          e = d
    
          if ( (2.0d0 * p < 3.0d0 * m * q - abs(tol * q)) .and. (p < abs(0.5 * s * q)) ) then
                 d = p / q
          else
                 e = m
                 d = e
          endif
       endif 

       sa = sb
       fa = fb

       if (tol < abs(d)) then
            sb = sb + d
       else if (m > 0.0d0) then
            sb = sb + tol
       else
            sb = sb - tol
       endif
  
       fb = func(sb) 
      
       if ( (fb > 0.0d0 .and. fc > 0.0d0 )   .or.  (fb .le. 0.0d0 .and. fc .le. 0.0d0) ) then
          c = sa
          fc = fa
          e = sb - sa
          d = e
       endif

    end do ! end while

    z = sb
    return

  end subroutine rootfinding_zero_brent


  ! note: making this recursive since this subroutine will call itself indirectly in some cases.
  recursive subroutine rootfinding_illinois(z, zmin, zmax, tolerance, iter_max, return_code, func)
    implicit none
    !> root
    double precision, intent(out):: z
    !> the range of the root
    double precision, intent(in) :: zmin, zmax
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximun iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = the root is not bracketed;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
    end interface

    integer :: it_root
    integer :: side
    double precision :: zm, zp
    double precision :: f, fm, fp
    double precision :: z_new
    
    side = 0
    return_code = -1

    zm = zmin
    zp = zmax
    fm = func(zm)
    fp = func(zp)

    ! check if the bounds are already the root
    if ( dabs(fm) == 0.0d0 ) then
       z = zm
       return_code = 0
       return
    else if ( dabs(fp) == 0.0d0 ) then
       z = zp
       return_code = 0
       return
    else if ( fm * fp > 0.0d0 ) then
       ! check if the root is bracketed
       !write(*,*) "Warning in illinois: the root is not bracketed!"
       !write(*,*) "root range = ", zm, zp
       !write(*,*) "f(root) range = ", fm, fp

!       if ( dabs(fp) > dabs(fm) ) then
!          z = zm
!          f = fm
!       else
!          z = zp
!          f = fp
!       end if
! 
!       z_new = tolerance
!       do it_root=1,10
!          z_new = z_new * 1.0d1
!          if ( dabs(f) <= z_new) then
!             return_code = 0
!             return
!          end if
!       end do
       return_code = 3
       return
       !call mpistop("error in illinois: the root is not bracketed!")
    end if

    do it_root = 1, iter_max

       z = (fm * zp - fp * zm) / (fm - fp)
       f = func(z)
       
       if ( dabs(zp-zm) <= (tolerance * 0.5d0 * dabs(zp+zm)) .or. &
            f == 0.0d0 ) then
          return_code = 0
          return
       end if

       if ( (f * fp) > 0.0d0 ) then
          !f and fp have same sign, copy z to zp
          zp = z
          fp = f
          if (side == 1) fm = fm * 0.5d0
          side = 1
       else if ( (f * fm) > 0.0d0 ) then
          !f and fm have same sign, copy z to zm
          zm = z
          fm = f
          if (side == 2) fp = fp * 0.5d0
          side = 2
       else !it looks like zero
          return_code = 0
          return
       end if

       if ( isnan(z) ) then
          return_code = 2
          return
       end if

    end do
    ! if the root is not found
    return_code = 1
    !call mpistop("fail to find the root")
  end subroutine rootfinding_illinois

  subroutine rootfinding_newton_raphson(z, tolerance, iter_max, return_code, func, dev_func)
    !> input as initial guess, output as root
    double precision, intent(inout):: z
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximun iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
       function dev_func(z)
          implicit none
          double precision :: dev_func
          double precision, intent(in) :: z
       end function dev_func
    end interface

    integer :: it_root
    double precision :: f, df

    return_code = -1

    do it_root = 1, iter_max

       if ( isnan(z) ) then
          return_code = 2
          return
       end if
       
       f  = func(z)
       df = dev_func(z)

       if ( dabs(f)  <  tolerance ) then
          ! the root is found!
          return_code = 0
          return
       end if

       z = z - f / df
    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_newton_raphson

  subroutine rootfinding_constrained_newton_raphson(z, zm, zp, tolerance, iter_max, return_code, func, dev_func)
    !> input as initial guess, output as root
    double precision, intent(inout):: z
    !> the range of the root
    double precision, intent(inout) :: zm, zp
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximum iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = the root is not bracketed;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
       function dev_func(z)
          implicit none
          double precision :: dev_func
          double precision, intent(in) :: z
       end function dev_func
    end interface

    integer :: it_root
    double precision :: f, df, z_new, stopping
    double precision :: fm, fp

    return_code = -1

    fm = func(zm)
    fp = func(zp)
    
    ! check if the bounds are already the root
    if ( dabs(fm) <= tolerance ) then
       z = zm
       return_code = 0
       return
    else if ( dabs(fp) <= tolerance ) then
       z = zp
       return_code = 0
       return
    else if ( fm * fp > 0.0d0 ) then
       ! check if the root is bracketed
       write(*,*) "root range = ", zm, zp
       write(*,*) "f(root) range = ", fm, fp
       return_code = 3
       return
       !call mpistop("error in constrained Newton Raphson: the root is not bracketed!")
    end if

    ! initial guess
    z = ( zp + zm ) * 0.5d0
    z_new = huge(0.0d0)
    stopping = huge(0.0d0)

    do it_root = 1, iter_max
       
       f  = func(z)
       df = dev_func(z)

       if ( dabs(f)  <=  tolerance .or. &
            dabs( z_new - z )  <=  dabs( z * tolerance) ) then
          return_code = 0
          return
       end if

       ! Newton-Raphson step
       z_new = z - f / df

       if ( z_new > zp .or. z_new < zm ) then
          ! the root is outside the range, use bisection here
          z_new = ( zp + zm ) * 0.5d0
          f  = func(z_new)
          if ( f * fm > 0.0d0 ) then
             zm = z_new
             fm = f
          else
             zp = z_new
             fp = f
          end if
          z = z_new !KEN
       end if

       if ( isnan(z_new) ) then
          ! if z_new is NaN, we take the previous z, and quit the iteration
          return_code = 0
          return
       else
          z = z_new
       end if

    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_constrained_newton_raphson

subroutine rootfinding_constrained_newton_raphson_M1( &
     z, zm, zp, tolerance, iter_max, return_code, func, dev_func)

  implicit none

  ! Input / output
  double precision, intent(inout) :: z
  double precision, intent(inout) :: zm, zp
  double precision, intent(in)    :: tolerance
  integer,          intent(in)    :: iter_max
  integer,          intent(out)   :: return_code

  interface
     function func(x)
       implicit none
       double precision :: func
       double precision, intent(in) :: x
     end function func

     function dev_func(x)
       implicit none
       double precision :: dev_func
       double precision, intent(in) :: x
     end function dev_func
  end interface

  ! Locals
  integer :: iter
  double precision :: xl, xh, x0, x1
  double precision :: fl, fh
  double precision :: f0, df0

  return_code = -1

  ! Initial bracket
  xl = zm
  xh = zp

  fl = func(xl)
  fh = func(xh)

  ! Check initial bracketing (Fortran safety, same as your original code)
  if ( fl * fh > 0.0d0 ) then
     return_code = 3
     return
  end if

  ! Initial guess: midpoint
  x0 = 0.5d0 * (xl + xh)
  x1 = x0

  do iter = 1, iter_max

     x0  = x1
     f0  = func(x0)
     df0 = dev_func(x0)

     ! NaN protection (Fortran-only, explicit)
     if ( isnan(f0) .or. isnan(df0) ) then
        return_code = 2
        z = x0
        return
     end if

     ! Newton step
     x1 = x0 - f0 / df0

     if ( isnan(x1) ) then
        return_code = 2
        z = x0
        return
     end if

     ! C++ behavior: if outside bracket, UPDATE BRACKET USING NEWTON STEP
     if ( x1 < xl .or. x1 > xh ) then
        if ( f0 * fh > 0.0d0 ) then
           xh = x1
           fh = f0
        else
           xl = x1
           fl = f0
        end if
     end if

     ! Convergence test (equivalent to stopif(x0,x1))
     if ( abs(x1 - x0) <= tolerance * max(1.0d0, abs(x0)) ) then
        z = x1
        return_code = 0
        return
     end if

  end do

  ! Not converged
  z = x1
  return_code = 1

end subroutine rootfinding_constrained_newton_raphson_M1



  subroutine rootfinding_constrained_newton_raphson_ken(z, zm, zp, tolerance, iter_max, return_code, func, dev_func)
  !> input as initial guess, output as root
  double precision, intent(inout):: z
  !> the range of the root
  double precision, intent(inout) :: zm, zp
  !> the tolerance
  double precision, intent(in) :: tolerance
  !> maximum iterations
  integer, intent(in)          :: iter_max
  !> return code: 
  !> -1 = nothing happened;
  !> 0 = the root is found;
  !> 1 = the root not found;
  !> 2 = the root is NaN;
  !> 3 = the root is not bracketed;
  integer, intent(out) :: return_code
  interface
     function func(z)
        implicit none
        double precision :: func
        double precision, intent(in) :: z
     end function func
     function dev_func(z)
        implicit none
        double precision :: dev_func
        double precision, intent(in) :: z
     end function dev_func
  end interface

  integer :: it_root
  double precision :: f, df, z_new
  double precision :: fm, fp

  return_code = -1

  fm = func(zm)
  fp = func(zp)
  
  ! check if the bounds are already the root
  if ( dabs(fm) <= tolerance ) then
     z = zm
     return_code = 0
     return
  else if ( dabs(fp) <= tolerance ) then
     z = zp
     return_code = 0
     return
  ! KEN this check is not needed
  !else if ( fm * fp > 0.0d0 ) then
     ! check if the root is bracketed
     !write(*,*) "root range = ", zm, zp
     !write(*,*) "f(root) range = ", fm, fp
     !return_code = 3
     !return
  end if

  ! initial guess
  z = ( zp + zm ) * 0.5d0

  do it_root = 1, iter_max
     
     f  = func(z)
     df = dev_func(z)

     if ( dabs(f) <= tolerance ) then
        return_code = 0
        return
     end if

     ! Newton-Raphson step
     z_new = z - f / df

     if ( isnan(z_new) ) then
        ! if z_new is NaN, we keep the previous z, and quit the iteration
        return_code = 2
        return
     end if

     if ( z_new > zp .or. z_new < zm ) then
        ! the root is outside the range, use bisection here
        z_new = ( zp + zm ) * 0.5d0
        f  = func(z_new)
        if ( f * fm > 0.0d0 ) then
           zm = z_new
           fm = f
        else
           zp = z_new
           fp = f
        end if
     end if

     ! check convergence based on change in z
     if ( dabs( z_new - z ) <= dabs( z * tolerance) ) then
        z = z_new
        return_code = 0
        return
     end if

     z = z_new

  end do
  
  ! if the root is not found
  return_code = 1
end subroutine rootfinding_constrained_newton_raphson_ken

  subroutine rootfinding_multid_newton_raphson(z, tolerance, iter_max, &
        return_code, fvec_in, fjac_in)
    use mod_lu
    !> input as initial guess, output as root
    double precision, dimension(:), intent(inout):: z
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximum iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    integer, intent(out) :: return_code
    interface
       function fvec_in(z_vec)
          implicit none
          double precision, dimension(:), intent(in)  :: z_vec
          double precision, dimension(1:size(z_vec))    :: fvec_in
       end function fvec_in
       function fjac_in(z_vec)
          implicit none
          double precision, dimension(:), intent(in)  :: z_vec
          double precision, dimension(1:size(z_vec),1:size(z_vec)) :: fjac_in
       end function fjac_in
    end interface

    integer                      :: it_root, n_dim
    integer                      :: err_flag, d_flag
    integer, dimension(size(z))  :: indx
    double precision, dimension(size(z), size(z)) :: jac
    double precision, dimension(size(z))          :: fvec, dz

    return_code = -1
    ! set dz as inf at the begining
    dz(:) = huge(0.0d0)

    n_dim = size(z)

    ! page 280 in F90
    do it_root = 1, iter_max
       
       fvec  = fvec_in( z )

       if ( maxval(dabs(fvec)) <=  tolerance .or. &
            maxval(dabs(dz)) <= tolerance ) then
          ! the root is found!
          return_code = 0
          return
       end if

       !jac= jacobian(z)
       dz(:) = -fvec(:)
       call ludcmp(jac, n_dim, indx, d_flag, err_flag)
       if (err_flag==1) then
          write(*,*) 'Rootfinding in Multi-D Newton-Raphson:', it_root
          write(*,*) 'roots:', z
          write(*,*) 'f(z):', fvec
          write(*,*) 'dz', dz
          write(*,*) 'jac(z):', jac
          write(*,*) '   '
          write(*,*) "Error in Multi-D Newton-Raphson: the jacobian is singular :/"
          return_code = 1 !fixme
          return
          !call mpistop("Error in Multi-D Newton-Raphson: the jacobian is singular :/")
       end if
       call lubksb(jac, n_dim, indx, dz)
       z = z + dz
    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_multid_newton_raphson

  subroutine rootfinding_global_multid_newton_raphson(n_dim, z, tolerance, iter_max, &
        return_code, fvec_in, fjac_in, fmin_in)
    use mod_lu
    use mod_lnsrch
    integer, intent(in)          :: n_dim
    !> input as initial guess, output as root
    double precision, dimension(1:n_dim), intent(inout):: z
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximum iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = reached local minimum;
    !> 4 = jacobian is singular;
    !> 5 = slope is larger or equals to zero;
    integer, intent(out) :: return_code
    interface
       function fvec_in(z_vec)
          implicit none
          double precision, dimension(:), intent(in)  :: z_vec
          double precision, dimension(1:size(z_vec))    :: fvec_in
       end function fvec_in
       function fjac_in(z_vec)
          implicit none
          double precision, dimension(:), intent(in)  :: z_vec
          double precision, dimension(1:size(z_vec),1:size(z_vec)) :: fjac_in
       end function fjac_in
       function fmin_in(z_vec)
          implicit none
          double precision :: fmin_in
          double precision, dimension(:), intent(in)  :: z_vec
       end function fmin_in
    end interface

    double precision, parameter                   :: tolz=epsilon(z), STPMX = 1.0d3
    integer                                       :: it_root
    integer                                       :: err_flag, d_flag
    integer, dimension(1:n_dim)                   :: indx
    double precision, dimension(1:n_dim)          :: fvec, dz, zold, g
    double precision, dimension(1:n_dim, 1:n_dim) :: jac
    double precision                              :: f, stpmax, fold
    logical                                       :: perform_lnsrch = .true.

    return_code = -1

    ! set dz as inf at the begining
    ! Initialize dz and zold to avoid any uninitialized values S
    dz(:) = huge(0.0d0)
    zold(:) = huge(0.0d0)

    !stpmax = huge(0.0d0) ! currently we dont limit the step length
    ! Set stepmax to limit step length:
    stpmax = STPMX * max( dot_product(z(:),z(:)), dble(n_dim))

    ! Initial function evaluation 
    f = fmin_in(z) ! this was here before

    do it_root = 1, iter_max
       ! Evaluate function vector
       fvec  = fvec_in( z )

       ! Check convergence criteria:
       !maxval( dabs(z(:)-zold(:))/max(dabs(z(:)), 1.0d0) ) <= tolz ) then
       if ( maxval(dabs(fvec)) <=  tolerance .or. &
            maxval(dabs(dz)) <= tolerance ) then
          ! the root is found!
          return_code = 0
          return
       end if

       ! Evaluate Jacobian matrix 
       ! Analytical jacobian: (does not work for kappa<0.5 or Q_ems<0.5)
        jac(:,:) = fjac_in(z(:))
       ! Numerical Jacobian:
       !call numerical_jacobian(fvec_in, z, jac)

       ! Compute gradient (check correctness)
       g = matmul(transpose(jac), fvec)  ! Gradient should be transpose(J) * fvec 
       !g(:) = matmul(fvec(:),jac(:,:))
       zold(:) = z(:)
       fold = f
       dz(:) = -fvec(:)  ! Negative of function values for direction

       ! LU decomposition of Jacobian
       call ludcmp(jac, n_dim, indx, d_flag, err_flag)
       if (err_flag==1) then
          return_code = 4
          return
          !call mpistop("Error in Multi-D Newton-Raphson: the jacobian is singular :/")
       end if
       call lubksb(jac, n_dim, indx, dz) ! Solve linear system for dz
       if(perform_lnsrch) then
           call lnsrch(zold, fold, g, dz, z, f, stpmax, err_flag, fmin_in) ! Perform line search 
           ! lnsrch returns new z
           if ( any(isnan(z(:))) ) then
              return_code = 2
              return
           end if
    
           ! Error handling for line search
           if (err_flag==1) then
              !write(*,*) "Warning in Global Multi-D Newton-Raphson: reached local minimum :/"
              return_code = 3
              return
           else if (err_flag==2) then
!              fold=dot_product(g(:),dz(:))
!              write(*,*) 'g dot p >= 0'
!              write(*,*) 'it', it_root
!              write(*,*) 'g dot p=', fold
!              write(*,*) 'fold=', f
!              write(*,*) 'g=', g
!              write(*,*) 'dz=',dz
!              write(*,*) 'z=',z
!              write(*,*) 'fvec=',fvec
!              write(*,*) 'jfac=',jac
              !call mpistop("Error in Multi-D Newton-Raphson: slope is larger then or equals to zero"
              return_code = 5 !2
              return
           end if
       else
            z = z + dz
       end if 

    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_global_multid_newton_raphson

  subroutine numerical_jacobian(fvec_in, z, jac)
   use mod_lu
   implicit none
   interface
       function fvec_in(z_vec)
           implicit none
           double precision, dimension(:), intent(in)  :: z_vec
           double precision, dimension(1:size(z_vec)) :: fvec_in
       end function fvec_in
   end interface
   double precision, dimension(:), intent(in)        :: z
   double precision, dimension(1:size(z),1:size(z)), intent(out) :: jac
   double precision, dimension(size(z))              :: fvec1, fvec2, z_temp
   double precision                                  :: h
   integer                                           :: i, j, n

   n = size(z)
   h = sqrt(epsilon(1.0d0))

   do i = 1, n
       z_temp = z
       z_temp(i) = z_temp(i) + h
       fvec1 = fvec_in(z_temp)

       z_temp(i) = z_temp(i) - 2*h
       fvec2 = fvec_in(z_temp)

       jac(:,i) = (fvec1 - fvec2) / (2*h)
   end do
end subroutine numerical_jacobian



  {#IFDEF ROOTFIND_POPULATION
  ! GA rootfinding, not yet implemented
  subroutine genetic_algorithm_rootfinding(n_dim, z, tolerance, iter_max, return_code, fvec_in, fjac_in, fmin_in)
   use mod_lu
   use mod_lnsrch
   integer, intent(in)          :: n_dim
   !> input as initial guess, output as root
   double precision, dimension(1:n_dim), intent(inout):: z
   !> the tolerance
   double precision, intent(in) :: tolerance
   !> maximum iterations
   integer, intent(in)          :: iter_max
   integer, intent(out) :: return_code
   interface
      function fvec_in(z_vec)
         implicit none
         double precision, dimension(:), intent(in)  :: z_vec
         double precision, dimension(1:size(z_vec))    :: fvec_in
      end function fvec_in
      function fjac_in(z_vec)
         implicit none
         double precision, dimension(:), intent(in)  :: z_vec
         double precision, dimension(1:size(z_vec),1:size(z_vec)) :: fjac_in
      end function fjac_in
      function fmin_in(z_vec)
         implicit none
         double precision :: fmin_in
         double precision, dimension(:), intent(in)  :: z_vec
      end function fmin_in
   end interface
   ! Initialize population
   double precision, dimension(pop_size, n_dim) :: population
   double precision, dimension(pop_size) :: fitness
   
   call initialize_population(population)
   
   ! Main GA loop
   do while (not converged)
       ! Evaluate fitness
       do i = 1, pop_size
           z = population(i, :)
           fitness(i) = fmin_in(z)
       end do
       
       ! Selection, crossover, and mutation
       call perform_selection(population, fitness)
       call perform_crossover(population)
       call perform_mutation(population)
       
       ! Apply Newton-Raphson to refine solutions
       do i = 1, pop_size
           z = population(i, :)
           call rootfinding_global_multid_newton_raphson(n_dim, z, tolerance, iter_max, return_code, fvec_in, fjac_in, fmin_in)
           if (return_code == 0) exit  ! Stop if root found
       end do
   end do
end subroutine genetic_algorithm_rootfinding
}


subroutine generate_candidate(z_current, z_new, temperature)
   implicit none
   double precision, dimension(:), intent(in) :: z_current
   double precision, dimension(:), intent(out) :: z_new
   double precision, intent(in) :: temperature
   integer :: n_dim, i
   double precision :: perturbation
   double precision, parameter :: perturbation_scale = 0.1d0  ! Scale factor for perturbations
   double precision :: rand_num

   n_dim = size(z_current)
  
   ! Copy current solution to new candidate
   z_new = z_current
  
   ! Generate a perturbation vector
   do i = 1, n_dim
       ! Generate a random number between 0 and 1
       call random_number(rand_num)
       ! Perturbation is scaled by temperature
       perturbation = perturbation_scale * temperature * (2.0d0 * rand_num - 1.0d0)
       z_new(i) = z_current(i) + perturbation
   end do

end subroutine generate_candidate

subroutine simulated_annealing_rootfinding(n_dim, z, tolerance, iter_max, return_code, fvec_in, fjac_in, fmin_in)
   use mod_lu
   use mod_lnsrch
   integer, intent(in)          :: n_dim
   !> input as initial guess, output as root
   double precision, dimension(1:n_dim), intent(inout):: z
   !> the tolerance
   double precision, intent(in) :: tolerance
   !> maximum iterations
   integer, intent(in)          :: iter_max
   integer, intent(out) :: return_code
   interface
      function fvec_in(z_vec)
         implicit none
         double precision, dimension(:), intent(in)  :: z_vec
         double precision, dimension(1:size(z_vec))    :: fvec_in
      end function fvec_in
      function fjac_in(z_vec)
         implicit none
         double precision, dimension(:), intent(in)  :: z_vec
         double precision, dimension(1:size(z_vec),1:size(z_vec)) :: fjac_in
      end function fjac_in
      function fmin_in(z_vec)
         implicit none
         double precision :: fmin_in
         double precision, dimension(:), intent(in)  :: z_vec
      end function fmin_in
   end interface
   ! Initialize temperature and current solution
   double precision :: initial_temp,temp, temp_min, alpha, rand_num
   double precision, dimension(1:size(z))  :: z_new,f_old,f_new
   integer :: i

   initial_temp = 1.0d0

   temp = initial_temp
   temp_min = 1.0d-8
   alpha = 0.9d0  ! Cooling rate
   
   ! Main SA loop
   do while (temp > temp_min)
       ! Generate a new candidate solution
       call generate_candidate(z, z_new, temp)
       f_old = fmin_in(z)
       f_new = fmin_in(z_new)
       call random_number(rand_num)
       ! Decide whether to accept the new solution
       do i =1, n_dim
         if (f_new(i) < f_old(i) .or. exp((f_old(i) - f_new(i)) / temp) > rand_num) then
            if(z_new(i)>0) then
              z(i) = z_new(i)
            end if
         end if
       end do
       
       ! Cool down the temperature
       temp = alpha * temp
   end do
   
   ! Apply Newton-Raphson to refine the solution
   call rootfinding_global_multid_newton_raphson(n_dim, z, tolerance, iter_max, return_code, fvec_in, fjac_in, fmin_in)

end subroutine simulated_annealing_rootfinding


subroutine multistart_rootfinding(n_dim, z, tolerance, iter_max, return_code, fvec_in, fjac_in, fmin_in, generate_random_start)
   use mod_lu
   use mod_lnsrch
   integer, intent(in)          :: n_dim
   !> input as initial guess, output as root
   double precision, dimension(1:n_dim), intent(inout):: z
   !> the tolerance
   double precision, intent(in) :: tolerance
   !> maximum iterations
   integer, intent(in)          :: iter_max
   integer, intent(out) :: return_code
   interface
      function fvec_in(z_vec)
         implicit none
         double precision, dimension(:), intent(in)  :: z_vec
         double precision, dimension(1:size(z_vec))    :: fvec_in
      end function fvec_in
      function fjac_in(z_vec)
         implicit none
         double precision, dimension(:), intent(in)  :: z_vec
         double precision, dimension(1:size(z_vec),1:size(z_vec)) :: fjac_in
      end function fjac_in
      function fmin_in(z_vec)
         implicit none
         double precision :: fmin_in
         double precision, dimension(:), intent(in)  :: z_vec
      end function fmin_in
      function generate_random_start(z_vec)
         implicit none
         double precision, dimension(:), intent(inout)  :: z_vec
         double precision, dimension(1:size(z_vec))    :: generate_random_start
      end function generate_random_start
   end interface
   ! Number of starting points
   integer, parameter :: num_starts = 10
   double precision, dimension(1:n_dim)  :: z_start
   double precision, dimension(1:n_dim)  :: rand_value
   integer :: i
   
   do i = 1, num_starts
       rand_value = generate_random_start(z_start)
       call rootfinding_global_multid_newton_raphson(n_dim, z_start, tolerance, iter_max, return_code, fvec_in, fjac_in, fmin_in)
       if (return_code == 0) then
           z = z_start
           return
       end if
   end do
   
   return_code = 1  ! If no root found after all starts
end subroutine multistart_rootfinding


end module mod_rootfinding
