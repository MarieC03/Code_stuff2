!> Module for handling problematic values in simulations, such as negative
!> pressures
module mod_small_values

  implicit none
!  private
  public

  !> Average over this many cells in each direction
  integer, public :: small_values_daverage = 1

  public :: prim_NaN_checker
  public :: small_values_error
!  public :: small_values_average

contains

  subroutine Simple_check_data_correctness_pt(w, x, cons_tmp, Location)
   use mod_eos
   include 'amrvacdef.f'
   double precision, intent(in) :: w(1:nw), cons_tmp(1:nw)
   double precision, intent(in) :: x(1:ndim)
   character(len=*), intent(in) :: Location
   logical :: stopflag 
   integer :: i, j

   stopflag = .false.

   {#IFDEF DY_SP
         ! Not diagonal ghostzones
         !if ( (x(1) .ge. 0.0d0 .and. x(2) < 0.0d0) .or. (x(1) < 0.0d0 .and. x(2) .ge. 0.0d0)) then
          if (w(psi_metric_) .lt. 1.0d0 .or. w(alp_metric_) .gt. 1.0d0) then
             write(*,*)  'w(d_), w(rho_)'
             write(*,*)  w(d_), w(rho_)
             write(*,*) 'w(u1_), w(u2_), w(u3_), w(T_eps_)'
             write(*,*) w(u1_), w(u2_), w(u3_), w(T_eps_)
             write(*,*) 'w(beta_metric1_), w(beta_metric2_), w(beta_metric3_), w(psi_metric_), w(alp_metric_), x(1:ndim)'
             write(*,*) w(beta_metric1_), w(beta_metric2_), w(beta_metric3_), w(psi_metric_), w(alp_metric_), x(1:ndim)
             write(*,*) 'Location at  ', Location
             call mpistop ('Psi_metric < 1 or alp_metric > 1 at this location')
          endif
         !endif
       
   }

    do i = 1, nw
    if (w( i) .ne. w( i)) then
       write(*,*) 'w( i) , i'
       write(*,*) w( i) , i
         stopflag = .true.
    endif
    enddo

    if (stopflag) then
      write(*,*)  'cons_tmp(d_), cons_tmp(s1_), cons_tmp(s2_),&
                   cons_tmp(s3_),cons_tmp(tau_)'
      write(*,*)  cons_tmp(d_), cons_tmp(s1_), cons_tmp(s2_),&
                  cons_tmp(s3_),cons_tmp(tau_)

      write(*,*) 'Location at  ', Location
      call mpistop('NaN for w in this location')
    endif

         ! Not diagonal ghostzones
    !if ( (x(1) .ge. 0.0d0 .and. x(2) < 0.0d0) .or. (x(1) < 0.0d0 .and. x(2) .ge. 0.0d0)) then
      if (eos_type == tabulated) then
         if (w(rho_) .gt. eos_rhomax) then
             write(*,*) 'w(rho_), eos_rhomax, x(1:ndim)'
             write(*,*) w(rho_), eos_rhomax, x(1:ndim)
             write(*,*) 'Location at  ', Location
           call mpistop("rho > eosmax in this location")
         endif
         if (w(ye_)  .lt. eos_yemin) then
             write(*,*) 'w(ye_), eos_yemin, x(1:ndim)'
             write(*,*) w(ye_), eos_yemin, x(1:ndim)
             write(*,*) 'Location at  ', Location
           call mpistop("ye < eosmin in this location")
         endif
      endif
    
  end subroutine Simple_check_data_correctness_pt
  !> Returns 0 in argument flag where values are ok
  subroutine prim_NaN_checker(ixI^L, ixO^L, w, x)
   include 'amrvacdef.f'

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    integer                      :: flag(ixI^S)
    integer                      :: iw, ix^D

    ! reset the value
    flag(ixO^S) = 0

    do iw=1, nw
       {do ix^D = ixO^LIM^D \}
          if ( w(ix^D,iw) /= w(ix^D,iw) ) flag(ix^D) = iw
       {end do^D&\}
    end do

    if ( any(flag(ixO^S) /= 0) ) &
       call NaN_values_error(w, x, ixI^L, ixO^L, flag)
  end subroutine prim_NaN_checker

  subroutine NaN_values_error(w, x, ixI^L, ixO^L, w_flag)
include 'amrvacdef.f'

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer, intent(in)          :: w_flag(ixI^S)
    integer                      :: ix_bad(ndim), iw

    ix_bad = maxloc(w_flag(ixO^S)) + [ ixOmin^D-1 ]

      write(*,*) "Error: NaN value of " !trim(prim_names(maxval(w_flag(ixO^S))))
!      write(*,*) "Iteration: ", it, " Time: ", t
      write(*,*) "Location: ", x({ix_bad(^D)}, :)
      write(*,*) "Cell number: ", ix_bad(:)
      do iw = 1, nw
         write(*, *) iw, ": ", &
              w({ix_bad(^D)}, iw)
      end do
      write(*,*) "Saving status at the previous time step"
     stop " crashed NaN_values_err"
  end subroutine NaN_values_error

  !> fixme: needed to be refined
  subroutine small_values_error(w, x, ixI^L, ixO^L, w_flag, subname)
include 'amrvacdef.f'

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer, intent(in)          :: w_flag(ixI^S)
    integer                      :: ix_bad(ndim), iw
    character(len=*), intent(in) :: subname

    ix_bad = maxloc(w_flag(ixO^S)) + [ ixOmin^D-1 ]

      write(*,*) "Error: small value of,  encountered when call ", subname
!      write(*,*) "Iteration: ", it, " Time: ", global_time
      write(*,*) "Location: ", x({ix_bad(^D)}, :)
      write(*,*) "Cell number: ", ix_bad(:)
      do iw = 1, nw
         write(*, *) iw, ": ", &
              w({ix_bad(^D)}, iw)
      end do
      write(*,*) "Saving status at the previous time step"
     stop " crashed small_values_err"
  end subroutine small_values_error

!  subroutine small_values_average(ixI^L, ixO^L, w, x, w_flag)
!include 'amrvacdef.f'
!
!    integer, intent(in)             :: ixI^L, ixO^L
!    integer, intent(in)             :: w_flag(ixI^S)
!    double precision, intent(inout) :: w(ixI^S, 1:nw)
!    double precision, intent(in)    :: x(ixI^S, 1:ndim)
!    integer                         :: iw, kxO^L, ix^D, i
!
!    {do ix^DB= ixO^LIM^DB\}
!
!    ! point with local failure identified by w_flag
!    if (w_flag(ix^D) /= 0) then
!      ! verify in cube with border width small_values_daverage the presence of
!      ! cells where all went ok
!      do i = 1, max(small_values_daverage, 1)
!        {kxOmin^D= max(ix^D-i, ixOmin^D);
!        kxOmax^D= min(ix^D+i, ixOmax^D);\}
!
!        ! in case cells are fine within smaller cube than 
!        ! the userset small_values_daverage: use that smaller cube
!        if (any(w_flag(kxO^S) == 0)) exit
!      end do
!
!      if (any(w_flag(kxO^S) == 0)) then
!        ! within surrounding cube, cells without problem were found
!
!        ! faulty cells are corrected by averaging here
!        ! only average those which were ok and replace faulty cells
!        do iw = 1, nw
!          if (small_values_fix_iw(iw)) then
!            w(ix^D, iw) = sum(w(kxO^S, iw), w_flag(kxO^S) == 0)&
!                 / count(w_flag(kxO^S) == 0)
!          end if
!        end do
!      else
!        write(*,*) "no cells without error were found in cube of size", & 
!             small_values_daverage
!        write(*,*) "at location:", x(ix^D, 1:ndim)
!        write(*,*) "at index:", ix^D
!        write(*,*) "w_flag(ix^D):", w_flag(ix^D)
!        write(*,*) "Saving status at the previous time step"
!     stop " crashed small_values_average"
!      end if
!    end if
!    {enddo^D&\}
!
!  end subroutine small_values_average

end module mod_small_values
