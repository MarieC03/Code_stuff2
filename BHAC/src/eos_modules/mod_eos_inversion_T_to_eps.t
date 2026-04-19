!> Module for eos
module mod_eos_inversion_T_to_eps
  use mod_eos_tabulated_parameters
  use mod_eos
  implicit none
  public

  double precision, allocatable :: logeps_table(:)
  double precision, allocatable :: eos_eps_tables(:,:,:,:)
  integer, parameter            :: nvars_eps = 1
  double precision              :: logepsmin
  double precision              :: logepsmax
 
  integer                       :: neps
  contains
  
  subroutine activate_inversion_T_to_eps
    use mod_eos_interpolation
    implicit none
    double precision :: logrho, logeps, logtemp, ye, logeps_bins
    integer          :: irho, itemp, iye, ieps
    double precision :: rho, dummy, temp, eps

    logepsmin = log10(eos_epsmin+energy_shift)
    logepsmax = log10(eos_epsmax+energy_shift)
    neps = ntemp

    write(*,*) 'logepsmin, logepsmax'
    write(*,*) logepsmin, logepsmax
!    if (mype == 0) then
!       write(*,*) 'You activated the inversion of tables from T to eps'
!       write(*,*) 'It is being constructed'
!    endif 

    allocate(eos_eps_tables(nrho,neps,nye,nvars_eps))
    allocate(logeps_table(neps))

    logeps_table(1) = logepsmin 
    logeps_bins     = (logepsmax - logepsmin)/neps 
    do ieps = 2, neps
       logeps_table(ieps) = logeps_table(ieps-1) + logeps_bins
    enddo

    do ieps = 1, neps
      write(*,*) logeps_table(ieps), ieps, 'logeps_table(ieps) ieps'
    enddo

    
!    do irho = 1, nrho
!    do iye  = 1, nye 
!    do itemp = 1, ntemp  ! neps = ntemp
!       rho = 1.6d-4
!       ye  = 0.2d0
!       temp = 100.0d0
!       call eos_get_eps_one_grid(dummy,rho,eps,temp,ye)
!       write(*,*) 'rho, ye, temp, eps'
!       write(*,*) rho, ye, temp, eps
!stop
!    enddo
!    enddo
!    enddo

    do irho = 1, nrho
    do iye  = 1, nye 
    do ieps = 1, neps  ! neps = ntemp

       logrho   = log10(1.6d-4)
       !logrho   = logrho_table(irho)
!       log_temp = logtemp_table(itemp)
       ye       = 0.2d0
       !ye       = ye_table(iye)
       !logeps   = logeps_table(100)
       logeps   = log10(0.609527d0 + energy_shift)
 

!       call rootfinding_brent(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)

       call get_T_one_grid_from_inversion(logrho, logeps, ye, logtemp)
!    
!       call intep3d(log_rho, log_temp, ye, &
!              eps, eos_tables(:,:,:,i_logenergy), &
!              nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
!    
!       eps = dlog10(eps) - energy_shift
!       eps_table(itemp) = eps
    enddo
    enddo
    enddo

stop
  end subroutine activate_inversion_T_to_eps

  subroutine get_T_one_grid_from_inversion(logrho, logeps, Ye, logtemp)
    use mod_eos
    use mod_rootfinding
    use mod_eos_interpolation
    double precision, intent(in)      :: logrho, logeps, Ye
    double precision, intent(inout)   :: logtemp
    integer                           :: error_code 
    double precision                  :: logtemp_min, logtemp_max
    error_code = -1

    logtemp_min = logtemp_table(1)
    logtemp_max = logtemp_table(ntemp)

    call rootfinding_brent(logtemp, logtemp_min, logtemp_max, 1.0d-15, 100, error_code, func_eps_of_temp1)
    select case (error_code)
    case (-1) ! nothing happened
       call mpistop("Error in Func_eps_of_T: " // &
           "have you ever attemp to find the root?")
    case (2) ! z is NaN
       call mpistop("Error in Func_eps_of_T: NaN")
    case (1)
       call mpistop("Fail to find the root for Func_eps_of_T")
    case(3)
    write(*,*) 'func_eps_of_temp1(logtemp_min), func_eps_of_temp1(logtemp_max)'
    write(*,*) func_eps_of_temp1(logtemp_min), func_eps_of_temp1(logtemp_max)
    write(*,*) 'logrho, logeps, ye, logtemp'
    write(*,*) logrho, logeps, ye, logtemp
       call mpistop("Cannot bracket the root for Func_eps_of_T")
    end select
    write(*,*) 'func_eps_of_temp1(logtemp_min), func_eps_of_temp1(logtemp_max)'
    write(*,*) func_eps_of_temp1(logtemp_min), func_eps_of_temp1(logtemp_max)
    write(*,*) 'logrho, logeps, ye, logtemp'
    write(*,*) logrho, logeps, ye, logtemp
stop

    contains
    double precision function func_eps_of_temp1(log_temp_in)
      double precision, intent(in) :: log_temp_in
      call intep3d(logrho, log_temp_in, ye, &
               func_eps_of_temp1, eos_tables(:,:,:,i_logenergy), &
                nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
      func_eps_of_temp1 = logeps - func_eps_of_temp1
      !func_eps_of_temp1 = 1.0d0 - func_eps_of_temp1 / logeps
    end function

    double precision function func_eps_of_temp2(log_temp_in)
      double precision, intent(in) :: log_temp_in
      call intep3d(logrho, log_temp_in, ye, &
               func_eps_of_temp2, eos_tables(:,:,:,i_logenergy), &
                nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
      func_eps_of_temp2 = 1.0d0 - func_eps_of_temp2 / logeps
    end function
   

  end subroutine get_T_one_grid_from_inversion

end module mod_eos_inversion_T_to_eps
