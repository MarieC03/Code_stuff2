!> Module for eos
module mod_eos_polytrope


  implicit none
  public


contains

  !> Read this module's parameters from a file
  subroutine eos_polytrope_read_params()
!    character(len=*), intent(in) :: files(:)
!    integer                      :: n

!    namelist /eos_polytrope_list/ eos_gamma, eos_adiab

!    do n = 1, size(files)
!       open(unitpar, file=trim(files(n)), status="old")
!       read(unitpar, eos_polytrope_list, end=111)
!111    close(unitpar)
!    end do

  end subroutine eos_polytrope_read_params

  subroutine eos_polytrope_activate()
    use mod_eos
    include 'amrvacdef.f'

    eos_type = polytrope
    ! check if the parameters are make sences
    if (eos_gamma <= 0.0d0) call mpistop ("Error: eos_gamma <= 0")
    if (eos_adiab < 0.0d0) call mpistop  ("Error: eos_adiab < 0")

    eos_get_pressure_one_grid         => polytrope_get_pressure_one_grid
    eos_get_eps_one_grid              => polytrope_get_eps_one_grid
    eos_get_cs2_one_grid              => polytrope_get_cs2_one_grid

    call eos_check
    if (mype == 0 ) then 
        write(*,*) 'eos_gamma:', eos_gamma
        write(*,*) 'eos_adiab:', eos_adiab
        print*, '****************************************'
        print*, '-------Polytrope EOS activated----------'
        print*, '****************************************'
    endif

  end subroutine eos_polytrope_activate

  subroutine polytrope_get_pressure_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(inout) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: ye
    double precision, intent(inout), optional :: temp

   ! if (rho<small_rho_thr) then
   !    call atmo_get_pressure_one_grid(prs,rho,eps)
   !    return
   ! end if
    prs = eos_adiab * rho**eos_gamma

  end subroutine polytrope_get_pressure_one_grid

  subroutine polytrope_get_eps_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in), optional :: temp, ye
    double precision, intent(inout) :: eps

   ! if (rho<small_rho_thr) then
   !    call atmo_get_eps_one_grid(prs,rho,eps)
   !    return
   ! end if
    !eps = prs / rho / ( eos_gamma - 1.0d0 )
! expression without using prs to cal
    eps = eos_adiab * rho**( eos_gamma - 1.0d0 ) &
          / ( eos_gamma - 1.0d0 )

  end subroutine polytrope_get_eps_one_grid

  subroutine polytrope_get_cs2_one_grid(cs2,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: ye
    double precision, intent(inout), optional :: temp

    double precision             :: prs, h

   ! if (rho<small_rho_thr) then
   !    call atmo_get_cs2_one_grid(cs2,rho,eps)
   !    return
   ! end if

    prs = eos_adiab * rho**eos_gamma
    h = 1.0d0 + eps + prs/rho

    ! use prs and h
    cs2= eos_gamma * prs / ( rho * h )

    ! use prs only
    !cs2= eos_gamma*( eos_gamma - 1.0d0 ) * prs &
    !     /( rho*( eos_gamma - 1.0d0 ) + eos_gamma*prs )

    ! strictly require cs2 is physical
    cs2 = max(0.d0, min(1.d0, cs2))

  end subroutine polytrope_get_cs2_one_grid

end module mod_eos_polytrope
