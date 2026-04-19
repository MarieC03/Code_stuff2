!> Module for eos
module mod_eos_idealgas


  implicit none
  public


contains

  !> Read this module's parameters from a file
!  subroutine eos_idealgas_read_params(files)

!  end subroutine eos_idealgas_read_params

  subroutine eos_idealgas_activate()
    use mod_eos
    include 'amrvacdef.f'

    eos_type = idealgas
    ! check if the parameters are make sences
    if (eos_gamma <= 0.0d0) call mpistop ("Error: eos_gamma <= 0")
    if (eos_adiab < 0.0d0) call mpistop  ("Error: eos_adiab < 0")

    eos_get_pressure_one_grid         => idealgas_get_pressure_one_grid
    eos_get_eps_one_grid              => idealgas_get_eps_one_grid
    eos_get_cs2_one_grid              => idealgas_get_cs2_one_grid


    if (mype == 0 ) then
        write(*,*) 'eos_gamma:', eos_gamma 
        write(*,*) 'eos_adiab:', eos_adiab
        write(*,*) 'small_rho:', small_rho
        write(*,*) 'small_press:', small_press
        write(*,*) 'small_eps:', small_eps
        print*, '****************************************'
        print*, '-------Idealgas EOS activated----------'
        print*, '****************************************'
    endif

!write(*,*) "activated idealgas law!"
!write(*,*) "eos_gamma, eos_adiab"
!write(*,*) eos_gamma, eos_adiab
  call eos_check

  end subroutine eos_idealgas_activate

  subroutine idealgas_get_pressure_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    include 'amrvacdef.f'
  !  implicit none
    
    double precision, intent(inout) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: ye
    double precision, intent(inout), optional :: temp

    if (.not. old_bhac_safety) then
      if (rho<small_rho_thr) then
         call atmo_get_pressure_one_grid(prs,rho,eps)
         return
      end if
    endif

    prs = ( eos_gamma - 1.0d0 ) * rho * eps

  end subroutine idealgas_get_pressure_one_grid

  subroutine idealgas_get_eps_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    include 'amrvacdef.f'
    !implicit none
    
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(inout) :: eps
    double precision, intent(in), optional :: temp, ye

    if (.not. old_bhac_safety) then
      if (rho<small_rho_thr) then
         call atmo_get_eps_one_grid(prs,rho,eps)
         return
      end if
    endif
    eps = prs / rho / ( eos_gamma - 1.0d0 )

  end subroutine idealgas_get_eps_one_grid

  subroutine idealgas_get_cs2_one_grid(cs2,rho,eps,temp,ye)

    use mod_eos
    include 'amrvacdef.f'
    !implicit none
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: ye
    double precision, intent(inout), optional :: temp

    double precision             :: prs
    double precision             :: dpde
    double precision             :: dpdrho
    double precision             :: rhoeps

    if (.not. old_bhac_safety) then
      if (rho<small_rho_thr) then
       call atmo_get_cs2_one_grid(cs2,rho,eps)
         return
      end if
    endif

    !prs = ( eos_gamma - 1.0d0 ) * rho * eps
    !dpde = (eos_gamma - 1.0d0 ) * rho
    !dpdrho = (eos_gamma - 1.0d0 ) * eps
    !cs2= dpdrho+dpde*prs/rho**2
    !cs2 = cs2 / ( (1.0d0 + prs/rho) + eps )

    rhoeps = rho * eps
    prs    = ( eos_gamma - 1.0d0 ) * rhoeps
    cs2    = eos_gamma * prs / ( rho + rhoeps + prs )


  end subroutine idealgas_get_cs2_one_grid


end module mod_eos_idealgas
