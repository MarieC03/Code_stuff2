!> Module for eos
module mod_eos

  use mod_eos_tabulated_parameters
  public


  ! code test
  !integer                           :: EOS_call_iteration
  !integer                           :: c2p_rootfinding_call_iteration
  !integer                           :: c2p_rootfinding_call_iteration_masfunconly
  integer :: total_cell_prev
  logical :: first_time = .true.

  integer, parameter                :: polytrope = 1
  integer, parameter                :: idealgas = 2
  integer, parameter                :: hybrid = 3
  integer, parameter                :: tabulated = 4

  
  character(128)                    :: atmo_type
  integer, parameter                :: beta_eqm = 5

  double precision                  :: atmo_gamma = 2.0d0
  double precision                  :: atmo_adiab = 1.0d2

  ! these two for idealgas and polytrope
  double precision                  :: eos_gamma = 2.0d0
  double precision                  :: eos_gamma_th = 1.0d0
  double precision                  :: eos_adiab = 1.0d2

  double precision, parameter       :: bigdouble_eos   = 1.0d+99 
  double precision, parameter       :: smalldouble_eos = 1.0d-99
  ! min-max values for tabulated eos:
  double precision                  :: eos_rhomin  = smalldouble_eos
  double precision                  :: eos_rhomax  = bigdouble_eos
  double precision                  :: eos_yemin   = smalldouble_eos
  double precision                  :: eos_yemax   = bigdouble_eos
  double precision                  :: eos_tempmin = smalldouble_eos
  double precision                  :: eos_tempmax = bigdouble_eos
  double precision                  :: eos_epsmin  = 0.0d0
  double precision                  :: eos_epsmax  = bigdouble_eos
  double precision                  :: eos_hmin    = 1.0d0   ! for tabeos, find from the table
  double precision                  :: eos_prsmin  = 0.0d0
  double precision                  :: eos_prsmax  = bigdouble_eos

  logical                           :: corrupted_compose = .false.

  !> The smallest allowed density
  integer                      :: atmo_eos_type = 1


  double precision             :: small_rho_fac = 0.1d0
!fixme
!need to change
  double precision             :: small_rho_thr = 0.0d0 ! = 1.0d-19
  double precision             :: small_rho = smalldouble_eos
  !> The smallest allowed eps
  double precision             :: small_eps = smalldouble_eos
  !> The smallest allowed press
  double precision             :: small_press = smalldouble_eos
  !> The smallest allowed ethalpy
  double precision             :: small_h = smalldouble_eos
  !> The biggest allowed electron fraction
  double precision             :: big_ye = bigdouble_eos
  !> The smallest allowed temp
  double precision             :: small_temp = smalldouble_eos

  !> The smallest allowed conserved density D
  double precision             :: small_D = smalldouble_eos
  !> The smallest allowed conserved variables tau ~ eps * rho
  double precision             :: small_tau = smalldouble_eos

  double precision             :: atmo_cs2 = 0.0d0


! for tabeos_beta_eqm_atm 
  double precision            :: ye_beta_eqm_atm
  double precision            :: eps_beta_eqm_atm
  double precision            :: press_beta_eqm_atm
  double precision            :: temp_beta_eqm_atm

  ! constants
  double precision, parameter :: ggrav              = 6.673d-8
  double precision, parameter :: c_cgs              = 2.99792458d10


  double precision, parameter :: mev_to_erg         = 1.60217733d-6
  double precision, parameter :: erg_to_mev         = 6.24150636d5
  double precision, parameter :: rho_gf             = 1.61930347d-18
  double precision, parameter :: press_gf           = 1.80171810d-39
  double precision, parameter :: eps_gf             = 1.11265006d-21
  double precision, parameter :: time_gf            = 2.03001708d+05
  double precision, parameter :: mass_gf            = 5.02765209d-34
  double precision, parameter :: length_gf          = 6.77140812d-06
  double precision, parameter :: energy_gf          = 5.59424238d-55
  double precision, parameter :: lum_gf             = 2.7556091d-60
  double precision, parameter :: MeV_gf      = MeV_to_erg * energy_gf

  double precision, parameter :: amu_cgs            = 1.66053873d-24
  double precision            :: massn_cgs          = 1.674927211d-24
  double precision, parameter :: amu_mev            = 931.49432d0
  double precision, parameter :: kb_erg             = 1.380658d-16
  double precision, parameter :: kb_mev             = 8.61738568d-11
  double precision, parameter :: temp_mev_to_kelvin = 1.1604447522806d10
  double precision, parameter :: planck             = 6.626176d-27
  double precision, parameter :: avo                = 6.0221367d23
  double precision, parameter :: hbarc_mevcm        = 1.97326966d-11
  double precision, parameter :: cm3_to_fm3         = 1.0d39
  double precision, parameter :: mp_mev             = 938.28             ! rest mass    in Mev/c^2
  double precision, parameter :: mn_mev             = 939.57             !
  double precision, parameter :: Qnp                = 1.293548d0      ! neutron proton mass difference in MeV


  !public

  !> describing the eos type of the simulation
  integer                              :: eos_type = -1
  character(128)                       :: eos_type_input = ""

  procedure(sub_get_pressure_one_grid), pointer        :: eos_get_pressure_one_grid        => null()
  procedure(sub_get_eps_one_grid), pointer             :: eos_get_eps_one_grid             => null()
  procedure(sub_get_cs2_one_grid), pointer             :: eos_get_cs2_one_grid             => null()

  procedure(sub_get_eps_range), pointer                :: eos_get_eps_range                => null()
  procedure(sub_get_temp_one_grid), pointer            :: eos_get_temp_one_grid            => null()
  procedure(sub_get_all_beta_eqm_one_grid), pointer    :: eos_get_all_beta_eqm_one_grid    => null()
  procedure(sub_eps_get_all_one_grid), pointer         :: eos_eps_get_all_one_grid         => null()
  procedure(sub_temp_get_all_one_grid), pointer        :: eos_temp_get_all_one_grid        => null()
!  procedure(sub_tabulated_interpolation), pointer      :: intep3d                          => null()
!  procedure(sub_tabulated_interpolation_many), pointer :: intep3d_many                     => null()

  abstract interface

     subroutine sub_get_eps_range(rho,epsmin,epsmax,ye)
       double precision, intent(in) :: rho
       double precision, intent(out):: epsmin, epsmax
       double precision, intent(in), optional :: ye
     end subroutine sub_get_eps_range

     subroutine sub_get_pressure_one_grid(prs,rho,eps,temp,ye)
       double precision, intent(inout) :: prs
       double precision, intent(in) :: rho
       double precision, intent(in) :: eps
       double precision, intent(in), optional :: ye
       double precision, intent(inout), optional :: temp
     end subroutine sub_get_pressure_one_grid
   
     subroutine sub_get_eps_one_grid(prs,rho,eps,temp,ye)
       double precision, intent(in) :: prs
       double precision, intent(in) :: rho
       double precision, intent(inout) :: eps
       double precision, intent(in), optional :: temp, ye
     end subroutine sub_get_eps_one_grid

     subroutine sub_get_cs2_one_grid(cs2,rho,eps,temp,ye)
       double precision, intent(inout) :: cs2
       double precision, intent(in) :: rho
       double precision, intent(in) :: eps
       double precision, intent(in), optional :: ye
       double precision, intent(inout), optional :: temp
     end subroutine sub_get_cs2_one_grid

     subroutine sub_get_temp_one_grid(rho,eps,temp,ye)
       double precision, intent(in) :: rho
       double precision, intent(in) :: ye
       double precision, intent(inout) :: temp, eps
     end subroutine sub_get_temp_one_grid

     subroutine sub_get_all_beta_eqm_one_grid(rho,temp,ye,eps,prs)
       double precision, intent(in)    :: rho, temp
       double precision, intent(inout) :: ye
       double precision, intent(inout), optional :: eps, prs
     end subroutine sub_get_all_beta_eqm_one_grid

     subroutine sub_eps_get_all_one_grid(rho,eps,ye,temp,prs,ent,cs2,dedt,&
                 dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)

       double precision, intent(inout) :: rho, eps, ye
       double precision, intent(inout), optional :: ent, prs, temp, cs2, dedt
       double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
       double precision, intent(inout), optional :: mu_e,mu_n,mu_p,muhat,munu

     end subroutine sub_eps_get_all_one_grid

     subroutine sub_temp_get_all_one_grid(rho,temp,ye,eps,prs,ent,cs2,dedt,&
                 dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)
       double precision, intent(inout) :: rho, eps, temp
       double precision, intent(in) :: ye

       double precision, intent(inout), optional :: prs,ent,cs2,dedt,&
                   dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu

     end subroutine sub_temp_get_all_one_grid


  end interface

contains

  subroutine eos_check
    {#IFNDEF UNIT_TESTS
    if (eos_type == -1) call mpistop("Error: no eos module is loaded")
    ! Checks whether the required physics methods have been defined
    if (.not. associated(eos_get_pressure_one_grid)) &
         call mpistop("Error: eos_get_pressure_one_grid not defined")

    if (.not. associated(eos_get_eps_one_grid)) &
         call mpistop("Error: eos_get_eps_one_grid not defined")

    if (.not. associated(eos_get_cs2_one_grid)) &
         call mpistop("Error: eos_get_cs2_one_grid not defined")
    }
    if (.not. associated(eos_get_eps_range)) then
      {#IFNDEF UNIT_TESTS
       if (eos_type == tabulated) &
         call mpistop("Error: eos_get_eps_range is not defined but using tabulated eos")
         }
       eos_get_eps_range => unlimited_eps_range
    end if
  end subroutine eos_check

  subroutine eos_atmo_activate()

    integer                                 :: n
!    character(len=128), intent(in)          :: atmo_type
!    double precision , intent(in), optional :: atmo_gamma_in, atmo_adiab_in
 
     call eos_check
!    if (present(atmo_gamma_in)) then 
!        atmo_gamma = atmo_gamma_in
!    endif
!   
!    if (present(atmo_adiab_in)) then
!        atmo_adiab = atmo_adiab_in 
!    endif

    ! check if the parameters are make sences

    select case (atmo_type)
    case ('polytrope','idealgas')
        atmo_eos_type = polytrope
        {#IFNDEF UNIT_TESTS
        if (atmo_gamma <= 0.0d0) call mpistop ("Error: atmo_gamma <= 0 in polytrope atmo")
        if (atmo_adiab < 0.0d0) call mpistop  ("Error: atmo_adiab < 0 in polytrope atmo")
        }
    case ('beta_eqm')
        atmo_eos_type = beta_eqm
    case default
      error stop "this eos_type is not supported"
    end select

!    small_rho_thr = eos_rhomin
  end subroutine eos_atmo_activate

  subroutine eos_initialize_atmo()
    {#IFDEF UNIT_TESTS
    use amrvacdef
    }
    {#IFNDEF UNIT_TESTS
    include 'amrvacdef.f'
    }
!    double precision, intent(in)           :: small_rho_thr_in, small_rho_in, small_rho_fac_in
    double precision, parameter :: idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))   !this just like gr1d
    double precision, parameter :: idealgamma = 1.66666666666d0
    double precision :: eps_min, eps_max, ye_tmp
    {#IFNDEF UNIT_TESTS
    if (eos_type == -1) call mpistop("Error: atmosphere can only be initialized after activitating eos")
    }
!    small_rho     = small_rho_in
!    small_rho_fac = small_rho_fac_in 
!    small_rho_thr = small_rho_thr_in * small_rho
    small_rho     = small_rho_thr * small_rho_fac
   
    select case (atmo_eos_type)
    case(polytrope)
    
        ! make sure atmosphere value is larger or equal to smalldouble_eos
        if (small_rho <= smalldouble_eos) then
           small_rho     = smalldouble_eos
           small_rho_thr = small_rho / small_rho_fac
        end if

       ! small_press = smalldouble_eos
       ! small_eps = smalldouble_eos
       ! small_press = atmo_adiab * small_rho**atmo_gamma
       ! small_eps = atmo_adiab * small_rho**(atmo_gamma - 1.0d0)/(atmo_gamma - 1.0d0)

       small_press = atmo_adiab * small_rho**atmo_gamma
       small_eps   = atmo_adiab * small_rho**( atmo_gamma - 1.0d0 ) &
                     / ( atmo_gamma - 1.0d0 )
       atmo_cs2 = atmo_gamma * ( atmo_gamma - 1.0d0 ) * small_press &
            / ( small_rho * ( atmo_gamma - 1.0d0 ) + atmo_gamma * small_press )


       ! make sure all the allowed min value are equal to atmosphere
       eos_rhomin = small_rho
       eos_epsmin = small_eps
       small_D = small_rho_thr
       small_tau = small_rho_thr * small_eps

    case(beta_eqm)
       small_temp    = max(small_temp, eos_tempmin)
       ! beta eqm atmo
       call eos_get_all_beta_eqm_one_grid(small_rho, small_temp,  &
                                          ye_beta_eqm_atm, eps=eps_beta_eqm_atm, prs=press_beta_eqm_atm)
 
       if (corrupted_compose) then ! eos_rhomin, eos_tempmin doesnot mean lowest eps in some compose tables
            ! so, by inputting small_rho_thr to modify the eos_epsmin, hmin
            ! fixme: hmin may be not the min for this eos_epsmin, what should I do
            eos_epsmax = 10.0d0**maxval(eos_tables(:,:,:,i_logenergy)) - energy_shift
 
            call eos_get_all_beta_eqm_one_grid(small_rho, eos_tempmin, ye_tmp)
            call eos_temp_get_all_one_grid(small_rho_thr, eos_tempmin, ye_tmp,&
                                 eos_epsmin, prs=eos_prsmin)
            eos_hmin = 1.0d0 + eos_epsmin + eos_prsmin/eos_rhomax
            eos_hmin = min(1.0d0, eos_hmin)
            !stop 'warning the EOS you are using has corrupted data, &
            !     if you still want to go, pls comment this line.'
       endif

       small_eps   = max(eos_epsmin, eps_beta_eqm_atm)
       small_press = press_beta_eqm_atm
       big_ye      = min(ye_beta_eqm_atm, eos_yemax)
       !atmo_cs2    = 0.0d0  !fixme: pls find this by beta_eqm
       call eos_get_cs2_one_grid(atmo_cs2,small_rho,small_eps,small_temp,big_ye)
       small_D     = small_rho_thr
       small_tau   = -1.0d10   !  should not be higher than eos_epsmin
       ! check small values within the bounds of eos
       {#IFNDEF UNIT_TESTS
          if (small_rho .lt. eos_rhomin .or. small_rho_thr .lt. eos_rhomin) then
             call mpistop("small_rho .lt. eos_rhomin .or. small_rho_thr .lt. eos_rhomin in initialize_atmo")
          endif
          if (big_ye .gt. eos_yemax) then
             call mpistop("(big_ye .gt. eos_yemax in initialize_atmo)")
          endif
          if (small_eps .lt. eos_epsmin) then
             call mpistop("(small_eps .lt. eos_epsmin in initialize_atmo)")
          endif
          if (small_temp .lt. eos_tempmin) then
             call mpistop("(small_temp .lt. eos_tempmin in initialize_atmo)")
          endif
          if (small_press .le. 0.0d0) then
             call mpistop("(small_press .le. 0.0d0 in initialize_atmo)")
          endif
        }
    end select
 

   if (mype == 0) then
    write(*,*) "################################################"
    write(6,*) "Not through set_eps_max_min"
    write(6,*) "Done reading eos tables ,  all in code unit"
    write(*,*)  eos_yemin, eos_yemax, "eos_yemin  and max"
    write(*,*)  eos_tempmin, eos_tempmax, "eos_tempmin  and max"
    write(*,*)  eos_rhomin, eos_rhomax, "eos_rhomin  and max"
    write(*,*)  eos_epsmin, eos_epsmax, "eos_epsmin  and max"
    write(*,*)  eos_hmin,            'eos_hmin'
    write(*,*) "################################################"

    print*,'About atmosphere: '
    print*,'atmo type: ', atmo_type
    write(*,*) small_rho_thr, 'small_rho_thr'
    write(*,*) small_rho, 'small_rho'
    write(*,*) small_eps, 'small_eps'
    write(*,*) small_press, 'small_press'
    write(*,*) big_ye, 'big_ye'
    write(*,*) small_temp, 'small_temp'
    write(*,*) atmo_cs2, 'atmo_cs2' 
    write(*,*) energy_shift, 'energy shift'
    write(*,*) "################################################"
    write(*,*) 'small_D'
    write(*,*) 'small_tau'
    write(*,*) small_D
    write(*,*) small_tau
   endif



  end subroutine eos_initialize_atmo

  subroutine atmo_get_pressure_one_grid(prs,rho,eps)
    implicit none
    double precision, intent(inout) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    prs = small_press
  end subroutine atmo_get_pressure_one_grid

  subroutine atmo_get_h_one_grid(h,prs,rho,eps)
    !specific enthaply atmo
    implicit none
    double precision, intent(in) :: prs, rho, eps
    double precision, intent(inout) :: h
    h = small_h
  end subroutine atmo_get_h_one_grid

  subroutine atmo_get_eps_one_grid(prs,rho,eps)
    implicit none
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(inout) :: eps
    eps = small_eps
  end subroutine atmo_get_eps_one_grid

  ! tabeos, atmo cs2 is not zero: fixme
  subroutine atmo_get_cs2_one_grid(cs2,rho,eps)
    implicit none
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision             :: prs
    double precision             :: dpde, h
    double precision             :: dpdrho
     cs2 = atmo_cs2
  end subroutine atmo_get_cs2_one_grid

  subroutine unlimited_eps_range(rho, epsmin, epsmax, ye)
     double precision, intent(in)    :: rho
     double precision, intent(out)   :: epsmin, epsmax
     double precision, intent(in), optional :: ye
     ! if it is not using tabulated eos
     epsmin = small_eps
     epsmax = bigdouble_eos
  end subroutine unlimited_eps_range

end module mod_eos
