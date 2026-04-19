!> Module for eos
module mod_eos_piecewise
  use mod_eos

  implicit none
  private

  double precision            :: eos_adiab_0 = 8.9505832D-2
  double precision            :: eos_gamma_0 = 1.35692395d0

  double precision            :: eos_gamma_1 = 1.3569d0 !1.3569
  double precision            :: eos_gamma_2 = 2.664d0  !2.664
  double precision            :: eos_gamma_3 = 2.194d0  !2.194
  double precision            :: eos_adiab_1 = 8.951d-2 !0.08951
  double precision            :: eos_adiab_2 = 1.371d4  !13710
  double precision            :: eos_adiab_3 = 4.836d2  !483.6
  double precision            :: eos_rho_0   = 1.0d-30
  double precision            :: eos_rho_1   = 1.07896d-4 
  double precision            :: eos_rho_2   = 8.11905d-4 

  !============================================================
  !> piecewise polytrope related parameters
  !============================================================
  integer, parameter          :: pw_type = 0

  !> options of pw eos
  integer, parameter          :: pw_user = 0
  integer, parameter          :: pw_GNH3 = 1
  integer, parameter          :: pw_H4   = 2
  integer, parameter          :: pw_ALF2 = 3
  integer, parameter          :: pw_SLy  = 4
  integer, parameter          :: pw_APR4 = 5

  ! if use pw_user, values have to be set here
  double precision            :: pw_adiab(0:3)
  double precision            :: pw_gamma(0:3)
  double precision            :: pw_rho(0:3)
  double precision            :: pw_a(0:3)

  ! GNH3
  double precision, parameter :: GNH3_gamma(1:3) = (/ 2.664d0, 2.194d0, 2.304d0 /)
  double precision, parameter :: GNH3_rho(0:2)   = (/ 1.078957079d-4, 8.119048974d-4, 1.619963244d-3 /)
  ! H4
  double precision, parameter :: H4_gamma(1:3)   = (/ 2.909d0, 2.246d0, 2.144d0 /)
  double precision, parameter :: H4_rho(0:2)     = (/ 1.43735d-4, 8.119048974d-4, 1.619963244d-3 /)
  ! ALF2
  double precision, parameter :: ALF2_gamma(1:3) = (/ 4.070d0, 2.411d0, 1.890d0 /)
  double precision, parameter :: ALF2_rho(0:2)   = (/ 3.1534d-4, 8.1147d-4, 1.6191d-3 /)
  ! SLy
  double precision, parameter :: SLy_gamma(1:3)  = (/ 3.005d0, 2.988d0, 2.851d0 /)
  double precision, parameter :: SLy_rho(0:2)    = (/ 2.3674d-4, 8.1147d-4, 1.6191d-3 /)
  ! APR4
  double precision, parameter :: APR4_gamma(1:3) = (/ 2.830d0, 3.445d0, 3.348d0 /)
  double precision, parameter :: APR4_rho(0:2)   = (/ 2.4479d-4, 8.1147d-4, 1.6191d-3 /)

  public :: eos_piecewise_activate
  public :: piecewise_get_eps

contains

  subroutine eos_piecewise_activate()
    use mod_global_parameters
    double precision            :: pw_prs(0:3)
    double precision            :: pw_eps(0:3)
    double precision            :: pw_h(0:3)
    integer                     :: n, i
    character(len=name_len)     :: pw_type = ""

    namelist /eos_piecewise_list/ pw_type, eos_rho_0, eos_rho_1, eos_rho_2, &
                eos_gamma_0, eos_gamma_1, eos_gamma_2, eos_gamma_3, &
                eos_adiab_0, eos_adiab_1, eos_adiab_2, eos_adiab_3

    ! Read this module's parameters from a file
    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status="old")
       read(unitpar, eos_piecewise_list, end=111)
111    close(unitpar)
    end do

    eos_type = piecewise

    pw_adiab(0) = eos_adiab_0
    pw_gamma(0) = eos_gamma_0
    if (pw_type == "user_define") then
       pw_adiab(1) = eos_adiab_1
       pw_adiab(2) = eos_adiab_2
       pw_adiab(3) = eos_adiab_3
       pw_gamma(1) = eos_gamma_1
       pw_gamma(2) = eos_gamma_2
       pw_gamma(3) = eos_gamma_3
       pw_rho(0)   = eos_rho_0
       pw_rho(1)   = eos_rho_1
       pw_rho(2)   = eos_rho_2
    else
       select case (pw_type)
       case ('GNH3') 
          pw_gamma(1:3) = GNH3_gamma
          pw_rho(0:2)   = GNH3_rho
       case ('H4') 
          pw_gamma(1:3) = H4_gamma
          pw_rho(0:2)   = H4_rho
       case ('ALF2') 
          pw_gamma(1:3) = ALF2_gamma
          pw_rho(0:2)   = ALF2_rho
       case ('SLy') 
          pw_gamma(1:3) = SLy_gamma
          pw_rho(0:2)   = SLy_rho
       case ('APR4') 
          pw_gamma(1:3) = APR4_gamma
          pw_rho(0:2)   = APR4_rho
       case default
         error stop "this pw_type is not supported"
       end select
    end if
    ! whatever larger the pw_rho2 is fine
    pw_rho(3) = maxval(pw_rho(0:2)) * 5.0d0

    ! get K_i
    do i = 1, 3
       pw_adiab(i) = pw_adiab(i-1)*pw_rho(i-1)**pw_gamma(i-1) &
                     / pw_rho(i-1)**pw_gamma(i)
    end do

    ! get a_i and e_i
    pw_eps(0) = pw_rho(0) &
      + pw_adiab(0) * pw_rho(0)**pw_gamma(0) &
       / (pw_gamma(0)-1.0d0)
    pw_a(0) = 0.0d0
    do i = 1, 3
      pw_a(i) = pw_eps(i-1) / pw_rho(i-1) - 1.0d0 &
        - pw_adiab(i) / (pw_gamma(i)-1.0d0) &
        * pw_rho(i-1)**(pw_gamma(i)-1.0d0)
      pw_eps(i) = (1.0d0+pw_a(i)) * pw_rho(i) &
        + pw_adiab(i) * pw_rho(i)**pw_gamma(i) &
         / (pw_gamma(i)-1.0d0)
    end do

    ! update the rest of the variables
    pw_prs = pw_adiab * pw_rho**pw_gamma
    pw_h = 1.0d0 + pw_a + pw_gamma/(pw_gamma-1.0d0) & 
            * pw_adiab * pw_rho**(pw_gamma-1.0d0)
    ! note that the pw_eps above is actually
    ! the energy density in Read et. al. 2009
    ! which is (1+eps)*rho in our notation,
    ! where eps is the internal energy density.
    ! here we transform it back to ours
    pw_eps = pw_eps/pw_rho - 1.0d0

    ! check if the parameters are make sences
    if (any(pw_gamma <= 0.0d0)) call mpistop("Error: some pw_gamma <= 0")
    if (any(pw_adiab <  0.0d0)) call mpistop("Error: some pw_adiab < 0")

    eos_setup_atmo_eos => piecewise_setup_atmo_eos
    eos_get_pressure   => piecewise_get_pressure
    eos_get_e          => piecewise_get_e
    eos_get_cs2        => piecewise_get_cs2
    eos_get_all        => piecewise_get_all

  end subroutine eos_piecewise_activate

  subroutine piecewise_setup_atmo_eos()
    atmo_adiab = pw_adiab(0)
    atmo_gamma = pw_gamma(0)
  end subroutine piecewise_setup_atmo_eos

  integer function pw_i(rho_in)
    implicit none
    double precision, intent(in) :: rho_in
    if ( rho_in < pw_rho(0) ) then
       pw_i = 0
    else if ( rho_in < pw_rho(1) ) then
       pw_i = 1
    else if ( rho_in < pw_rho(2) ) then
       pw_i = 2
    !else if ( pw_rho(2) <= rho_in ) then
    else
       pw_i = 3
    end if
  end function

  subroutine piecewise_get_pressure(prs,rho,e,ye)
    implicit none
    double precision, intent(inout) :: prs
    double precision, intent(in)    :: rho, e, ye
    integer                         :: i
    if (rho<small_rho_thr) then
       call atmo_get_pressure(prs,rho,e,ye)
       return
    end if
    i = pw_i(rho)
    prs = pw_adiab(i) * rho**pw_gamma(i)
  end subroutine piecewise_get_pressure

  subroutine piecewise_get_e(rho,ye,eps,e)
    implicit none
    double precision, intent(in)    :: rho, ye, eps
    double precision, intent(inout) :: e
    if (rho<small_rho_thr) then
       call atmo_get_e(rho,ye,eps,e)
       return
    end if
    e = rho * eps
  end subroutine piecewise_get_e

  subroutine piecewise_get_cs2(cs2,rho,e,ye)
    implicit none
    double precision, intent(inout) :: cs2
    double precision, intent(in)    :: rho, e, ye
    double precision                :: prs, rhoh
    integer                         :: i
    if (rho<small_rho_thr) then
       call atmo_get_cs2(cs2,rho,e,ye)
       return
    end if
    i = pw_i(rho)
    prs = pw_adiab(i) * rho**pw_gamma(i)
    rhoh = rho + e + prs
    ! use prs and h
    cs2 = pw_gamma(i) * prs / rhoh
  end subroutine piecewise_get_cs2

  !> Get all eos var from (rho, e, ye)
  subroutine piecewise_get_all(rho,e,ye,rhoeps,prs,cs2,ent)
    implicit none
    double precision, intent(in)              :: rho, e, ye
    double precision, intent(inout), optional :: rhoeps, prs, cs2, ent
    if (present(rhoeps)) rhoeps = e
    if (present(prs)) call eos_get_pressure(prs,rho,e,ye)
    if (present(cs2)) call eos_get_cs2(cs2,rho,e,ye)
    if (present(ent)) ent = 0.0d0
  end subroutine piecewise_get_all

  subroutine piecewise_get_eps(prs,rho,eps)
    implicit none
    double precision, intent(in)    :: prs
    double precision, intent(in)    :: rho
    double precision, intent(inout) :: eps
    integer                         :: i

    if (rho<small_rho_thr) then
       call atmo_get_eps_one_grid(prs,rho,eps)
       return
    end if
    i = pw_i(rho)
    eps  = pw_a(i) + prs / rho / (pw_gamma(i)-1.0d0)
  end subroutine piecewise_get_eps

end module mod_eos_piecewise
