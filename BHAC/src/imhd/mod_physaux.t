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
! mod_physaux to perform some useful computations, e.g. obtain b2, the square 
! of the co-moving magnetic field.
! You could also make a subroutine to calculate the current and other things...
! Oliver Porth, 2016-01-27
!=============================================================================
module mod_physaux
  implicit none

contains

  !=============================================================================
  subroutine get_b2(ixI^L,ixO^L,w,x,b2,w_is_primitive)
    ! Obtains b**2, the invariant co-moving magnetic field.

    use mod_metric, only: lower3
    include 'amrvacdef.f' 

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,nw), intent(in)      :: w
    double precision, dimension(ixI^S,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixI^S), intent(out)        :: b2
    ! .. local ..
    logical                                                :: is_primitive
    double precision                                       :: bD(ixI^S,1:^NC)
    double precision                                       :: tmp(ixI^S), BdotV(ixI^S)
    !-----------------------------------------------------------------------------
    {#IFDEF DY_SP
      call mpistop('define DY_SP should not use get_b2, you should use imhd_get_intermediate_variables')
    }

    if(.not.present(w_is_primitive)) then
       is_primitive = .false.  ! By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if
    is_primitive = .true.

    call lower3(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),bD)
    !if (is_primitive) then
    !   BdotV(ixO^S) = {^C& bD(ixO^S,^C)*w(ixO^S,u^C_)+}
    !   if (useprimitiveRel) BdotV(ixO^S) = BdotV(ixO^S)/w(ixO^S,lfac_)
    !else
    !   BdotV(ixO^S)= ({^C& w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_)
    !end if
       BdotV(ixO^S) = {^C& bD(ixO^S,^C)*w(ixO^S,u^C_)+}
       if (useprimitiveRel) BdotV(ixO^S) = BdotV(ixO^S)/w(ixO^S,lfac_)

    tmp(ixO^S) = {^C& w(ixO^S,b^C_)*bD(ixO^S,^C)+}

    b2(ixO^S) = tmp(ixO^S)/w(ixO^S,lfac_)**2 + BdotV(ixO^S)**2

  end subroutine get_b2

  subroutine get_b2_pt(ix^D,w_pt,x_pt,b2)
    ! Obtains b**2, the invariant co-moving magnetic field.

    use mod_metric
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ix^D
    double precision, dimension(1:nw), intent(in)          :: w_pt
    double precision, dimension(1:ndim), intent(in)        :: x_pt
    double precision, intent(out)                          :: b2
    !type(metric), pointer                                  :: m
    ! .. local ..
    double precision                                       :: bD(1:^NC)
    double precision                                       :: tmp, BdotV
    !-----------------------------------------------------------------------------
    {#IFDEF DY_SP
      call mpistop('define DY_SP should not use get_b2_pt, you should use imhd_get_intermediate_variables')
    }

    call lower3_ixD(ix^D,myM,w_pt(b1_:b^NC_),bD)
    BdotV = {^C& bD(^C)*w_pt(u^C_)+}
    if (useprimitiveRel) BdotV = BdotV/w_pt(lfac_)

    tmp = {^C& w_pt(b^C_)*bD(^C)+}

    b2 = tmp/w_pt(lfac_)**2 + BdotV**2

  end subroutine get_b2_pt
 
  !=============================================================================
  subroutine get_u4(ixI^L,ixO^L,w,x,u4,w_is_primitive)
    ! Obtains the contravariant four-velocity u4
    ! from conservative or primitive variables

    use mod_metric, only: raise3, square3u
    include 'amrvacdef.f' 

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,nw), intent(in)      :: w
    double precision, dimension(ixI^S,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixI^S,0:ndir), intent(out) :: u4
    ! .. local ..
    logical                                                :: is_primitive
    double precision                                       :: v(ixI^S)
    double precision, dimension(ixI^S)                     :: a, b, c, d
    double precision, dimension(ixI^S)                     :: VdotB, sqrB
    double precision, dimension(ixI^S,1:ndir)              :: bD, sU
    !-----------------------------------------------------------------------------
    {#IFDEF DY_SP
      call mpistop('define DY_SP should not use get_u4, you should use imhd_get_intermediate_variables')
    }

    if(.not.present(w_is_primitive)) then
       is_primitive = .false.  ! By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    is_primitive = .true.
    !if (is_primitive) then
    !   {^C& u4(ixO^S,^C) = w(ixO^S,u^C_) - w(ixO^S,lfac_)/myM%alpha(ixO^S)*myM%beta(^C)%elem(ixO^S) \}       
    !else
    !   VdotB(ixO^S)= ({^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_)
    !   call square3u(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),sqrB(ixI^S))
    !   call raise3(ixI^L,ixO^L,myM,w(ixI^S,s1_:s^NC_),sU)

    !   {^C&
    !      v(ixO^S)=(sU(ixO^S,^C)+VdotB(ixO^S)*w(ixO^S,b0_+^C))/ &
    !        (w(ixO^S,xi_)+sqrB(ixO^S))
    !   u4(ixO^S,^C) = w(ixO^S,lfac_) * &
    !        ( v(ixO^S) - myM%beta(^C)%elem(ixO^S)/myM%alpha(ixO^S) )
    !   \}
    !end if
       {^C& u4(ixO^S,^C) = w(ixO^S,u^C_) - w(ixO^S,lfac_)/myM%alpha(ixO^S)*myM%beta(^C)%elem(ixO^S) \}       

    ! Obtain the time-component from u_alpha * u**alpha = -1

    a(ixO^S) = myM%g(0,0)%elem(ixO^S)
    b(ixO^S) = 2.0d0*({^C& u4(ixO^S,^C)*myM%g(^C,0)%elem(ixO^S)|+})
    c(ixO^S) = 1.0d0 &
         + ({^C& u4(ixO^S,^C)**2*myM%g(^C,^C)%elem(ixO^S)|+}) &
         + 2.0d0 * ({^NOONEC myM%g(1,2)%elem(ixO^S) * u4(ixO^S,1) * u4(ixO^S,2) &
         {#IFDEF C3 + myM%g(1,3)%elem(ixO^S)*u4(ixO^S,1) * u4(ixO^S,3) &
         + myM%g(2,3)%elem(ixO^S)*u4(ixO^S,2) * u4(ixO^S,3)}})
    d(ixO^S) = b(ixO^S)**2 - 4.0d0*a(ixO^S)*c(ixO^S)

    u4(ixO^S,0) = -half/a(ixO^S)*(b(ixO^S)+sqrt(d(ixO^S)))
    
  end subroutine get_u4
  !=============================================================================
  subroutine get_b4(ixI^L,ixO^L,w,x,b4,w_is_primitive)
    ! Obtains the contravariant four-magnetic field in the fluid frame (little b) 
    ! from conservative or primitive variables
    !=============================================================================
    use mod_metric, only: lower3
    include 'amrvacdef.f' 

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,nw), intent(in)      :: w
    double precision, dimension(ixI^S,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixI^S,0:ndir), intent(out) :: b4
    ! .. local ..
    logical                                                :: is_primitive
    double precision, dimension(ixI^S,1:ndir)              :: vd, bD
    double precision, dimension(ixI^S,0:ndir)              :: u4
    double precision, dimension(ixI^S)                     :: VdotB, sqrB
    !-----------------------------------------------------------------------------
    {#IFDEF DY_SP
      call mpistop('define DY_SP should not use get_b4, you should use imhd_get_intermediate_variables')
    }

    if(.not.present(w_is_primitive)) then
       is_primitive = .false.  ! By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    is_primitive = .true.

    !if (is_primitive) then
    !   call lower3(ixI^L,ixO^L,myM,w(ixI^S,u1_:u^NC_),vd)
    !   {^C& vd(ixO^S,^C) = vd(ixO^S,^C)/w(ixO^S,lfac_) \}
    !else
    !   
    !   VdotB(ixO^S) = ({^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_)
    !   call lower3(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),bD)
    !   sqrB(ixO^S)  = {^C&w(ixO^S,b^C_)*bD(ixO^S,^C)+}

    !   {^C& 
    !   vd(ixO^S,^C) = (w(ixO^S,s^C_)+VdotB(ixO^S)*bd(ixO^S,^C))/ &
    !   (w(ixO^S,xi_)+sqrB(ixO^S))
    !   \}
    !   
    !end if ! is_primitive
       call lower3(ixI^L,ixO^L,myM,w(ixI^S,u1_:u^NC_),vd)
       {^C& vd(ixO^S,^C) = vd(ixO^S,^C)/w(ixO^S,lfac_) \}

    b4(ixO^S,0) = w(ixO^S,lfac_) * ({^C& w(ixO^S,b^C_)*vd(ixO^S,^C)|+})/myM%alpha(ixO^S)
    call get_u4(ixI^L,ixO^L,w,x,u4,w_is_primitive=is_primitive)
    {^C& b4(ixO^S,^C) = (w(ixO^S,b^C_)+myM%alpha(ixO^S)*b4(ixO^S,0)*u4(ixO^S,^C))/w(ixO^S,lfac_)\}

    
  end subroutine get_b4
  !=============================================================================
  subroutine get_TEM4_ud(ixI^L,ixO^L,w,x,tem4,w_is_primitive)
    ! Obtains the electromagnetic part of the four-stress-energy tensor 
    ! One index up, one index down.  
    !=============================================================================
    use mod_metric, only: lower4
    include 'amrvacdef.f' 

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,nw), intent(in)      :: w
    double precision, dimension(ixI^S,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixI^S,0:ndir,0:ndir), intent(out) :: tem4
    ! .. local ..
    logical                                                :: is_primitive
    double precision, dimension(ixI^S,0:ndir)              :: u4, b4, u4d, b4d
    double precision, dimension(ixI^S)                     :: b2
    integer                                                :: imu, inu
    !=============================================================================
    {#IFDEF DY_SP
      call mpistop('define DY_SP should not use get_TEM4_ud, you should use imhd_get_intermediate_variables')
    }

    if(.not.present(w_is_primitive)) then
       is_primitive = .false.  ! By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    is_primitive = .true.

    call get_u4(ixI^L,ixO^L,w,x,u4,w_is_primitive=is_primitive)
    call get_b4(ixI^L,ixO^L,w,x,b4,w_is_primitive=is_primitive)
    call get_b2(ixI^L,ixO^L,w,x,b2,w_is_primitive=is_primitive)

    call lower4(ixI^L,ixO^L,myM,u4,u4D)
    call lower4(ixI^L,ixO^L,myM,b4,b4D)

    do imu = 0, ndir
       do inu = 0, ndir
          tem4(ixO^S,imu,inu) = b2(ixO^S)*u4(ixO^S,imu)*u4d(ixO^S,inu) &
               + half*b2(ixO^S)*kr(imu,inu) &
               - b4(ixO^S,imu)*b4d(ixO^S,inu)
       end do
    end do

  end subroutine get_TEM4_ud
  !=============================================================================
  subroutine get_TPAKE4_ud(ixI^L,ixO^L,w,x,tpake4,w_is_primitive)
    ! Obtains the free particle energy part of the four-stress-energy tensor 
    ! as defined by Mc Kinney et al. (2012), MNRAS 423, p. 3083, Eq. (6)
    ! One index up, one index down.  
    !=============================================================================
    use mod_metric, only: lower4
    include 'amrvacdef.f' 

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,nw), intent(in)      :: w
    double precision, dimension(ixI^S,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixI^S,0:ndir,0:ndir), intent(out) :: tpake4
    ! .. local ..
    logical                                                :: is_primitive
    double precision, dimension(ixI^S,0:ndir)              :: u4, u4d
    double precision, dimension(ixI^S)                     :: rho
    integer                                                :: imu, inu
    !=============================================================================
    {#IFDEF DY_SP
      call mpistop('define DY_SP should not use get_TPAKE4_ud, you should use imhd_get_intermediate_variables')
    }

    if(.not.present(w_is_primitive)) then
       is_primitive = .false.  ! By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    is_primitive = .true.

    !if (is_primitive) then
    !   rho(ixO^S) = w(ixO^S,rho_)
    !else
    !   rho(ixO^S) = w(ixO^S,d_)/w(ixO^S,lfac_)
    !end if
    rho(ixO^S) = w(ixO^S,rho_)

    call get_u4(ixI^L,ixO^L,w,x,u4,w_is_primitive=is_primitive)
    call lower4(ixI^L,ixO^L,myM,u4,u4D)

    do imu = 0, ndir
       do inu = 0, ndir
          tpake4(ixO^S,imu,inu) = (u4d(ixO^S,inu)+one*kr(0,inu)) &
               * rho(ixO^S)*u4(ixO^S,imu)
       end do
    end do

  end subroutine get_TPAKE4_ud
  !=============================================================================
  subroutine get_TEN4_ud(ixI^L,ixO^L,w,x,ten4,w_is_primitive)
    ! Obtains the enthalpy part of the four-stress-energy tensor 
    ! as defined by Mc Kinney et al. (2012), MNRAS 423, p. 3083, Eq. (6)
    ! One index up, one index down.  
    !=============================================================================
    use mod_eos
    use mod_metric, only: lower4
    use mod_imhd_con2prim

    include 'amrvacdef.f' 

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,nw), intent(in)      :: w
    double precision, dimension(ixI^S,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixI^S,0:ndir,0:ndir), intent(out) :: ten4
    ! .. local ..
    logical                                                :: is_primitive
    double precision                                       :: wprim(ixI^S,1:nw), rho_local, temp_local
    double precision, dimension(ixI^S,0:ndir)              :: u4, u4d
    double precision, dimension(ixI^S)                     :: p, rhoh, e, eps_tmp
    integer                                                :: imu, inu, ix^D
    !=============================================================================
    {#IFDEF DY_SP
      call mpistop('define DY_SP should not use get_TEN4_ud, you should use imhd_get_intermediate_variables')
    }

    if(.not.present(w_is_primitive)) then
       is_primitive = .false.  ! By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    is_primitive = .true.

    {do ix^D = ixO^LIM^D \}
        if (eos_type == tabulated) then
          temp_local = w(ix^D,T_eps_)
          rho_local  = w(ix^D,rho_)
          call eos_temp_get_all_one_grid(rho_local,temp_local,w(ix^D,ye_),&
                                         eps_tmp(ix^D))
        else
          eps_tmp(ix^D) = w(ix^D, T_eps_)
        endif
    {enddo^D&\}

    p(ixO^S) = w(ixO^S, pp_)
    rhoh(ixO^S) = w(ixO^S, rho_) * (1.0d0 + eps_tmp(ixO^S) + w(ixO^S, pp_)/w(ixO^S, rho_))

    e(ixO^S) = rhoh(ixO^S) - w(ixO^S,rho_)

    call get_u4(ixI^L,ixO^L,w,x,u4,w_is_primitive=is_primitive)
    call lower4(ixI^L,ixO^L,myM,u4,u4D)

    do imu = 0, ndir
       do inu = 0, ndir
          ten4(ixO^S,imu,inu) = e(ixO^S)*u4(ixO^S,imu)*u4d(ixO^S,inu) &
               + kr(imu,inu) * p(ixO^S)
       end do
    end do

  end subroutine get_TEN4_ud
  !=============================================================================
end module mod_physaux
!=============================================================================
