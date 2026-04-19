! ====================================================================== 
! ====================================================================== 
!                 Module containing closure-related 
!                 subroutine for M1 radiation transport
! ======================================================================
! ======================================================================  
module mod_m1_closure
  ! ====================================================================== 
  {#IFDEF UNIT_TESTS 
  use mod_m1_tests
  } !end IFDEF UNIT_TESTS
  implicit none
  {#IFNDEF UNIT_TESTS
  integer, parameter  :: m1_npress = ((^NC)**2-(^NC))/2 + (^NC)
  integer, parameter  :: m1_nvars = 1 + 1 + (^NC)
    ! Parameters :
  ! IMPORTANT: this parameter indicates 
  ! the number of m1 specific variables per species 
  ! stored in the cons portion of w.
  ! N + E + F_i 
  }
  ! ====================================================================== 
  ! Just avoiding division by zero 
  double precision, parameter :: M1_TINY=1.0d-45
  ! ====================================================================== 

  ! ====================================================================== 
  ! Collection of helper variables used 
  ! inside the m1 closure for various purposes
  ! NB: if not specified by a "_low" 
  !     all indices are upper.
  ! ====================================================================== 
  type m1_closure_helpers
    ! Number of independent components of press tensor 
    ! integer, parameter :: m1_npress = ((^NC)**2-(^NC))/2 + (^NC)
    ! pieces of H_mu H^mu
    ! Pieces of J 
    double precision :: B0, Bthin, BThick
    ! Norm of F and scalar products 
    double precision :: E, F, F2, Fv, fhatv, Fden
    ! Helpful vectors 
    double precision, dimension(^NC) :: F_low, F_up, fhat_low, fhat_up
    double precision, dimension(^NC) :: H_low, H_up
    ! Fluid variables 
    double precision :: v2, W
    double precision, dimension(^NC) :: vel, vlow
    ! Eddington factor and helpers
    double precision :: zeta, chi, dthick, dthin
    ! Fluid frame moments 
    double precision :: J, Gamma, H0 !Hi
    ! root finding helpers
    double precision :: E_closure,HH0, HHthickthick, HHthinthin, HHthin, HHthick, HHthickthin
    ! Jacobian coefficients 
    double precision :: an, av, aF, an_thick, av_thick, aF_thick, an_thin, av_thin, af_thin
  end type m1_closure_helpers
  ! ====================================================================== 
  ! ====================================================================== 
  public
  ! ====================================================================== 
  ! ====================================================================== 
  procedure(sub_m1_closure_func), pointer :: m1_closure_func          => null()
  ! ====================================================================== 
  procedure(sub_m1_closure_deriv), pointer :: m1_closure_deriv        => null()
  ! ====================================================================== 
  ! ====================================================================== 
  abstract interface
  ! ======================================================================    

     ! ====================================================================== 
     function sub_m1_closure_func(zeta)
       double precision :: sub_m1_closure_func
       double precision, intent(in) :: zeta
     end function sub_m1_closure_func
     ! ====================================================================== 
     function sub_m1_closure_deriv(zeta)
       double precision :: sub_m1_closure_deriv
       double precision, intent(in) :: zeta
     end function sub_m1_closure_deriv
     ! ====================================================================== 
  end interface
! ====================================================================== 
! ====================================================================== 
contains
  ! ====================================================================== 
  subroutine m1_closure_type_activate()
    !use mod_m1_minerbo_closure
    !use mod_m1_levermore_closure
    use mod_m1_closure_functions
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"  
    }
    select case(m1_closure_type)
      case("Minerbo")
        m1_closure_func    => minerbo_closure
        m1_closure_deriv   => minerbo_deriv
      case("Levermore")
        m1_closure_func    => levermore_closure
        m1_closure_deriv   => levermore_deriv
      !>by default minerbo
    end select
  end subroutine m1_closure_type_activate


  ! ====================================================================== 
  subroutine m1_closure_activate()
    call m1_closure_type_activate
    if( .not. associated(m1_closure_func) ) then
      {#IFNDEF UNIT_TESTS
       call mpistop("You need to provide a closure for M1 to work")
      }{#IFDEF UNIT_TESTS
       stop "You need to provide a closure for M1 to work"
      }       
    end if
  end subroutine m1_closure_activate
  ! ====================================================================== 
  !< Update the closure on the whole grid. Rootfinding to find 
  !< zeta is optional, Gamma and Press are optional arguments 
  !< that, when passed, are filled with the N normalisation 
  !< and the pressure tensor, respectively. This routine 
  !< is for a single species, since it is meant to be called 
  !< from top level m1 routines contained in mod_m1.
  subroutine m1_update_closure(metricM1,wprim,x,ixI^L,ixO^L,species,get_zeta,Gamma,Hup,Jrad,Press,chi,W,vel,zetaOut,FilePrint)
    use mod_m1_metric_interface
    use, intrinsic :: ieee_arithmetic 
    {#IFNDEF UNIT_TESTS
    use mod_variables
    include "amrvacdef.f"  
    }
    ! ====================================================================== 
    ! ====================================================================== 
    integer, intent(in) :: ixI^L,ixO^L,species
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: wprim(ixI^S,1:nw)
    ! Optional return arguments 
    double precision, dimension(ixI^S), intent(out), optional :: Gamma 
    double precision, dimension(ixI^S,m1_npress), intent(out), optional &
        :: Press
    double precision, dimension(ixI^S,^NC), intent(out), optional :: Hup 
    double precision, dimension(ixI^S,^NC), intent(out), optional :: vel  
    double precision, dimension(ixI^S), intent(out), optional     :: Jrad  
    double precision, dimension(ixI^S), intent(out), optional     :: chi  
    double precision, dimension(ixI^S), intent(out), optional     :: W 
    double precision, dimension(ixI^S), intent(out), optional     :: zetaOut
    logical, intent(in), optional :: get_zeta
    integer, optional :: FilePrint 
    ! ====================================================================== 
    ! Metric and closure helpers 
    type(m1_metric_helper), intent(in)   :: metricM1 ! KEN changed intent to in
    type(m1_closure_helpers) :: stateM1 
    ! Grid loop index 
    integer :: ix^D 
    integer :: tensidx, i, j
    logical :: get_zeta_impl 
    double precision                   :: JT, W2, W3, lfac_tmp 
	  double precision, dimension(^NC)   :: HT_up
    integer :: nidx, eidx, fidx^C, zidx 
    logical :: velocity_set = .false. 
    double precision :: vel_value = 1.0d-10
    ! ====================================================================== 
    ! ====================================================================== 
    if( present(get_zeta ) ) then 
      get_zeta_impl = get_zeta 
    else 
      get_zeta_impl = .true. 
    endif
    
    {#IFNDEF UNIT_TESTS
    ! Resolve indices 
    if( species .gt. ^NS ) then
      call mpistop("In closure, got species > NS !")
    endif 
    }
    ! Note: nspecies in [1,3], for first species ^KSP=1
    ! nidx = nrad1_ +(spec-1)*5 = 1+ (1-1)*5 = 1 = nrad1_ , for species = 1 
    ! nidx = nrad1_ +(spec-1)*5 = 1+ (2-1)*5 = 6 = nrad1_ +5 , for species = 2 
    nidx = nrad1_ + (species-1) * m1_nvars 
    eidx = erad1_ + (species-1) * m1_nvars 
    {^C& fidx^C = frad1^C_ + (species-1) * m1_nvars \} 

    {#IFDEF UNIT_TESTS2
      call fill_metric(metricM1)
    }

    ! Grid loop
    {^D& do ix^D=ixOmin^D,ixOmax^D \}

      stateM1%E = wprim(ix^D,eidx)
      {^C& stateM1%F_low(^C) = wprim(ix^D, fidx^C) \} 
      {^C& stateM1%vel(^C)   = wprim(ix^D,u0_+^C) \}  ! this is W v 
      if(velocity_set) then
        {^C&  stateM1%vel(^C) = vel_value \}
      end if 
      stateM1%zeta = 0.0d0 

      call m1_update_closure_ixD(stateM1,metricM1,ix^D,get_zeta_impl)

      if( present(zetaOut) ) then 
         zetaOut(ix^D)= stateM1%zeta
      end if

      ! Fill Gamma normalisation if present in input 
      if( present(Gamma) ) then 
        Gamma(ix^D) = stateM1%Gamma 
      endif 
      ! Fill pressure tensor if present in input 
      if( present(Press) ) then 
      W2 = stateM1%W**2
      W3 = stateM1%W**3 
      ! J in the optically thick limit 
      JT = 3.0d0/(2.0d0*W2 + 1.0d0)*((2.0d0*W2-1.0d0)*stateM1%E-2.0d0*W2*stateM1%Fv)
      ! H_i in the optically thick limit
      {^C& HT_up(^C) = stateM1%F_up(^C)/(stateM1%W+M1_TINY) + stateM1%vel(^C)*stateM1%W/(2.0d0*W2+1.0d0)&
                          *((4.0d0*W2+1.0d0)*stateM1%Fv-4.0d0*W2*stateM1%E) \}
      tensidx= 1
      do i=1, ^NC
        do j=i, ^NC
            Press(ix^D,tensidx) = stateM1%dthick * ( JT/3.0d0 * (4.0d0*W2 * stateM1%vel(i)*stateM1%vel(j) &
                + metricM1%gammaUPij(ix^D,tensidx) ) + stateM1%W*(HT_up(i)*stateM1%vel(j) + HT_up(j)*stateM1%vel(i) ) )&
                + stateM1%dthin * stateM1%E* stateM1%F_up(i)*stateM1%F_up(j)/(stateM1%F2+M1_TINY)

            tensidx = tensidx + 1 
        
        enddo
      enddo
      endif ! if( present(Press) ) 

      if( present(Hup) ) then 
        ! need extra array as stateM1 has no ix^D dimension
        call metricM1%raise_ixD(ix^D,stateM1%H_low, stateM1%H_up) 
        {^C& Hup(ix^D,^C) = stateM1%H_up(^C)  \} 
      endif 

      if( present(Jrad) ) then 
        Jrad(ix^D) = stateM1%J 
      endif

      if( present(chi) ) then 
        chi(ix^D) = stateM1%chi 
      endif  

      if( present(W) ) then 
        W(ix^D) = stateM1%W 
      endif 

      if( present(vel) ) then 
        vel(ix^D,:) = stateM1%vel 
        if(velocity_set) then
          vel(ix^D,:) = vel_value
        end if 
      endif 

      ! Also overwrite E and F_i since 
      ! limiting happens within closure 
      ! update procedure 
      wprim(ix^D,eidx) = stateM1%E
      {^C& wprim(ix^D, fidx^C) = stateM1%F_low(^C) \} 

    {^D& end do \}
    ! end of grid loop
    
  end subroutine m1_update_closure
!!  } ! end IFNDEF unit_tests
  ! ====================================================================== 
  
  ! ====================================================================== 
  ! Main (pointwise) closure update routine. This is pointwise 
  ! since it requires rootfinding. It has the option 
  ! of not updating zeta (in which case it simply returns
  ! the auxiliary variables).
  subroutine m1_update_closure_ixD(stateM1,metricM1,ix^D,get_zeta,get_vel_impl,use_init_zeta)
    ! ====================================================================== 
    ! ====================================================================== 
    use mod_m1_metric_interface
    use mod_m1_closure_functions
    use mod_rootfinding, only: rootfinding_constrained_newton_raphson, rootfinding_brent, rootfinding_constrained_newton_raphson_M1
    use, intrinsic :: ieee_arithmetic
    ! ====================================================================== 
    {#IFNDEF UNIT_TESTS 
    include "amrvacdef.f"
    }
    ! ====================================================================== 
    ! ====================================================================== 
    class(m1_closure_helpers), intent(inout) :: stateM1
    type(m1_metric_helper), intent(in)   :: metricM1 
    integer, intent(in) :: ix^D 
    logical, intent(in) :: get_zeta  
    logical, intent(in), optional :: get_vel_impl
    logical, intent(in), optional :: use_init_zeta
    ! ====================================================================== 
    !rootfinding helpers
    integer :: flag
    double precision :: fluxfact,W2, W3 
    double precision :: E1 ! internal E
    logical :: get_vel_W
    double precision :: Fmag, fac, Fmax2
    double precision :: fmaxfact = 0.999d0 !0.99999d9 !1.0d0 - 1.0d-10
    double precision :: epsF_rel = 1.0d-6
    double precision :: epsF_abs = 1.0d-30
    double precision :: Hmag2, Jsafe, scaleH, Jfloor
    double precision :: Hn_target,Hn_curr,vdotvlow,corr
    double precision, parameter :: epsH = 1.0d-12
    double precision, parameter :: tau_F = 1.0d-2
    double precision :: zm, zp 
    logical :: smallF, closure_bad
    logical :: use_init_zeta_impl
    double precision :: gUU_XX, gUU_XY, gUU_XZ, gUU_YY, gUU_YZ, gUU_ZZ
    double precision :: gDD_XX, gDD_XY, gDD_XZ, gDD_YY, gDD_YZ, gDD_ZZ


    if( present(get_vel_impl ) ) then 
      get_vel_W = get_vel_impl
    else 
      get_vel_W = .true. 
    endif

    if( present(use_init_zeta ) ) then 
      use_init_zeta_impl = use_init_zeta
    else 
      use_init_zeta_impl = .false. 
    endif

    gDD_XX = metricM1%gammaij(ix^D, 1)
    gDD_XY = metricM1%gammaij(ix^D, 2)
    gDD_XZ = metricM1%gammaij(ix^D, 3)
    gDD_YY = metricM1%gammaij(ix^D, 4)
    gDD_YZ = metricM1%gammaij(ix^D, 5)
    gDD_ZZ = metricM1%gammaij(ix^D, 6)

    gUU_XX = metricM1%gammaUPij(ix^D, 1)
    gUU_XY = metricM1%gammaUPij(ix^D, 2)
    gUU_XZ = metricM1%gammaUPij(ix^D, 3)
    gUU_YY = metricM1%gammaUPij(ix^D, 4)
    gUU_YZ = metricM1%gammaUPij(ix^D, 5)
    gUU_ZZ = metricM1%gammaUPij(ix^D, 6)

    ! sanitize incoming state
    if (.not. ieee_is_finite(stateM1%E)) stateM1%E = m1_E_atmo
    {^C& if (.not. ieee_is_finite(stateM1%F_low(^C))) stateM1%F_low(^C) = 0.0d0 \}
    {^C& if (.not. ieee_is_finite(stateM1%vel(^C)))   stateM1%vel(^C)   = 0.0d0 \}

    ! lower vel -> vlow
    stateM1%vlow(1) = gDD_XX*stateM1%vel(1) + gDD_XY*stateM1%vel(2) + gDD_XZ*stateM1%vel(3)
    stateM1%vlow(2) = gDD_XY*stateM1%vel(1) + gDD_YY*stateM1%vel(2) + gDD_YZ*stateM1%vel(3)
    stateM1%vlow(3) = gDD_XZ*stateM1%vel(1) + gDD_YZ*stateM1%vel(2) + gDD_ZZ*stateM1%vel(3)

    ! raise F_low -> F_up
    stateM1%F_up(1) = gUU_XX*stateM1%F_low(1) + gUU_XY*stateM1%F_low(2) + gUU_XZ*stateM1%F_low(3)
    stateM1%F_up(2) = gUU_XY*stateM1%F_low(1) + gUU_YY*stateM1%F_low(2) + gUU_YZ*stateM1%F_low(3)
    stateM1%F_up(3) = gUU_XZ*stateM1%F_low(1) + gUU_YZ*stateM1%F_low(2) + gUU_ZZ*stateM1%F_low(3)

    {#IFDEF METRIC_ZERO_COMP
    if(^ND .lt. 3) then
      stateM1%F_low(3) = 0.0d0
      stateM1%F_up(3) = 0.0d0
      if(^ND .lt. 2) then
        stateM1%F_low(2) = 0.0d0
        stateM1%F_up(2) = 0.0d0
      end if 
    end if 
    }

    ! This is a point by point loop since rootfinding is required 
    E1 = max(stateM1%E, m1_E_atmo)
    if (.not. ieee_is_finite(E1)) E1 = m1_E_atmo

    stateM1%F2 = {^C& stateM1%F_low(^C)*stateM1%F_up(^C) +}
    if (.not. ieee_is_finite(stateM1%F2)) then
      stateM1%F2 = 0.0d0
      {^C& stateM1%F_low(^C) = 0.0d0 \}
      {^C& stateM1%F_up(^C)  = 0.0d0 \}
    end if
    if (stateM1%F2 < 0.0d0) then
      if (dabs(stateM1%F2) < 1.0d-12*max(E1*E1, M1_TINY)) then
        stateM1%F2 = 0.0d0
      else
        stateM1%F2 = 0.0d0
        {^C& stateM1%F_low(^C) = 0.0d0 \}
        {^C& stateM1%F_up(^C)  = 0.0d0 \}
      end if
    end if
    !stateM1%F2 = max(stateM1%F2, M1_TINY) !0.0d0
    stateM1%F2 = max(stateM1%F2, 0.0d0)

    if(get_vel_W) then
       ! What we get as input inside of vel is Wv^i
       ! we need to convert it to the eulerian 3-velocity
       W2 = 1.0d0 + {^C&stateM1%vlow(^C)*stateM1%vel(^C)+}
       if (.not. ieee_is_finite(W2)) W2 = 1.0d0
       W2 = max(W2, 1.0d0)
       stateM1%W = dsqrt(W2)
       W3 = stateM1%W**3 
       {^C& 
       stateM1%vlow(^C) = stateM1%vlow(^C) / max(stateM1%W, 1.0d0)
       stateM1%vel(^C)  = stateM1%vel(^C)  / max(stateM1%W, 1.0d0)
       \} 
       stateM1%v2 = 1.0d0 - 1.0d0/max(stateM1%W**2, 1.0d0)
    else
       stateM1%v2 = {^C&stateM1%vlow(^C)*stateM1%vel(^C)+}
       if (.not. ieee_is_finite(stateM1%v2)) stateM1%v2 = 0.0d0
       stateM1%v2 = max(0.0d0, min(stateM1%v2, 1.0d0 - 1.0d-12))
       W2 = 1.0d0/max(1.0d0 - stateM1%v2, 1.0d-12)
       stateM1%W = dsqrt(W2)
       W3 = stateM1%W**3 
    end if 

    if (.not. ieee_is_finite(stateM1%W)) stateM1%W = 1.0d0
    if (.not. ieee_is_finite(stateM1%v2)) stateM1%v2 = 0.0d0
    stateM1%v2 = max(0.0d0, min(stateM1%v2, 1.0d0 - 1.0d-12))
    W2 = stateM1%W*stateM1%W

    ! enforce strict F^2 < E^2 admissibility
    Fmax2 = (fmaxfact*E1)*(fmaxfact*E1)
    if (stateM1%F2 > Fmax2) then
      fluxfact = dsqrt(Fmax2 / max(stateM1%F2, M1_TINY))
      {^C& stateM1%F_low(^C) = fluxfact * stateM1%F_low(^C) \}
      {^C& stateM1%F_up(^C)  = fluxfact * stateM1%F_up(^C)  \}
      stateM1%F2 = Fmax2
    end if

   ! Fmag = dsqrt(max(stateM1%F2, M1_TINY))
    Fmag = dsqrt(max(stateM1%F2, 0.0d0))
    stateM1%F = Fmag

    !---------- check if flux is small --------------
    !smallF = (Fmag <= tau_F*E1)
    smallF = (Fmag <= max(tau_F*E1, epsF_abs))
    if(smallF) then 
      ! treat radiation as isotropic and remove all hidden 1/F paths
      stateM1%F = 0.0d0
      stateM1%F2 = 0.0d0
      stateM1%Fv = 0.0d0
      stateM1%fhatv = 0.0d0
      stateM1%Fden = 1.0d0
      {^C& stateM1%F_low(^C) = 0.0d0 \}
      {^C& stateM1%F_up(^C)  = 0.0d0 \}
      {^C& stateM1%fhat_low(^C) = 0.0d0 \}
      {^C& stateM1%fhat_up(^C)  = 0.0d0 \}
    else  
      stateM1%Fden = max(Fmag, max(epsF_abs, epsF_rel*E1))
      {^C& stateM1%fhat_low(^C) = stateM1%F_low(^C) / stateM1%Fden \}
      {^C& stateM1%fhat_up(^C)  = stateM1%F_up(^C)  / stateM1%Fden \}

      stateM1%Fv    = (^C&stateM1%F_low(^C)*stateM1%vel(^C)+)    
      stateM1%fhatv = (^C&stateM1%fhat_low(^C)*stateM1%vel(^C)+)

      if (.not. ieee_is_finite(stateM1%Fv))    stateM1%Fv = 0.0d0
      if (.not. ieee_is_finite(stateM1%fhatv)) stateM1%fhatv = 0.0d0

      if (dabs(stateM1%fhatv) < 1.0d-15) stateM1%fhatv = 0.0d0
      stateM1%fhatv = max(-1.d0 + 1.d-12, min(1.d0 - 1.d-12, stateM1%fhatv))
    end if 

    stateM1%B0     = W2*( E1 - 2.0d0*stateM1%Fv)
    stateM1%Bthin  = W2*E1*stateM1%fhatv**2
    stateM1%BThick = (W2-1.0d0)/(1.0d0+2.0d0*W2)*(4.0d0*W2*stateM1%Fv+(3.0d0-2.0d0*W2)*E1)
    stateM1%an = stateM1%W*stateM1%B0 + stateM1%W*(stateM1%Fv-E1)
    stateM1%av = stateM1%W*stateM1%B0
    stateM1%aF = -stateM1%W

    stateM1%an_thick = stateM1%W*stateM1%BThick
    stateM1%av_thick = stateM1%W*stateM1%BThick + stateM1%W/(2.0d0*W2+1.0d0)*((3.0d0-2.0d0*W2)*E1&
        +(2.0d0*W2-1.0d0)*stateM1%Fv)

    stateM1%aF_thick = stateM1%W*stateM1%v2

    stateM1%an_thin = stateM1%W*stateM1%Bthin
    stateM1%av_thin = stateM1%an_thin
    stateM1%af_thin = stateM1%W*E1*stateM1%fhatv

    closure_bad = .false.
    if (.not. ieee_is_finite(stateM1%B0))     closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%Bthin))  closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%BThick)) closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%an))     closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%av))     closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%aF))     closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%an_thick)) closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%av_thick)) closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%aF_thick)) closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%an_thin))  closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%av_thin))  closure_bad = .true.
    if (.not. ieee_is_finite(stateM1%af_thin))  closure_bad = .true.

    if( get_zeta ) then 
      if (closure_bad) then
        stateM1%zeta = 0.0d0
      else if ( stateM1%v2 < 1.0d-15 ) then
        stateM1%zeta = min(dsqrt(max(stateM1%F2,0.0d0)/max(E1*E1, M1_TINY)),1.0d0)
      else
        ! Need rootfinding to determine zeta
        ! first we compute the pieces of H_sq
        stateM1%E_closure = max(E1, m1_E_atmo) + M1_TINY

        stateM1%HH0 = stateM1%av**2 *stateM1%v2 + stateM1%aF**2*stateM1%F2 + 2.0d0*stateM1%av*stateM1%aF*stateM1%Fv - stateM1%an**2

        stateM1%HHthickthick = stateM1%av_thick**2*stateM1%v2 + stateM1%aF_thick**2*stateM1%F2 &
              +2.0d0*stateM1%aF_thick*stateM1%av_thick*stateM1%Fv - stateM1%an_thick**2
        stateM1%HHthinthin = stateM1%av_thin**2*stateM1%v2 + stateM1%af_thin**2 + 2.0d0*stateM1%af_thin*stateM1%av_thin*stateM1%fhatv &
              - stateM1%an_thin**2

        stateM1%HHthin = 2.0d0*(stateM1%av*stateM1%av_thin*stateM1%v2 + stateM1%aF*stateM1%af_thin * stateM1%F &
              + stateM1%af_thin*stateM1%av*stateM1%Fv/stateM1%Fden + stateM1%av_thin*stateM1%aF*stateM1%Fv - stateM1%an_thin*stateM1%an)
        stateM1%HHthick = 2.0d0*(stateM1%av*stateM1%av_thick*stateM1%v2 + stateM1%aF*stateM1%aF_thick*stateM1%F2 &
              + stateM1%aF_thick*stateM1%av*stateM1%Fv + stateM1%av_thick*stateM1%aF*stateM1%Fv - stateM1%an_thick*stateM1%an)

        stateM1%HHthickthin = 2.0d0*(stateM1%av_thin*stateM1%av_thick*stateM1%v2 + stateM1%af_thin*stateM1%aF_thick * stateM1%F &
              + stateM1%af_thin*stateM1%av_thick*stateM1%fhatv + stateM1%av_thin*stateM1%aF_thick*stateM1%Fv - stateM1%an_thin*stateM1%an_thick)

        if (.not. ieee_is_finite(stateM1%HH0))          closure_bad = .true.
        if (.not. ieee_is_finite(stateM1%HHthickthick)) closure_bad = .true.
        if (.not. ieee_is_finite(stateM1%HHthinthin))   closure_bad = .true.
        if (.not. ieee_is_finite(stateM1%HHthin))       closure_bad = .true.
        if (.not. ieee_is_finite(stateM1%HHthick))      closure_bad = .true.
        if (.not. ieee_is_finite(stateM1%HHthickthin))  closure_bad = .true.

        if (closure_bad) then
          stateM1%zeta = 0.0d0
        else
          zm = 0.0d0
          zp = 1.0d0
          if(use_init_zeta_impl) then
            zm = max(0.0d0, stateM1%zeta - 0.01d0)
            zp = min(1.0d0, stateM1%zeta + 0.01d0)
          endif
          call rootfinding_constrained_newton_raphson_M1(stateM1%zeta, zm, zp, 1.0d-15, 100, flag, xi, dxidz)

          if (flag .ne. 0 .or. .not. ieee_is_finite(stateM1%zeta) .or. stateM1%zeta > 1.0d0 .or. stateM1%zeta < 0.0d0) then
            stateM1%zeta = 0.5d0
            call rootfinding_brent(stateM1%zeta, 0.0d0, 1.0d0, 1.0d-15, 100, flag, xi)
            if (flag .ne. 0 .or. .not. ieee_is_finite(stateM1%zeta)) then
              stateM1%zeta = 0.0d0
            end if
          end if
        end if
      end if
    end if 

    ! Finally fill in primitive array:
    ! Compute J, H_i and store zeta
    if (.not. ieee_is_finite(stateM1%zeta)) stateM1%zeta = 0.0d0
    stateM1%zeta = max(0.0d0, min(stateM1%zeta, 1.0d0))

    stateM1%chi = m1_closure_func(stateM1%zeta)
    if (.not. ieee_is_finite(stateM1%chi)) stateM1%chi = 1.0d0/3.0d0
    stateM1%chi = max(1.0d0/3.0d0, min(stateM1%chi, 1.0d0))

    stateM1%dthin = 1.5d0*stateM1%chi-0.5d0
    stateM1%dthick = 1.0d0 - stateM1%dthin

    stateM1%J = stateM1%B0 + stateM1%dthin*stateM1%Bthin + stateM1%dthick*stateM1%BThick
    Jfloor = max(M1_TINY, 1.0d-12*E1)
    if (.not. ieee_is_finite(stateM1%J)) stateM1%J = Jfloor

    {^C& stateM1%H_low(^C) = - (stateM1%av+stateM1%dthin*stateM1%av_thin+stateM1%dthick*stateM1%av_thick)*stateM1%vlow(^C) &
        -(stateM1%aF+stateM1%dthick*stateM1%aF_thick)*stateM1%F_low(^C) &
        -stateM1%af_thin*stateM1%dthin*stateM1%fhat_low(^C) \}

    ! raise H_low -> H_up
    stateM1%H_up(1) = gUU_XX*stateM1%H_low(1) + gUU_XY*stateM1%H_low(2) + gUU_XZ*stateM1%H_low(3)
    stateM1%H_up(2) = gUU_XY*stateM1%H_low(1) + gUU_YY*stateM1%H_low(2) + gUU_YZ*stateM1%H_low(3)
    stateM1%H_up(3) = gUU_XZ*stateM1%H_low(1) + gUU_YZ*stateM1%H_low(2) + gUU_ZZ*stateM1%H_low(3)

    if (stateM1%J < Jfloor) then
      stateM1%J = Jfloor
      stateM1%zeta = 0.0d0
      stateM1%chi = 1.0d0/3.0d0
      stateM1%dthin = 0.0d0
      stateM1%dthick = 1.0d0
      {^C& stateM1%H_low(^C) = 0.0d0 \}
      {^C& stateM1%H_up(^C)  = 0.0d0 \}
    end if

    Hmag2 = {^C& stateM1%H_low(^C)*stateM1%H_up(^C) +}
    if (.not. ieee_is_finite(Hmag2)) then
      Hmag2 = 0.0d0
      {^C& stateM1%H_low(^C) = 0.0d0 \}
      {^C& stateM1%H_up(^C)  = 0.0d0 \}
    end if
    Hmag2 = max(Hmag2, 0.0d0)
    Jsafe = max(stateM1%J, M1_TINY)

    if (Hmag2 >= (1.d0 - epsH)**2 * Jsafe*Jsafe) then
      scaleH = (1.d0 - epsH) * Jsafe / dsqrt(max(Hmag2, M1_TINY))
      {^C& stateM1%H_low(^C) = scaleH*stateM1%H_low(^C) \}
      {^C& stateM1%H_up(^C)  = scaleH*stateM1%H_up(^C)  \}
    end if

    stateM1%Gamma = stateM1%W*(E1-stateM1%Fv)/max(stateM1%J, M1_TINY)
    if (.not. ieee_is_finite(stateM1%Gamma)) stateM1%Gamma = 1.0d0
    stateM1%Gamma = max(stateM1%Gamma, M1_TINY)
    stateM1%Gamma = min(stateM1%Gamma, 1.0d12)

    stateM1%E = E1

    ! ====================================================================== 
    !                END OF CLOSURE CALCULATION 
    ! ====================================================================== 

    contains 

    function xi(z)
      implicit none 
      double precision :: xi
      double precision, intent(in) :: z
      double precision :: J_closure, H_sq
      double precision :: chi, dthick, dthin 
      chi = m1_closure_func(z)
      if (.not. ieee_is_finite(chi)) chi = 1.0d0/3.0d0
      chi = max(1.0d0/3.0d0, min(chi, 1.0d0))
      dthin = 1.5d0 * chi - 0.5d0
      dthick = 1.0d0 - dthin
      J_closure = stateM1%B0 + dthin*stateM1%Bthin + dthick*stateM1%BThick
      H_sq      = stateM1%HH0 + dthin*dthin*stateM1%HHthinthin &
	      + dthick*dthick*stateM1%HHthickthick &
          + dthin*dthick * stateM1%HHthickthin + dthick * stateM1%HHthick + dthin  * stateM1%HHthin

      if ((.not. ieee_is_finite(J_closure)) .or. (.not. ieee_is_finite(H_sq))) then
        xi = 1.0d99
      else
        xi = ( J_closure**2*z**2 - H_sq ) / max(stateM1%E_closure**2, M1_TINY)
      end if
    end function 

    function dxidz(z)
      double precision :: dxidz
      double precision, intent(in) :: z
      double precision :: J_closure, H_sq
      double precision :: chi, dthin, dthick
      double precision :: dchidz, d_dthin_dz, d_dthick_dz
      double precision :: dJdz, dHsqddthick, dHsqddthin

      chi = m1_closure_func(z)
      if (.not. ieee_is_finite(chi)) chi = 1.0d0/3.0d0
      chi = max(1.0d0/3.0d0, min(chi, 1.0d0))
      dthin = 1.5d0 * chi - 0.5d0
      dthick = 1.0d0 - dthin

      J_closure = stateM1%B0 + dthin*stateM1%Bthin + dthick*stateM1%BThick
      H_sq      = stateM1%HH0 + dthin*dthin*stateM1%HHthinthin &
	      + dthick*dthick*stateM1%HHthickthick &
          + dthin*dthick * stateM1%HHthickthin + dthick * stateM1%HHthick + dthin  * stateM1%HHthin

      dchidz = m1_closure_deriv(z)
      if (.not. ieee_is_finite(dchidz)) dchidz = 0.0d0
      d_dthin_dz = 1.5d0*dchidz
      d_dthick_dz = -1.5d0*dchidz

      dJdz = stateM1%Bthin*d_dthin_dz + stateM1%BThick * d_dthick_dz

      dHsqddthick = stateM1%HHthick + stateM1%HHthickthin*dthin + 2.0d0 * stateM1%HHthickthick * dthick
      dHsqddthin  = stateM1%HHthin + stateM1%HHthickthin*dthick + 2.0d0 * stateM1%HHthinthin * dthin

      dxidz = 2.0d0*J_closure*dJdz*z**2 + 2.0d0*J_closure**2*z &
           - d_dthin_dz*dHsqddthin - d_dthick_dz * dHsqddthick

      if ((.not. ieee_is_finite(dxidz)) .or. (.not. ieee_is_finite(J_closure)) .or. (.not. ieee_is_finite(H_sq))) then
        dxidz = 0.0d0
      else
        dxidz = dxidz / max(stateM1%E_closure * stateM1%E_closure, M1_TINY)
      end if
    end function dxidz
    ! ====================================================================== 
    ! ====================================================================== 
  end subroutine m1_update_closure_ixD
  ! ====================================================================== 

! ====================================================================== 
! ====================================================================== 
! ====================================================================== 
end module mod_m1_closure
! ====================================================================== 
!                            END OF FILE 
! ====================================================================== 

