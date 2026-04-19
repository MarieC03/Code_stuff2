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
  ! Metric components for Cylindrified modified Kerr-Schild coordinates -coord=cmks
  !
  !
  ! Oliver Porth
  ! 17.01.2017 : started with module
  ! 18.03.2017 : Fully implemented
  !=============================================================================
  
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------
  
  character*20, parameter              :: coord="cmks"
  logical, save                        :: init_from_g4 = .false.  ! We provide the three-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put e.g. the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: R0_=1, h_=2, Rmax_=3, deltaR_=4, thetamax_=5
  integer, parameter                   :: ncoordpar = 5
  double precision, save               :: coordpar(ncoordpar)
  
  interface rks
     module procedure rks_d, rks_c
  end interface rks
  
  interface thetaks
     module procedure thetaks_d, thetaks_c
  end interface thetaks
  
  interface r
     module procedure r_d, r_c
  end interface r
  
  interface theta
     module procedure theta_d, theta_c
  end interface theta

contains
!=============================================================================
  subroutine init_coord
    
  end subroutine init_coord
  !=============================================================================
  double complex function r_c(s)

    double complex, intent(in)             :: s
    !-----------------------------------------------------------------------------

    r_c = coordpar(R0_) + exp(s)
    
  end function r_c
  !=============================================================================
  double precision function r_d(s)

    double precision, intent(in)            :: s
    !-----------------------------------------------------------------------------

    r_d = coordpar(R0_) + exp(s)
    
  end function r_d
  !=============================================================================
  double complex function theta_c(vartheta)

    double complex, intent(in)             :: vartheta
    !-----------------------------------------------------------------------------

    theta_c = vartheta + 0.5d0 * coordpar(h_) * sin(2.0d0*vartheta)
    
  end function theta_c
  !=============================================================================
  double precision function theta_d(vartheta)

    double precision, intent(in)             :: vartheta
    !-----------------------------------------------------------------------------

    theta_d = vartheta + 0.5d0 * coordpar(h_) * sin(2.0d0*vartheta)
    
  end function theta_d
  !=============================================================================
  double complex function rks_c(s,vartheta)

    double complex, intent(in)             :: s, vartheta
    ! .. local  ..
    double complex                         :: mytheta, myr
    double complex                         :: step, dstepdrho, dstepdrhoc
    double complex                         :: rhoc, drhocdthetaks
    !-----------------------------------------------------------------------------

    myr = r(s); mytheta = theta(vartheta)
    call get_rhoc(mytheta,rhoc,drhocdthetaks)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    rks_c = sqrt(myr**2*Cos(mytheta)**2+(rhoc+myr*step  &
         - 1.0d0*rhoc*step)**2*Sin(mytheta)**2)
    
  end function rks_c
  !=============================================================================
  double precision function rks_d(s,vartheta)

    double precision, intent(in)           :: s, vartheta
    ! .. local  ..
    double complex                         :: mytheta, myr
    double complex                         :: step, dstepdrho, dstepdrhoc
    double complex                         :: rhoc, drhocdthetaks
    !-----------------------------------------------------------------------------

    myr = r(s); mytheta = theta(vartheta)
    call get_rhoc(mytheta,rhoc,drhocdthetaks)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    rks_d = sqrt(myr**2*Cos(mytheta)**2+(rhoc+myr*step  &
         - 1.0d0*rhoc*step)**2*Sin(mytheta)**2)

  end function rks_d
  !=============================================================================
  double complex function thetaks_c(s,vartheta)

    !-----------------------------------------------------------------------------
    ! Returns the KS-theta from the current coordinates s, vartheta
    ! Since Atan2 does not work for complex arguments, need to identify the
    ! quadrant.
    ! Works as long as vartheta in [-pi/2,3pi/2]
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    
    double complex, intent(in)             :: s, vartheta
    ! .. local  ..
    double complex                         :: mytheta, myr
    double complex                         :: step, dstepdrho, dstepdrhoc
    double complex                         :: rhoc, drhocdthetaks
    !-----------------------------------------------------------------------------

    myr = r(s); mytheta = theta(vartheta)
    call get_rhoc(mytheta,rhoc,drhocdthetaks)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)

    thetaks_c = my_atan( ((rhoc+myr*step - 1.0d0*rhoc*step)*Sin(mytheta)) &
         / (myr*Cos(mytheta)) )

    if (real(vartheta) .gt. dpi/2.0d0 .and. real(vartheta) .le. dpi) then
       thetaks_c = thetaks_c + dpi
    else if (real(vartheta) .le. -dpi/2.0d0 .or. real(vartheta) .gt. dpi) then
       thetaks_c = thetaks_c - dpi
    end if
    
  end function thetaks_c
  !=============================================================================
  double precision function thetaks_d(s,vartheta)
    
    !-----------------------------------------------------------------------------
    ! Returns the KS-theta from the current coordinates s, vartheta
    ! Since Atan2 does not work for complex arguments, need to identify the
    ! quadrant.
    ! Works as long as vartheta in [-pi/2,3pi/2]
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    
    double precision, intent(in)           :: s, vartheta
    ! .. local  ..
    double complex                         :: mytheta, myr
    double complex                         :: step, dstepdrho, dstepdrhoc
    double complex                         :: rhoc, drhocdthetaks
    !-----------------------------------------------------------------------------

    myr = r(s); mytheta = theta(vartheta)
    call get_rhoc(mytheta,rhoc,drhocdthetaks)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)

    thetaks_d = atan(dble( ((rhoc+myr*step - 1.0d0*rhoc*step)*Sin(mytheta)) &
         / (myr*Cos(mytheta)) ))
    
    if (vartheta .gt. dpi/2.0d0 .and. vartheta .le. dpi) then
       thetaks_d = thetaks_d + dpi
    else if (vartheta .le. -dpi/2.0d0 .or. vartheta .gt. dpi) then
       thetaks_d = thetaks_d - dpi
    end if

  end function thetaks_d
  !=============================================================================
  subroutine KSToCoord_point(rKS,thetaKS,s,vartheta)

    !-----------------------------------------------------------------------------
    ! Returns the current coordinates s, vartheta
    ! from the KS coordinates rks, thetaks
    ! Need to do 2D Newton-Raphson.
    !-----------------------------------------------------------------------------

    use mod_rtsafe2D
    include 'amrvacdef.f'

    double precision, intent(in)           :: rKS, thetaKS
    double precision, intent(out)          :: s, vartheta
    ! .. local ..
    integer                                :: niter, ierror
    !-----------------------------------------------------------------------------

    ! Initial guess:
    s = log(rks - coordpar(R0_)); vartheta = thetaks

    MAXIT = 30
    call rtsafe2D(frks,fthetaks,0.0d0,bigdouble,-2.0d0*dpi,2.0d0*dpi,absaccnr,s,vartheta,niter,ierror)
    
    if (ierror .ne. 0) call mpistop('KSToCoord_point: Does not converge')

  contains
    !=============================================================================
    subroutine frks(s,vartheta,f,dfdx,dfdy,return_derivatives)
      double precision, intent(in)        :: s, vartheta
      double precision, intent(out)       :: f, dfdx, dfdy
      logical, intent(in)                 :: return_derivatives
      ! .. local ..
      double complex                         :: mytheta, myr
      double complex                         :: step, dstepdrho, dstepdrhoc
      double complex                         :: rhoc, drhocdtheta
      !-----------------------------------------------------------------------------

      myr = r(cmplx(s, 0.0d0, kind(1.d0))); mytheta = theta(cmplx(vartheta, 0.0d0, kind(1.d0)))
      call get_rhoc(mytheta,rhoc,drhocdtheta)
      call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)

      f = -1.0d0*rKS+sqrt(myr**2*Cos(mytheta)**2+(rhoc+myr*step  &
           - 1.0d0*rhoc*step)**2*Sin(mytheta)**2)

      if (return_derivatives) then 
         dfdx = (exp(s)*(myr*Cos(mytheta)**2+(dstepdrho*(myr-1.0d0*rhoc)+  &
              step)*(rhoc+myr*step-  &
              1.0d0*rhoc*step)*Sin(mytheta)**2))/Sqrt(myr**2*Cos(mytheta)**2  &
              +(rhoc+myr*step-1.0d0*rhoc*step)**2*Sin(mytheta)**2)

         dfdy = ((1+coordpar(h_)*Cos(2.0d0*vartheta))*Sin(mytheta)*((myr-  &
              1.0d0*rhoc)*(-1+step)*(myr+rhoc+myr*step-  &
              1.0d0*rhoc*step)*Cos(mytheta)+drhocdtheta*(1+dstepdrhoc*(myr  &
              -1.0d0*rhoc)-1.0d0*step)*(rhoc+myr*step-  &
              1.0d0*rhoc*step)*Sin(mytheta)))/Sqrt(myr**2*Cos(mytheta)**2+  &
              (rhoc+myr*step-1.0d0*rhoc*step)**2*Sin(mytheta)**2)
      end if

    end subroutine frks
    !=============================================================================
    subroutine fthetaks(s,vartheta,f,dfdx,dfdy,return_derivatives)
      double precision, intent(in)        :: s, vartheta
      double precision, intent(out)       :: f, dfdx, dfdy
      logical, intent(in)                 :: return_derivatives
      ! .. local ..
      double complex                         :: mytheta, myr
      double complex                         :: step, dstepdrho, dstepdrhoc
      double complex                         :: rhoc, drhocdtheta
      !-----------------------------------------------------------------------------

      myr = r(cmplx(s, 0.0d0, kind(1.d0))); mytheta = theta(cmplx(vartheta, 0.0d0, kind(1.d0)))
      call get_rhoc(mytheta,rhoc,drhocdtheta)
      call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)

      f = -1.0d0*thetaKS+atan2(dble( (rhoc+myr*step-  &
           1.0d0*rhoc*step)*Sin(mytheta) ), dble( myr*Cos(mytheta)) )

      if (return_derivatives) then
         dfdx = (exp(s)*(dstepdrho*myr*(myr-1.0d0*rhoc)+rhoc*(-1+  &
              step))*Cos(mytheta)*Sin(mytheta))/(myr**2*Cos(mytheta)**2+  &
              (rhoc*(-1+step)-1.0d0*myr*step)**2*Sin(mytheta)**2)

         dfdy = (myr*(1+coordpar(h_)*Cos(2.0d0*vartheta))*(2.0d0*(rhoc+  &
              myr*step-rhoc*step)+drhocdtheta*(1+dstepdrhoc*(myr-  &
              1.0d0*rhoc)-  &
              1.0d0*step)*Sin(2.0d0*mytheta)))/(2.0d0*(myr**2*Cos(mytheta)**2  &
              +(rhoc*(-1+step)-1.0d0*myr*step)**2*Sin(mytheta)**2))
      end if

    end subroutine fthetaks
    !=============================================================================
  end subroutine KSToCoord_point
  !=============================================================================
  subroutine get_rhoc(mytheta,rhoc,drhocdtheta)

    include 'amrvacdef.f'
    
    double complex, intent(in)             :: mytheta
    double complex, intent(out)            :: rhoc, drhocdtheta
    double precision, parameter            :: mythetamin=dpi/40.0d0
    !-----------------------------------------------------------------------------

    rhoc = (coordpar(Rmax_)*r(xprobmin1)*Sin(coordpar(thetamax_)))/(r(xprobmin1) &
         * Sin(coordpar(thetamax_))  &
         + (coordpar(Rmax_)  &
         - 1.0d0*r(xprobmin1))*sqrt(mythetamin**2  &
         + Sin(mytheta)**2))

    drhocdtheta = (coordpar(Rmax_)*r(xprobmin1)*(-1.0d0*coordpar(Rmax_)  &
         + r(xprobmin1))*Cos(mytheta)*Sin(mytheta)*Sin(coordpar(thetamax_))) &
         / (Sqrt(mythetamin**2  &
         + Sin(mytheta)**2)*(r(xprobmin1)*Sin(coordpar(thetamax_))  &
         + (coordpar(Rmax_)  &
         - 1.0d0*r(xprobmin1))*sqrt(mythetamin**2  &
         + Sin(mytheta)**2))**2)
    
  end subroutine get_rhoc
  !=============================================================================
  subroutine get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)

    double complex, intent(in)               :: myr, rhoc
    double complex, intent(out)              :: step, dstepdrho, dstepdrhoc
    double precision, parameter              :: epsilon = 1.0d-12
    !-----------------------------------------------------------------------------

    if (abs(1.0d0/(exp((myr-rhoc)/coordpar(deltaR_)+1.0d0))) .gt. epsilon) then

       step       = 1.0d0 - 1.0d0 / (exp((myr-rhoc)/coordpar(deltaR_)) + 1.0d0)
       
       dstepdrho  = 1.0d0/(2.0d0*coordpar(deltaR_)+2.0d0*coordpar(deltaR_)*my_cosh((myr  &
            - rhoc)/coordpar(deltaR_)))
       
       dstepdrhoc = 1.0d0/(-2.0d0*coordpar(deltaR_)-2.0d0*coordpar(deltaR_)*my_cosh((myr  &
            - rhoc)/coordpar(deltaR_)))

    else ! Need to give these to avoid floating cancelation error:
       
       step = 1.0d0
       dstepdrho = 0.0d0
       dstepdrhoc = 0.0d0
       
    end if
    
  end subroutine get_step
  !=============================================================================
  double complex function alpha(x^D)

    !-----------------------------------------------------------------------------
    ! This is the lapse.  It is a complex function because I use the
    ! complex step numerical derivative.
    ! Its best to use this like 
    ! myalpha = alpha({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
    ! where the position is converted to complex numbers as expected.
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    
    double complex, intent(in)             :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta)
    mythetaks = thetaks(s,vartheta)
    
    alpha = sqrt(1.0d0-(4.0d0*eqpar(m_)*myrks)/(eqpar(a_)**2+  &
         2.0d0*myrks*(2.0d0*eqpar(m_)+myrks)+  &
         eqpar(a_)**2*Cos(2.0d0*mythetaks)))

  end function alpha
  !=============================================================================
  double complex function betar(x^D)

    !-----------------------------------------------------------------------------
    ! This is the contravariant radial shift beta^r.
    ! It is a complex function because I use the
    ! complex step numerical derivative.
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    !
    double complex, intent(in)             :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks, myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta); mythetaks = thetaks(s,vartheta)
    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    betar = (2.0d0*eqpar(m_)*exp(-s)*myr*myrks*sqrt(myr**2*Cos(mytheta)**2  &
         + (rhoc*(-1.0d0+step)  &
         - 1.0d0*myr*step)**2*Sin(mytheta)**2)*(2*(rhoc+myr*step  &
         - rhoc*step)+drhocdtheta*(1+dstepdrhoc*(myr-rhoc)  &
         - step)*Sin(2*mytheta)))/((myrks*(2*eqpar(m_)+myrks)  &
         + eqpar(a_)**2*Cos(mythetaks)**2)*(-2*(dstepdrho*myr*(myr  &
         - rhoc)+rhoc*(-1  &
         + step))*Cos(mytheta)*Sin(mytheta)**2*(-(myr**2*Cos(mytheta))  &
         + (rhoc+myr*step-rhoc*step)*((rhoc+myr*step  &
         - rhoc*step)*Cos(mytheta)+drhocdtheta*(1+dstepdrhoc*(myr-rhoc)  &
         - step)*Sin(mytheta)))+myr*(myr*Cos(mytheta)**2  &
         + (dstepdrho*(myr-rhoc)+step)*(rhoc+myr*step  &
         - rhoc*step)*Sin(mytheta)**2)*(2*(rhoc+myr*step-rhoc*step)  &
         + drhocdtheta*(1+dstepdrhoc*(myr-rhoc)-step)*Sin(2*mytheta))))

  end function betar
  !=============================================================================
  double complex function betatheta(x^D)
    
    !-----------------------------------------------------------------------------
    ! This is the contravariant theta shift beta^\theta.
    ! It is a complex function because I use the
    ! complex step numerical derivative.
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    !
    double complex, intent(in)             :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks, myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta); mythetaks = thetaks(s,vartheta)
    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    betatheta = (-4.0d0*eqpar(m_)*myrks*(dstepdrho*myr*(myr-rhoc)+rhoc*(-1  &
         + step))*Cos(mytheta)*Sin(mytheta)*sqrt(myr**2*Cos(mytheta)**2  &
         + (rhoc*(-1.0d0+step)  &
         - 1.0d0*myr*step)**2*Sin(mytheta)**2))/((myrks*(2*eqpar(m_)  &
         + myrks)+eqpar(a_)**2*Cos(mythetaks)**2)*(1  &
         + coordpar(h_)*Cos(2*vartheta))*(-2*(dstepdrho*myr*(myr  &
         - rhoc)+rhoc*(-1  &
         + step))*Cos(mytheta)*Sin(mytheta)**2*(-(myr**2*Cos(mytheta))  &
         + (rhoc+myr*step-rhoc*step)*((rhoc+myr*step  &
         - rhoc*step)*Cos(mytheta)+drhocdtheta*(1+dstepdrhoc*(myr-rhoc)  &
         - step)*Sin(mytheta)))+myr*(myr*Cos(mytheta)**2  &
         + (dstepdrho*(myr-rhoc)+step)*(rhoc+myr*step  &
         - rhoc*step)*Sin(mytheta)**2)*(2*(rhoc+myr*step-rhoc*step)  &
         + drhocdtheta*(1+dstepdrhoc*(myr-rhoc)-step)*Sin(2*mytheta))))
    
  end function betatheta
  !=============================================================================
  double complex function gammarr(x^D)
    
    include 'amrvacdef.f'
    
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks, myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta); mythetaks = thetaks(s,vartheta)
    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    gammarr = (exp(2.0d0*s)*((dstepdrho*myr*(myr-1.0d0*rhoc)+rhoc*(-1+  &
         step))**2*Cos(mytheta)**2*(myrks**2+  &
         eqpar(a_)**2*Cos(mythetaks)**2)*Sin(mytheta)**2+(1+  &
         (2.0d0*eqpar(m_)*myrks)/(myrks**2+  &
         eqpar(a_)**2*Cos(mythetaks)**2))*(myr**2*Cos(mytheta)**2+  &
         (rhoc*(-1+step)-  &
         1.0d0*myr*step)**2*Sin(mytheta)**2)*(myr*Cos(mytheta)**2+  &
         (dstepdrho*(myr-1.0d0*rhoc)+step)*(rhoc+myr*step-  &
         1.0d0*rhoc*step)*Sin(mytheta)**2)**2))/(myr**2*Cos(mytheta)**2  &
         +(rhoc*(-1+step)-1.0d0*myr*step)**2*Sin(mytheta)**2)**2
    
  end function gammarr
  !=============================================================================
  double complex function gammartheta(x^D)
    
    include 'amrvacdef.f'
    
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks, myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta); mythetaks = thetaks(s,vartheta)
    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    gammartheta = (exp(s)*(1+  &
         coordpar(h_)*Cos(2.0d0*vartheta))*Sin(mytheta)*(2.0d0*(1+  &
         (2*eqpar(m_)*myrks)/(myrks**2+  &
         eqpar(a_)**2*Cos(mythetaks)**2))*(myr**2*Cos(mytheta)**2+  &
         (rhoc*(-1+step)-  &
         myr*step)**2*Sin(mytheta)**2)*(myr*Cos(mytheta)**2+  &
         (dstepdrho*(myr-rhoc)+step)*(rhoc+myr*step-  &
         rhoc*step)*Sin(mytheta)**2)*(-(myr**2*Cos(mytheta))+(rhoc+  &
         myr*step-rhoc*step)*((rhoc+myr*step-rhoc*step)*Cos(mytheta)+  &
         drhocdtheta*(1+dstepdrhoc*(myr-rhoc)-step)*Sin(mytheta)))+  &
         myr*(dstepdrho*myr*(myr-1.0d0*rhoc)+rhoc*(-1+  &
         step))*Cos(mytheta)*(myrks**2+  &
         eqpar(a_)**2*Cos(mythetaks)**2)*(2.0d0*(rhoc+myr*step-  &
         rhoc*step)+drhocdtheta*(1+dstepdrhoc*(myr-1.0d0*rhoc)-  &
         1.0d0*step)*Sin(2.0d0*mytheta))))/(2.0d0*(myr**2*Cos(mytheta)**2  &
         +(rhoc*(-1+step)-1.0d0*myr*step)**2*Sin(mytheta)**2)**2)
    
  end function gammartheta
  !=============================================================================
  double complex function gammarphi(x^D)
    
    include 'amrvacdef.f'
    
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks, myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta); mythetaks = thetaks(s,vartheta)
    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    gammarphi = (-1.0d0*eqpar(a_)*exp(s)*(1+  &
         (2*eqpar(m_)*myrks)/(myrks**2+  &
         eqpar(a_)**2*Cos(mythetaks)**2))*(myr*Cos(mytheta)**2+  &
         (dstepdrho*(myr-rhoc)+step)*(rhoc+myr*step-  &
         rhoc*step)*Sin(mytheta)**2)*Sin(mythetaks)**2) &
         /Sqrt(myr**2*Cos(mytheta)**2  &
         +(rhoc*(-1+step)-myr*step)**2*Sin(mytheta)**2)
    
  end function gammarphi
  !=============================================================================
  double complex function gammathetatheta(x^D)
    
    include 'amrvacdef.f'
    
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks, myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta); mythetaks = thetaks(s,vartheta)
    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)

    gammathetatheta = ((1+coordpar(h_)*Cos(2.0d0*vartheta))**2*(4.0d0*(1+  &
         (2*eqpar(m_)*myrks)/(myrks**2+  &
         eqpar(a_)**2*Cos(mythetaks)**2))*Sin(mytheta)**2*(myr**2*Cos(mytheta)**2  &
         +(rhoc*(-1+step)-  &
         myr*step)**2*Sin(mytheta)**2)*(myr**2*Cos(mytheta)-(rhoc+  &
         myr*step-rhoc*step)*((rhoc+myr*step-rhoc*step)*Cos(mytheta)+  &
         drhocdtheta*(1+dstepdrhoc*(myr-rhoc)-step)*Sin(mytheta)))**2  &
         +myr**2*(myrks**2+  &
         eqpar(a_)**2*Cos(mythetaks)**2)*(2.0d0*(rhoc+myr*step-  &
         rhoc*step)+drhocdtheta*(1+dstepdrhoc*(myr-1.0d0*rhoc)-  &
         1.0d0*step)*Sin(2.0d0*mytheta))**2))/(4.0d0*(myr**2*Cos(mytheta)**2  &
         +(rhoc*(-1+step)-1.0d0*myr*step)**2*Sin(mytheta)**2)**2)

  end function gammathetatheta
  !=============================================================================
  double complex function gammathetaphi(x^D)
    
    include 'amrvacdef.f'
    
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks, myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta); mythetaks = thetaks(s,vartheta)
    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    gammathetaphi = (-1.0d0*eqpar(a_)*(1+(2*eqpar(m_)*myrks)/(myrks**2  &
         + eqpar(a_)**2*Cos(mythetaks)**2))*(1  &
         + coordpar(h_)*Cos(2*vartheta))*Sin(mytheta)*(-(myr**2*Cos(mytheta))  &
         + (rhoc+myr*step-rhoc*step)*((rhoc+myr*step  &
         - rhoc*step)*Cos(mytheta)+drhocdtheta*(1+dstepdrhoc*(myr-rhoc)  &
         - step)*Sin(mytheta)))*Sin(mythetaks)**2)/Sqrt(myr**2*Cos(mytheta)**2 &
         + (rhoc*(-1+step)-myr*step)**2*Sin(mytheta)**2)
    
  end function gammathetaphi
  !=============================================================================
  double complex function gammaphiphi(x^D)
    
    include 'amrvacdef.f'
    
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: s, vartheta
    double complex :: myrks, mythetaks, myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    s = x1
    {^IFZIN
    vartheta = x^Z
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }

    myrks     = rks(s,vartheta); mythetaks = thetaks(s,vartheta)
    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)

    gammaphiphi = Sin(mythetaks)**2*(eqpar(a_)**2+myrks**2+  &
         (2.0d0*eqpar(a_)**2*eqpar(m_)*myrks*Sin(mythetaks)**2)/(myrks**2  &
         +eqpar(a_)**2*Cos(mythetaks)**2))

  end function gammaphiphi
  !=============================================================================
  subroutine get_sqrtgamma_analytic(x^D,sqrtgamma,is_analytic)

    !-----------------------------------------------------------------------------
    ! You can specify sqrtgamma here.
    ! If you don't want to, set is_analytic = .false.
    ! In that case, get_sqrtgamma() will calculate from the metric,
    ! which can be a major slowdown.
    ! init_metric() will set sqrtgamma_is_analytic from the return-value
    ! of "is_analytic".  
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    !-----------------------------------------------------------------------------

    if(present(is_analytic)) is_analytic = .false.

    sqrtgamma = 0.0d0 ! to avoid compiler-warnings.
    
  end subroutine get_sqrtgamma_analytic
  !=============================================================================
  subroutine get_alpha(x^D,myalpha,iszero,dalphadj_iszero,dalphadj,jdir)

    !-----------------------------------------------------------------------------
    ! This is the lapse.  Optional parameter is true if lapse is
    ! identically zero (does not really make sense)
    ! Optional parameters jdir and dalphadj request derivatives
    ! \partial_j \alpha ; j=jdir
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    double precision, intent(in)                     :: x^D
    integer, optional, intent(in)                    :: jdir
    double precision, intent(out)                    :: myalpha
    logical, optional, intent(out)                   :: iszero, dalphadj_iszero
    double precision, optional, intent(out)          :: dalphadj
    ! .. local ..
    !-----------------------------------------------------------------------------
    if(present(dalphadj) .and. .not. present(jdir) .or. &
         present(dalphadj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_alpha: derivatives requested without direction or output-slot given.")

    if(present(iszero)) iszero = .false.

    myalpha = alpha({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

    if (present(jdir)) then
       if (jdir .eq. 1 {^IFZIN .or. jdir .eq. ^Z} ) then
          if (present(dalphadj)) then
             call complex_derivative(x^D,alpha,jdir,dalphadj)
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .false.
          end if
       else
          if (present(dalphadj)) then
             dalphadj = 0.0d0
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .true.
          end if
       end if
    end if

  end subroutine get_alpha
  !=============================================================================
  subroutine get_beta(idir,x^D,mybeta,iszero,dbetaidj_iszero,dbetaidj,jdir)

    !-----------------------------------------------------------------------------
    ! This is the (contravariant!!) shift vector.
    ! The optional argument iszero is true if shift-component is 
    ! identically zero.
    ! if requested, dbetaidj is the derivative of the contravariant shift.
    ! \partial_j \beta^i ; i=idir, j=jdir
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    integer, intent(in)                      :: idir
    double precision, intent(in)             :: x^D
    integer, optional, intent(in)            :: jdir
    double precision, intent(out)            :: mybeta
    logical, optional, intent(out)           :: iszero, dbetaidj_iszero
    double precision, optional, intent(out)  :: dbetaidj
    ! .. local ..
    !-----------------------------------------------------------------------------
    if(present(dbetaidj) .and. .not. present(jdir) .or. &
         present(dbetaidj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_beta: derivatives requested &
         &without direction or output-slot given.")

    if (idir .eq. 1) then

       !==============================
       ! betar component
       !==============================

       if(present(iszero)) iszero = .false.
       mybeta = betar({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

       if (present(jdir)) then
          if (jdir .eq. 1) then ! dbetar/dr
             if (present(dbetaidj)) call complex_derivative(x^D,betar,jdir,dbetaidj)
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .false.
          else if (jdir .eq. ^Z) then ! dbetar/dtheta
             {^IFZIN
             if (present(dbetaidj)) call complex_derivative(x^D,betar,jdir,dbetaidj)
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .false.
             }{^IFZOUT
             if (present(dbetaidj)) dbetaidj = 0.0d0
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.
             }
          else ! dbetar/dphi
             if (present(dbetaidj)) dbetaidj = 0.0d0
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.
          end if
       end if

    else if (idir .eq. ^Z) then

       !==============================
       ! betatheta component:
       !==============================

       if(present(iszero)) iszero = .false.
       mybeta = betatheta({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
       
       if (present(jdir)) then
          if (jdir .eq. 1) then ! dbetatheta/dr
             if (present(dbetaidj)) call complex_derivative(x^D,betatheta,jdir,dbetaidj)
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .false.
          else if (jdir .eq. ^Z) then ! dbetatheta/dtheta
             {^IFZIN
             if (present(dbetaidj)) call complex_derivative(x^D,betatheta,jdir,dbetaidj)
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .false.
             }{^IFZOUT
             if (present(dbetaidj)) dbetaidj = 0.0d0
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.
             }
          else ! dbetatheta/dphi
             if (present(dbetaidj)) dbetaidj = 0.0d0
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.
          end if
       end if
    else 

       !==============================
       ! betaphi component:
       !==============================

       mybeta = zero
       if(present(iszero))           iszero          = .true.
       if (present(dbetaidj))        dbetaidj        = 0.0d0
       if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.

    end if

  end subroutine get_beta
  !=============================================================================
  subroutine get_g_component(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir)

    !-----------------------------------------------------------------------------
    ! This is at the heart of the scheme: Set the three-metric gamma.
    ! Indices of the metric are down (covariant) g_{i j}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dgdk and kdir request derivatives of the metric
    ! \partial_k g_{ij} ; i=iin, j=jin, k=kdir
    ! The optional argument dgdk_iszero is true if the element is identically zero.
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: g
    logical, optional, intent(out)           :: iszero, dgdk_iszero
    double precision, optional, intent(out)  :: dgdk
    ! .. local ..
    integer                                  :: i,j
    !-----------------------------------------------------------------------------
    if(present(dgdk) .and. .not. present(kdir) .or. &
         present(dgdk_iszero) .and. .not. present(kdir)) &
         call mpistop("get_g_component: derivatives requested without &
         &direction or output-slot given.")

    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    select case(i)
    case(1)
       ! ==============================
       ! r-something coordinate
       ! ==============================
       select case(j)
       case(1) ! rr component
          
          if(present(iszero))       iszero      = .false.
          g = gammarr({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammarr,kdir,dgdk)
             else 
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if
          
       case(^Z) ! rtheta component
          
          if(present(iszero))       iszero      = .false.
          g = gammartheta({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammartheta,kdir,dgdk)
             else
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if
          
       case(^PHI) ! rphi component

          if(present(iszero))       iszero      = .false.
          g = gammarphi({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammarphi,kdir,dgdk)
             else
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if
          
       end select

    case(^Z)
       ! ==============================
       ! theta-something coordinate
       ! ==============================
       select case(j)
       case(1) ! thetar component
          if(present(iszero))       iszero      = .false.
          g = gammartheta({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammartheta,kdir,dgdk)
             else 
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if

       case(^Z) ! thetatheta component

          if(present(iszero))       iszero      = .false.
          g = gammathetatheta({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammathetatheta,kdir,dgdk)
             else
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if

       case(^PHI) ! thetaphi component

          if(present(iszero))       iszero      = .false.
          g = gammathetaphi({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammathetaphi,kdir,dgdk)
             else
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if

       end select
       
    case(^PHI)
       ! ==============================
       ! phi-something coordinate
       ! ==============================
       select case(j)
       case(1) ! phir component
          if(present(iszero))       iszero      = .false.
          g = gammarphi({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammarphi,kdir,dgdk)
             else 
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if

       case(^Z) ! phitheta component

          if(present(iszero))       iszero      = .false.
          g = gammathetaphi({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammathetaphi,kdir,dgdk)
             else
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if

       case(^PHI) ! phiphi component

          if(present(iszero))       iszero      = .false.
          g = gammaphiphi({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

          if (present(kdir)) then
             if (kdir .eq. 1 {^IFZIN .or. kdir .eq. ^Z} ) then
                if(present(dgdk_iszero))  dgdk_iszero = .false.
                if (present(dgdk)) call complex_derivative(x^D,gammaphiphi,kdir,dgdk)
             else
                if(present(dgdk_iszero))  dgdk_iszero = .true.
                if (present(dgdk)) dgdk = 0.0d0
             end if
          end if

       end select
    end select

    {^IFZOUT ! We are actually in the r-phi plane and dtheta is zero:
    if (present(kdir)) then
       if (kdir .eq. ^Z) then
          if(present(dgdk_iszero))  dgdk_iszero = .true.
          dgdk = 0.0d0
       end if
    end if
    }

  end subroutine get_g_component
  !=============================================================================
  double precision function outerhorizon()

    !-----------------------------------------------------------------------------
    ! Get the postition of the horizon on the equatorial plane.
    ! vartheta=thetaKS=pi/2
    ! Need to do 1D Newton-Raphson
    !-----------------------------------------------------------------------------

    use mod_rtsafe2D, only: MAXIT, rt1D
    include 'amrvacdef.f'
    
    ! .. local ..
    double precision                                       :: rksout
    integer                                                :: ierror
    !-----------------------------------------------------------------------------
    
    rksout       = eqpar(m_) + sqrt(eqpar(m_)**2 - eqpar(a_)**2)
    outerhorizon = log(rksout-coordpar(R0_))

    MAXIT = 30
    call rt1D(func_outerhorizon,outerhorizon,absaccnr,ierror)

    if (ierror .ne. 0) call mpistop('outerhorizon: Does not converge')
    
  contains
    !=============================================================================
    subroutine func_outerhorizon(s,f,df)

      double precision, intent(in)      :: s
      double precision, intent(out)     :: f, df
      ! .. local ..
      double complex     :: rhoc, drhocdthetaks
      double complex     :: step, dstepdrho, dstepdrhoc
      !-----------------------------------------------------------------------------

      call get_rhoc(cmplx(dpi/2.0d0, 0.0d0, kind(1.d0)),rhoc,drhocdthetaks)
      call get_step(cmplx(r(s), 0.0d0, kind(1.d0)),rhoc,step,dstepdrho,dstepdrhoc)

      f  = -rksout+abs(-1.0d0*rhoc*(-1.0d0+step)+r(s)*step)
      df = (exp(s)*(dstepdrho*(r(s)-1.0d0*rhoc)+step)*(rhoc+  &
           r(s)*step-1.0d0*rhoc*step))/abs(rhoc+r(s)*step-  &
           1.0d0*rhoc*step)
      
    end subroutine func_outerhorizon
    !=============================================================================
  end function outerhorizon
  !=============================================================================
  subroutine BLToCoord(ixI^L,ixO^L,xBL,xCoord)

    !-----------------------------------------------------------------------------
    ! Transforms the Boyer-Lindquist spatial coordinates to the current
    ! cylindrified modified KS coordinates.
    !-----------------------------------------------------------------------------

    use mod_transform, only: BLToKS
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    !-----------------------------------------------------------------------------

    call BLToKS(ixI^L,ixO^L,xBL,xKS)
    call KSToCoord(ixI^L,ixO^L,xKS,xCoord)

  end subroutine BLToCoord
  !=============================================================================
  subroutine KSToCoord(ixI^L,ixO^L,xKS,xCoord)

    !-----------------------------------------------------------------------------
    ! Transforms the Kerr-Schild spatial coordinates to the current
    ! cylindrified modified KS coordinates.
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    integer                                                :: ix^D
    double precision                                       :: s, vartheta
    !-----------------------------------------------------------------------------

    {do ix^DB=ixOmin^DB,ixOmax^DB\}

      call KSToCoord_point(xKS(ix^D,1),{^IFZIN xKS(ix^D,^Z)} {^IFZOUT dpi/2.0d0}, &
           s,vartheta)

      xCoord(ix^D,1) = s
      {^IFZIN
      xCoord(ix^D,^Z) = vartheta
      }
      {^IFPHIIN
      xCoord(ix^D,^PHI) = xKS(ix^D,^PHI)
      }

    {end do\}

  end subroutine KSToCoord
  !=============================================================================
  subroutine CoordToKS(ixI^L,ixO^L,xCoord,xKS)

    !-----------------------------------------------------------------------------
    ! Transforms the cylindrified modified KS coordinates to
    ! ordinary Kerr-Schild coordinates.
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xKS
    ! .. local ..
    integer                                                :: ix^D
    double precision                                       :: s, vartheta
    !-----------------------------------------------------------------------------

    {do ix^DB=ixOmin^DB,ixOmax^DB\}

      s = xCoord(ix^D,1)
      {^IFZIN
      vartheta = xCoord(ix^D,^Z)
      }{^IFZOUT
      vartheta = dpi/2.0d0
      }

      xKS(ix^D,1) = rks(s,vartheta)
      {^IFZIN
      xKS(ix^D,^Z) = thetaks(s,vartheta)
      }
      {^IFPHIIN
      xKS(ix^D,^PHI) = xCoord(ix^D,^PHI)
      }

    {end do\}

  end subroutine CoordToKS
  !=============================================================================
  subroutine CoordToBL(ixI^L,ixO^L,xCoord,xBL)

    use mod_transform, only: KSToBL
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xBL
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    !-----------------------------------------------------------------------------

    call CoordToKS(ixI^L,ixO^L,xCoord,xKS)
    call KSToBL(ixI^L,ixO^L,xKS,xBL)

  end subroutine CoordToBL
  !=============================================================================
  subroutine u4BLtoCoord(ixI^L,ixO^L,xBL,u4BL,u4Coord)

    !-----------------------------------------------------------------------------
    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to the current coordinates u4Coord.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    ! Position also given in BL coordinates.
    !-----------------------------------------------------------------------------
    
    use mod_transform, only: u4BLToKS, BLToKS
    include 'amrvacdef.f'

    integer                                                :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    double precision, dimension(ixI^S,0:^NC)               :: u4KS
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    !-----------------------------------------------------------------------------

    call u4BLtoKS(ixI^L,ixO^L,xBL,u4BL,u4KS)
    call BLToKS(ixI^L,ixO^L,xBL,xKS)
    
    call u4KStoCoord(ixI^L,ixO^L,xKS,u4KS,u4Coord)
    
  end subroutine u4BLtoCoord
  !=============================================================================
  subroutine u4KStoCoord(ixI^L,ixO^L,xKS,u4KS,u4Coord,J)

    !-----------------------------------------------------------------------------
    ! Transforms the (contravariant) four-velocity u4KS from Kerr-Schild coordinates
    ! to the current coordinates u4Coord.  
    ! Position also given in KS coordinates.
    !-----------------------------------------------------------------------------
    
    include 'amrvacdef.f'

    integer                                                :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4KS
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xCoord
    double precision, dimension(ixI^S,0:^NC,0:^NC)         :: Jac
    integer                                                :: ix^D, ix, jx
    double precision                                       :: s, vartheta
    double precision                                       :: rks, thetaks
    double complex                                         :: mytheta, myr
    double complex                                         :: step, dstepdrho, dstepdrhoc
    double complex                                         :: rhoc, drhocdtheta
    !-----------------------------------------------------------------------------

    ! Get the current coordinates from KS:
    call KSToCoord(ixI^L,ixO^L,xKS,xCoord)
    
    ! Assemble the Jacobian:
    Jac            = 0.0d0
    
    Jac(ixO^S,0,0) = 1.0d0

    {^IFPHI
    Jac(ixO^S,^PHI,^PHI) = 1.0d0
    }

    {do ix^DB=ixOmin^DB,ixOmax^DB\}

      s = xCoord(ix^D,1)
      {^IFZIN
      vartheta = xCoord(ix^D,^Z)
      }{^IFZOUT
      vartheta = dpi/2.0d0
      }
    
      myr = r(cmplx(s, 0.0d0, kind(1.d0))); mytheta = theta(cmplx(vartheta, 0.0d0, kind(1.d0)))
      call get_rhoc(mytheta,rhoc,drhocdtheta)
      call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
      Jac(ix^D,1,1) = (exp(-1.0d0*s)*myr*sqrt(myr**2*Cos(mytheta)**2+(rhoc*(-1.0d0  &
           +step)-1.0d0*myr*step)**2*Sin(mytheta)**2)*(2.0d0*(rhoc+  &
           myr*step-rhoc*step)+drhocdtheta*(1+dstepdrhoc*(myr-  &
           1.0d0*rhoc)-  &
           1.0d0*step)*Sin(2.0d0*mytheta)))/(-2.0d0*(dstepdrho*myr*(myr  &
           -rhoc)+rhoc*(-1+  &
           step))*Cos(mytheta)*Sin(mytheta)**2*(-(myr**2*Cos(mytheta))+  &
           (rhoc+myr*step-rhoc*step)*((rhoc+myr*step-  &
           rhoc*step)*Cos(mytheta)+drhocdtheta*(1+dstepdrhoc*(myr-rhoc)  &
           -step)*Sin(mytheta)))+myr*(myr*Cos(mytheta)**2+  &
           (dstepdrho*(myr-1.0d0*rhoc)+step)*(rhoc+myr*step-  &
           1.0d0*rhoc*step)*Sin(mytheta)**2)*(2.0d0*(rhoc+myr*step-  &
           rhoc*step)+drhocdtheta*(1+dstepdrhoc*(myr-1.0d0*rhoc)-  &
           1.0d0*step)*Sin(2.0d0*mytheta)))

      {^IFZ
      Jac(ix^D,1,^Z) = (-2.0d0*exp(-s)*Sin(mytheta)*(myr**2*Cos(mytheta)**2+  &
           (rhoc*(-1+step)-  &
           myr*step)**2*Sin(mytheta)**2)*(-(myr**2*Cos(mytheta))+(rhoc+  &
           myr*step-rhoc*step)*((rhoc+myr*step-rhoc*step)*Cos(mytheta)+  &
           drhocdtheta*(1+dstepdrhoc*(myr-rhoc)-  &
           step)*Sin(mytheta))))/(-2*(dstepdrho*myr*(myr-rhoc)+rhoc*(-1  &
           +step))*Cos(mytheta)*Sin(mytheta)**2*(-(myr**2*Cos(mytheta))  &
           +(rhoc+myr*step-rhoc*step)*((rhoc+myr*step-  &
           rhoc*step)*Cos(mytheta)+drhocdtheta*(1+dstepdrhoc*(myr-rhoc)  &
           -step)*Sin(mytheta)))+myr*(myr*Cos(mytheta)**2+  &
           (dstepdrho*(myr-rhoc)+step)*(rhoc+myr*step-  &
           rhoc*step)*Sin(mytheta)**2)*(2*(rhoc+myr*step-rhoc*step)+  &
           drhocdtheta*(1+dstepdrhoc*(myr-rhoc)-step)*Sin(2*mytheta)))

      Jac(ix^D,^Z,1) = (-2.0d0*(dstepdrho*myr*(myr-rhoc)+rhoc*(-1+  &
           step))*Cos(mytheta)*Sin(mytheta)*sqrt(myr**2*Cos(mytheta)**2  &
           +(rhoc*(-1.0d0+step)-  &
           1.0d0*myr*step)**2*Sin(mytheta)**2))/((1+  &
           coordpar(h_)*Cos(2*vartheta))*(-2*(dstepdrho*myr*(myr-  &
           rhoc)+rhoc*(-1+  &
           step))*Cos(mytheta)*Sin(mytheta)**2*(-(myr**2*Cos(mytheta))+  &
           (rhoc+myr*step-rhoc*step)*((rhoc+myr*step-  &
           rhoc*step)*Cos(mytheta)+drhocdtheta*(1+dstepdrhoc*(myr-rhoc)  &
           -step)*Sin(mytheta)))+myr*(myr*Cos(mytheta)**2+  &
           (dstepdrho*(myr-rhoc)+step)*(rhoc+myr*step-  &
           rhoc*step)*Sin(mytheta)**2)*(2*(rhoc+myr*step-rhoc*step)+  &
           drhocdtheta*(1+dstepdrhoc*(myr-rhoc)-step)*Sin(2*mytheta))))

      Jac(ix^D,^Z,^Z) = (-1.0d0*(rhoc*(dstepdrho*rhoc-step)*(-1+step)+  &
           dstepdrho*myr**2*step+myr*(1+dstepdrho*rhoc-  &
           2*dstepdrho*rhoc*step+step**2)-(rhoc*(dstepdrho*rhoc-  &
           step)*(-1+step)+dstepdrho*myr**2*step+myr*(-1+dstepdrho*rhoc  &
           -2*dstepdrho*rhoc*step+step**2))*Cos(2*mytheta)))/((1+  &
           coordpar(h_)*Cos(2*vartheta))*(-(dstepdrho*myr**2)-rhoc+  &
           dstepdrho*myr*rhoc-2*myr*step+rhoc*step+(dstepdrho*myr*(myr-  &
           rhoc)+rhoc*(-1+step))*Cos(2*mytheta)-drhocdtheta*(1+  &
           dstepdrhoc*myr-dstepdrhoc*rhoc-step)*Sin(2*mytheta)))
      }

      {end do\}
          
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled, now
    ! transform contravariant four-vector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    u4Coord(ixO^S,:) = zero
    do jx=0,ndir
       do ix=0,ndir
          u4Coord(ixO^S,ix) = u4Coord(ixO^S,ix) + Jac(ixO^S,ix,jx) * u4KS(ixO^S,jx)
       end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(J)) J = Jac

    
  end subroutine u4KStoCoord
  !=============================================================================
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)

    !-----------------------------------------------------------------------------
    ! Transforms any contravariant three vector u3Coord in the current coordinates
    ! to a vector in the Cartesian basis.
    ! 
    ! u3Cart = J u3Coord
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^NC), intent(in)   :: u3Coord
    double precision, dimension(ixI^S,1:^NC), intent(out)  :: u3Cart
    double precision, dimension(ixI^S,1:^NC,1:^NC), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixI^S,1:^NC,1:^NC)         :: Jac
    integer                                                :: zcart_, ycart_, ix,jx
    integer                                                :: ix^D
    double precision   :: s, vartheta, sph, cph
    double complex     :: myr, mytheta, rhoc, drhocdtheta, myrks, mythetaks
    double complex     :: step, dstepdrho, dstepdrhoc
    !-----------------------------------------------------------------------------

    if (ndim.eq.3.and.ndir.eq.3) then
       ycart_=2
       zcart_=3
    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. ^PHI .eq. 2) then
       ycart_=2
       zcart_=3
    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. ^Z .eq. 2) then
       ycart_=3
       zcart_=2
    else if (ndir.eq.2 .and. ^Z.eq.2) then
       ycart_=1 ! must catch ycart_=1 case and don't use
       zcart_=2
    else if (ndir.eq.2 .and. ^PHI.eq.2) then
       ycart_=2
       zcart_=1 ! must catch zcart_=1 case and don't use
    else if (ndir.eq.1) then
       ycart_=1
       zcart_=1
    else
       call mpistop("u3CoordToCart: unknown parameter combination of -d, -phi and -z!")
    end if

 
    Jac(ixO^S,:,:) = zero
    {do ix^DB=ixOmin^DB,ixOmax^DB\}

    s = x(ix^D,1)
    {^IFZIN
    vartheta = x(ix^D,^Z)
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }
    {^IFPHIIN
    sph = sin(x(ix^D,phi_))
    cph = cos(x(ix^D,phi_))
    }{^IFPHIOUT
    sph = zero
    cph = one
    }

    myr       = r(s);            mytheta   = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)
    
    ! Jxr:
    Jac(ix^D,1,1) = cph*(dstepdrho*exp(s)*myr-1.0d0*dstepdrho*exp(s)*rhoc+  &
         exp(s)*step)*Sin(mytheta)

    {^IFZ
    ! Jxtheta:
    Jac(ix^D,1,^Z) = -1.0d0*cph*(1+coordpar(h_)*Cos(2*vartheta))*((rhoc*(-1+  &
         step)-myr*step)*Cos(mytheta)+drhocdtheta*(-1+  &
         dstepdrhoc*(-myr+rhoc)+step)*Sin(mytheta))
    }
    {^IFPHI
    ! Jxphi:
    Jac(ix^D,1,^PHI) = sph*(rhoc*(-1+step)-1.0d0*myr*step)*Sin(mytheta)
    }
    
    if (ycart_.ne.1) then
       ! Jyr:
       Jac(ix^D,ycart_,1) = exp(s)*sph*(dstepdrho*(myr-1.0d0*rhoc)+step)*Sin(mytheta)
       {^IFZ
       ! Jytheta:
       Jac(ix^D,ycart_,^Z) = -1.0d0*sph*(1+coordpar(h_)*Cos(2*vartheta))*((rhoc*(-1+  &
            step)-myr*step)*Cos(mytheta)+drhocdtheta*(-1+  &
            dstepdrhoc*(-myr+rhoc)+step)*Sin(mytheta))
       }
       {^IFPHI
       ! Jyphi:
       Jac(ix^D,ycart_,^PHI) = cph*(rhoc+myr*step-1.0d0*rhoc*step)*Sin(mytheta)
       }       
    end if

    if (zcart_.ne.1) then
       ! Jzr:
       Jac(ix^D,zcart_,1) = exp(s)*Cos(mytheta)
       {^IFZ
       ! Jztheta:
       Jac(ix^D,zcart_,^Z) = -1.0d0*myr*(1+coordpar(h_)*Cos(2*vartheta))*Sin(mytheta)
       }
       {^IFPHI
       ! Jzphi: Always zero
       Jac(ix^D,zcart_,^PHI) = 0.0d0
       }
    end if
    
    {end do\}
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled, now
    ! transform contravariant three-vector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    u3Cart(ixO^S,:) = zero
    do ix=1,ndir
       do jx=1,ndir
          u3Cart(ixO^S,ix) = u3Cart(ixO^S,ix) + Jac(ixO^S,ix,jx) * u3Coord(ixO^S,jx)
       end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (present(J)) J = Jac

  end subroutine u3CoordToCart
  !=============================================================================
  subroutine CoordToCart(ixI^L,ixO^L,x,xCart)

    include 'amrvacdef.f'
    
    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCart
    ! .. local ..
    integer                                                :: zcart_, ycart_
    integer                                                :: ix^D
    double complex :: s, vartheta
    double complex :: myr, mytheta
    double complex :: step, dstepdrho, dstepdrhoc, rhoc, drhocdtheta
    double precision                                       :: sph, cph
    !-----------------------------------------------------------------------------

    if (ndim.eq.3) then
       ycart_=2
       zcart_=3
    else if (ndim.eq.2 .and. ^Z.eq.2) then
       ycart_=1 ! must catch ycart_=1 case and don't use
       zcart_=2
    else if (ndim.eq.2 .and. ^PHI.eq.2) then
       ycart_=2
       zcart_=1 ! must catch zcart_=1 case and don't use
    else if (ndim.eq.1) then
       ycart_=1
       zcart_=1
    else
       call mpistop("CoordToCart: unknown parameter combination of -d, -phi and -z!")
    end if
    !-----------------------------------------------------------------------------

    {do ix^DB=ixOmin^DB,ixOmax^DB\}

    s = x(ix^D,1)
    {^IFZIN
    vartheta = x(ix^D,^Z)
    }{^IFZOUT
    vartheta = dpi/2.0d0
    }
    {^IFPHIIN
    sph = sin(x(ix^D,phi_))
    cph = cos(x(ix^D,phi_))
    }{^IFPHIOUT
    sph = zero
    cph = one
    }

    myr     = r(s)
    mytheta = theta(vartheta)

    call get_rhoc(mytheta,rhoc,drhocdtheta)
    call get_step(myr,rhoc,step,dstepdrho,dstepdrhoc)

    xCart(ix^D,1) = (-rhoc*(-1.0d0+step)+myr*step)*Sin(mytheta)*cph
    if (ycart_.ne.1) &
         xCart(ix^D,ycart_) = (-1.0d0*rhoc*(-1+step)+myr*step)*Sin(mytheta)*sph
    if (zcart_.ne.1) &
         xCart(ix^D,zcart_) = {^IFZIN myr*Cos(mytheta)}{^IFZOUT 0.0d0}

    {end do\}
    
  end subroutine CoordToCart
  !=============================================================================
  ! Auxilliary functions
  !=============================================================================
  double complex function my_cosh(x)
    double complex, intent(in)       :: x
    !-----------------------------------------------------------------------------
    
    my_cosh = ( exp(x) + exp(-x) ) / 2
    
  end function my_cosh
  !=============================================================================
  double complex function my_atan(z)
    double complex, intent(in)       :: z
    ! .. local ..
    double complex, parameter        :: I = (0.0d0,1.0d0)
    !-----------------------------------------------------------------------------
    
    my_atan = 0.5d0*I * ( log(1.0d0 - I*z) - log(1.0d0 + I*z) )
    
  end function my_atan
  !=============================================================================
  ! DUMMIES
  !=============================================================================
  subroutine get_dgammainvdk(x^D,dgammainvdk,k)

    !-----------------------------------------------------------------------------
    ! Obtains the derivatives of the inverse _spatial_ metric
    ! \partial_k gamma^{ij} for a given k. 
    ! This is required to obtain derivatives e.g. of the contravariant shift
    ! and the lapse.  
    !-----------------------------------------------------------------------------

    use mod_lu
    include 'amrvacdef.f'

    double precision, intent(in)           :: x^D
    double precision, intent(out)          :: dgammainvdk(1:^NC,1:^NC)
    integer, intent(in)                    :: k
    ! .. local ..
    !-----------------------------------------------------------------------------

    call mpistop('get_dgammainvdk: not implemented and not needed!')

    dgammainvdk = 0.0d0 ! To avoid compiler warnings.
    
  end subroutine get_dgammainvdk
  !=============================================================================
  ! End of coordinate-specific definitions.
  !=============================================================================
