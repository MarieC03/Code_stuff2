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
  ! Metric components for modified Kerr-Schild coordinates -coord=mks
  !
  ! Check modKS.nb for the definitions.
  ! Oliver Porth
  ! 02.02.2017
  !=============================================================================
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------

  character*20, parameter              :: coord="mksnorm"
  logical, save                        :: init_from_g4 = .false.  ! We don't provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: ncoordpar = 0
  double precision, save               :: coordpar(ncoordpar)
  ! not yet implemented...

contains
!=============================================================================
  subroutine init_coord
    
  end subroutine init_coord

  !=============================================================================
  double complex function s(rks)

    include 'amrvacdef.f'
    double complex                               :: rks
    ! .. local ..
    !-----------------------------------------------------------------------------

    s = log( rks - eqpar(R0_) )

  end function s
  !=============================================================================
  double complex function vartheta(theta)

    include 'amrvacdef.f'
    double complex                               :: theta
    ! .. local ..
    double complex                               :: rtmp
    !-----------------------------------------------------------------------------

    if (eqpar(h_) .eq. zero) then 

       vartheta = theta
       return

    else

       rtmp = (9.0d0*eqpar(h_)**2*(dpi-2.0d0*theta)+  &
            sqrt(3.0d0)*eqpar(h_)*sqrt(-1.0d0*eqpar(h_)*((-4.0d0+  &
            eqpar(h_))*(dpi+2.0d0*dpi*eqpar(h_))**2+  &
            1.080d2*dpi*eqpar(h_)*theta-  &
            1.080d2*eqpar(h_)*theta**2)))**3.3333333333333330d-1

       vartheta = (dpi**6.6666666666666660d-1 &
            * (6.0d0*dpi**3.333333333333333d-1  &
            -(5.2414827884177940d0*dpi**6.6666666666666660d-1*(-1.0d0+  &
            eqpar(h_)))/rtmp-  &
            (2.28942848510666330d0*rtmp)/eqpar(h_)))/1.20d1

    end if

  end function vartheta
  !=============================================================================
  double complex function rks(s)

    include 'amrvacdef.f'
    double complex                               :: s
    ! .. local ..
    !-----------------------------------------------------------------------------

    rks = eqpar(R0_) + exp(s)

  end function rks
  !=============================================================================
  double complex function thetaks(vartheta)

    include 'amrvacdef.f'
    double complex                               :: vartheta

    ! .. local ..
    !-----------------------------------------------------------------------------

    thetaks = vartheta+(2.0d0*eqpar(h_)*vartheta*(-2.0d0*vartheta+  &
         dpi)*(-vartheta+dpi))/dpi**2

  end function thetaks
  !=============================================================================
  subroutine get_sqrtgamma_analytic(x^D,sqrtgamma,is_analytic)

    include 'amrvacdef.f'

    ! You can specify sqrtgamma here.
    ! If you don't want to, set is_analytic = .false.
    ! In that case, get_sqrtgamma() will calculate from the metric,
    ! which can be a major slowdown.
    ! init_metric() will set the flag sqrtgamma_is_analytic from the return-value
    ! of this subroutine.  

    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    ! .. local ..
    {^IFZIN
    double precision                                 :: mythetaks}
    !-----------------------------------------------------------------------------

    if(present(is_analytic)) is_analytic = .true.

    {^IFZIN
    mythetaks = thetaks(cmplx(x^Z, 0.0d0, kind(1.d0)))
    }
    
    {^IFZIN
    sqrtgamma = sqrt(exp(2.0d0*x1)*(dpi**2*(1.0d0+2.0d0*eqpar(h_))  &
         - 1.20d1*dpi*eqpar(h_)*x^Z  &
         + 1.20d1*eqpar(h_)*x^Z**2)**2*((exp(x1)  &
         + eqpar(R0_))*(exp(x1)+2.0d0*eqpar(m_)+eqpar(R0_))  &
         + eqpar(a_)**2*Cos(mythetaks)**2)*Sin(mythetaks)**2*(eqpar(a_)**2  &
         + (exp(x1)+eqpar(R0_))**2  &
         - 1.0d0*eqpar(a_)**2*Sin(mythetaks)**2))/dpi**2
    }{^IFZOUT
    sqrtgamma = sqrt(exp(2.0d0*x1)*(-1.0d0+eqpar(h_))**2*(exp(x1)  &
         + eqpar(R0_))**3*(exp(x1)+2.0d0*eqpar(m_)  &
         + eqpar(R0_)))
    }


  end subroutine get_sqrtgamma_analytic
  !=============================================================================
  subroutine get_alpha(x^D,myalpha,iszero,dalphadj_iszero,dalphadj,jdir)

    include 'amrvacdef.f'

    ! get the lapse.  Optional parameter is true if lapse is
    ! identically zero (does not really make sense)
    ! Optional parameters jdir and dalphadj request derivatives
    ! \partial_j \alpha ; j=jdir
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

    ! Derivative info requested
    if (present(jdir)) then
       if (jdir .eq. 1) then
          if (present(dalphadj)) call complex_derivative(x^D,alpha,jdir,dalphadj)
          if(present(dalphadj_iszero)) dalphadj_iszero = .false.
          {^IFZIN
       else if (jdir .eq. ^Z) then
          if (present(dalphadj)) call complex_derivative(x^D,alpha,jdir,dalphadj)
          if(present(dalphadj_iszero)) dalphadj_iszero = .false.      
          }
       else 
          if (present(dalphadj)) dalphadj = 0.0d0
          if (present(dalphadj_iszero)) dalphadj_iszero = .true.
       end if
    end if


  contains
    !=============================================================================
    double complex function alpha(x^D)
      ! this is the lapse.  It is a complex function because I use the
      ! complex step numerical derivative.
      ! Its best to use this like 
      ! myalpha = alpha({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
      ! where the position is converted to complex numbers as expected.
      !
      double complex, intent(in)             :: x^D
      double complex                         :: myrks, mythetaks
      !-----------------------------------------------------------------------------

      myrks = rks(x1)
      {^IFZIN
      mythetaks = thetaks(x^Z)
      }
      {^IFZOUT
      mythetaks = dpi/2.0d0
      }  

      {^IFZIN
      alpha = sqrt((myrks**2+  &
           eqpar(a_)**2*Cos(mythetaks)**2)/(myrks*(2.0d0*eqpar(m_)  &
           +myrks)+eqpar(a_)**2*Cos(mythetaks)**2))
      }{^IFZOUT
      alpha =sqrt(myrks/(2.0d0*eqpar(m_)+myrks))
      }

    end function alpha
    !=============================================================================
  end subroutine get_alpha
  !=============================================================================
  subroutine get_beta(idir,x^D,beta,iszero,dbetaidj_iszero,dbetaidj,jdir)

    include 'amrvacdef.f'

    ! get the (contravariant!!) shift vector.
    ! The optional argument iszero is true if shift-component is 
    ! identically zero.
    ! if requested, dbetaidj is the derivative of the contravariant shift.
    ! \partial_j \beta^i ; i=idir, j=jdir
    integer, intent(in)                      :: idir
    double precision, intent(in)             :: x^D
    integer, optional, intent(in)            :: jdir
    double precision, intent(out)            :: beta
    logical, optional, intent(out)           :: iszero, dbetaidj_iszero
    double precision, optional, intent(out)  :: dbetaidj
    ! .. local ..
    !-----------------------------------------------------------------------------
    if(present(dbetaidj) .and. .not. present(jdir) .or. &
         present(dbetaidj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_beta: derivatives requested &
         &without direction or output-slot given.")

    select case(idir)

    case(1)
       ! Betar
       beta = betar({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
       if(present(iszero)) iszero = .false.

       ! Derivative info requested
       if (present(jdir)) then
          if (jdir .eq. 1) then
             if (present(dbetaidj)) call complex_derivative(x^D,betar,jdir,dbetaidj)
             if(present(dbetaidj_iszero)) dbetaidj_iszero = .false.
             {^IFZIN
          else if (jdir .eq. ^Z) then
             if (present(dbetaidj)) call complex_derivative(x^D,betar,jdir,dbetaidj)
             if(present(dbetaidj_iszero)) dbetaidj_iszero = .false.      
             }
          else 
             if (present(dbetaidj)) dbetaidj = 0.0d0
             if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.
          end if
       end if

    case default
       ! All zeroes except r:
       beta = 0.0d0
       if(present(iszero))           iszero          = .true.
       if (present(dbetaidj))         dbetaidj       = 0.0d0
       if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.

    end select

  contains
    !=============================================================================
    double complex function betar(x^D)
      !
      double complex, intent(in)             :: x^D
      double complex                         :: myrks, mythetaks
      !-----------------------------------------------------------------------------

      myrks = rks(x1)
      {^IFZIN
      mythetaks = thetaks(x^Z)
      }
      {^IFZOUT
      mythetaks = dpi/2.0d0
      }  

      {^IFZIN
      betar = (2.0d0*eqpar(m_)*myrks)/(exp(x1)*(myrks*(2*eqpar(m_)  &
           +myrks)+eqpar(a_)**2*Cos(mythetaks)**2))
      }{^IFZOUT
      betar = (2.0d0*eqpar(m_))/(exp(x1)*(2*eqpar(m_)+myrks))
      }
      
    end function betar
    !=============================================================================
  end subroutine get_beta
  !=============================================================================
  subroutine get_g_component(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir)

    include 'amrvacdef.f'

    ! This is at the heart of the scheme: Set the (spatial) metric components here
    ! and only here...
    ! Indices of the metric are down (covariant) g_{ij}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dgdk and kdir request derivatives of the metric
    ! \partial_k g_{ij} ; i=iin, j=jin, k=kdir
    ! The optional argument dgdk_iszero
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
    if (iin .gt. jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    if(present(iszero)) iszero = .false.

    ! Non-diagonal!
    if (i .eq. 1 .and. j .eq. 1) then
       !grr:
       g = grr({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

       ! Derivative info requested
       if (present(kdir)) then
          if (kdir .eq. 1) then
             if (present(dgdk)) call complex_derivative(x^D,grr,kdir,dgdk)
             if(present(dgdk_iszero)) dgdk_iszero = .false.
             {^IFZIN
          else if (kdir .eq. ^Z) then
             if (present(dgdk)) call complex_derivative(x^D,grr,kdir,dgdk)
             if(present(dgdk_iszero)) dgdk_iszero = .false.                
             }
          else 
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.
          end if
       end if

    else if (i .eq. ^PHI .and. j .eq. ^PHI) then
       ! gphiphi
       g = gphph({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

       
       ! Derivative info requested
       if (present(kdir)) then
          if (kdir .eq. 1) then
             if (present(dgdk)) call complex_derivative(x^D,gphph,kdir,dgdk)
             if(present(dgdk_iszero)) dgdk_iszero = .false.
             {^IFZIN
          else if (kdir .eq. ^Z) then
             if (present(dgdk)) call complex_derivative(x^D,gphph,kdir,dgdk)
             if(present(dgdk_iszero)) dgdk_iszero = .false.                
             }
          else 
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.
          end if
       end if

       
    else if (i .eq. ^Z .and. j .eq. ^Z) then
       ! gthetatheta
       g = gthth({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})


       ! Derivative info requested
       if (present(kdir)) then
          if (kdir .eq. 1) then
             if (present(dgdk)) call complex_derivative(x^D,gthth,kdir,dgdk)
             if(present(dgdk_iszero)) dgdk_iszero = .false.
             {^IFZIN
          else if (kdir .eq. ^Z) then
             if (present(dgdk)) call complex_derivative(x^D,gthth,kdir,dgdk)
             if(present(dgdk_iszero)) dgdk_iszero = .false.                
             }
          else 
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.
          end if
       end if

       
    else if (i .eq. 1 .and. j .eq. ^PHI) then
       ! grphi
       g = grph({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
 

       ! Derivative info requested
       if (present(kdir)) then
          if (kdir .eq. 1) then
             if (present(dgdk)) call complex_derivative(x^D,grph,kdir,dgdk)
             if(present(dgdk_iszero)) dgdk_iszero = .false.
             {^IFZIN
          else if (kdir .eq. ^Z) then
             if (present(dgdk)) call complex_derivative(x^D,grph,kdir,dgdk)
             if(present(dgdk_iszero)) dgdk_iszero = .false.                
             }
          else 
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.
          end if
       end if

    else

       g = 0.0d0

       if(present(iszero))      iszero       = .true.
       if (present(dgdk))       dgdk         = 0.0d0
       if(present(dgdk_iszero)) dgdk_iszero  = .true.

    end if
  contains
    !=============================================================================
    double complex function grr(x^D)
      double complex, intent(in)          :: x^D
      double complex                      :: myrks, mythetaks
      !-----------------------------------------------------------------------------
      
      myrks = rks(x1)
      {^IFZIN
      mythetaks = thetaks(x^Z)
      }
      {^IFZOUT
      mythetaks = dpi/2.0d0
      }  


      {^IFZIN
      grr = exp(2.0d0*x1)*(1.0d0  &
           + (2.0d0*eqpar(m_)*myrks)/(myrks**2  &
           + eqpar(a_)**2*Cos(mythetaks)**2))
      }{^IFZOUT
      grr = (exp(2.0d0*x1)*(2.0d0*eqpar(m_)+myrks))/myrks
      }

    end function grr
    !=============================================================================
    double complex function gthth(x^D)
      double complex, intent(in)          :: x^D
      double complex                      :: myrks, mythetaks
      !-----------------------------------------------------------------------------

      myrks = rks(x1)
      {^IFZIN
      mythetaks = thetaks(x^Z)
      }
      {^IFZOUT
      mythetaks = dpi/2.0d0
      }  

      {^IFZIN
      gthth = (1.0d0+(2.0d0*eqpar(h_)*(dpi**2-6.0d0*dpi*x^Z  &
           + 6.0d0*x^Z**2))/dpi**2)**2*(myrks**2  &
           + eqpar(a_)**2*Cos(mythetaks)**2)
      }{^IFZOUT
      gthth = (-1.0d0+eqpar(h_))**2*myrks**2
      }

    end function gthth
    !=============================================================================
    double complex function gphph(x^D)
      double complex, intent(in)          :: x^D
      double complex                      :: myrks, mythetaks
      !-----------------------------------------------------------------------------

      myrks = rks(x1)
      {^IFZIN
      mythetaks = thetaks(x^Z)
      }
      {^IFZOUT
      mythetaks = dpi/2.0d0
      }  

      {^IFZIN
      gphph = Sin(mythetaks)**2*(eqpar(a_)**2+myrks**2  &
           + (2.0d0*eqpar(a_)**2*eqpar(m_)*myrks*Sin(mythetaks)**2)/(myrks**2  &
           + eqpar(a_)**2*Cos(mythetaks)**2))
      }{^IFZOUT
      ! equatorial plane
      gphph = myrks**2+(eqpar(a_)**2*(2.0d0*eqpar(m_)  &
           + myrks))/myrks
      }

    end function gphph
    !=============================================================================
    double complex function grph(x^D)
      double complex, intent(in)          :: x^D
      double complex                      :: myrks, mythetaks
      !-----------------------------------------------------------------------------

      myrks = rks(x1)
      {^IFZIN
      mythetaks = thetaks(x^Z)
      }
      {^IFZOUT
      mythetaks = dpi/2.0d0
      }  

      {^IFZIN
      grph = -1.0d0*exp(x1)*eqpar(a_)*(1.0d0  &
           + (2.0d0*eqpar(m_)*myrks)/(myrks**2  &
           + eqpar(a_)**2*Cos(mythetaks)**2))*Sin(mythetaks)**2
      }{^IFZOUT
      grph = (-1.0d0*exp(x1)*eqpar(a_)*(2.0d0*eqpar(m_)  &
           + myrks))/myrks
      }

    end function grph
    !=============================================================================
  end subroutine get_g_component
  !=============================================================================
  double precision function outerhorizon()
    
    include 'amrvacdef.f'
    double precision                       :: rksout
    !-----------------------------------------------------------------------------

    rksout = eqpar(m_) + sqrt(eqpar(m_)**2 - eqpar(a_)**2)
    outerhorizon = s(cmplx(rksout, 0.0d0, kind(1.d0)))
    
  end function outerhorizon
  !=============================================================================
  subroutine u4BLtoCoord(ixI^L,ixO^L,xBL,u4BL,u4Coord)

    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to the current (modified KS) coordinates u4Coord.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    !
    !
    ! utKS       = utBL + 2 M rBL / Delta * urBL
    ! urKS       = e**-s * urBL
    ! uthetaKS   = uthetaBL * dpi**2/(1.20d1*eqpar(h_)*varth(ixO^S)**2-  &
    ! 1.20d1*eqpar(h_)*varth(ixO^S)*dpi+(1+2.0d0*eqpar(h_))*dpi**2)
    ! uphiKS     = uphiBL + a/Delta * urBL
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    integer                                                :: ix^D
    double precision, dimension(ixI^S)                     :: Delta, twoMr, ems, varth
    !-----------------------------------------------------------------------------

    ems(ixO^S)   = one/(xBL(ixO^S,1) - eqpar(R0_))

    {^IFZIN
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    varth(ix^D) = vartheta(cmplx(xBL(ix^D,^Z), 0.0d0, kind(1.d0)))
    {end do \}
    }{^IFZOUT
    varth(ixO^S) = dpi/2.0d0
    }

    twoMr(ixO^S) = 2.0d0*eqpar(m_)*abs(xBL(ixO^S,1))
    Delta(ixO^S) = xBL(ixO^S,1)**2 - twoMr(ixO^S) + eqpar(a_)**2

    u4Coord(ixO^S,0) = u4BL(ixO^S,0) + twoMr(ixO^S)/Delta(ixO^S) * u4BL(ixO^S,1)
    u4Coord(ixO^S,1) = ems(ixO^S) * u4BL(ixO^S,1)

    {^IFZ
    u4Coord(ixO^S,^Z) = dpi**2/(1.20d1*eqpar(h_)*varth(ixO^S)**2-  &
         1.20d1*eqpar(h_)*varth(ixO^S)*dpi+(1+  &
         2.0d0*eqpar(h_))*dpi**2) * u4BL(ixO^S,^Z)
    }
    {^IFPHI
    u4Coord(ixO^S,^PHI) = u4BL(ixO^S,^PHI) + eqpar(a_)/Delta(ixO^S) * u4BL(ixO^S,1)
    }

  end subroutine u4BLtoCoord
  !=============================================================================
  subroutine BLToCoord(ixI^L,ixO^L,xBL,xCoord)

    ! Transforms the (spatial) Boyer-Lindquist coordinates to the current coordinates
    ! We don't return the new time-component although it changes.
    !
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
  subroutine CoordToBL(ixI^L,ixO^L,xCoord,xBL)

    ! Transforms the (spatial) modified KS coordinate to Boyer-Lindquist coordinates.
    ! We don't return the new time-component although it changes.
    !
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
  subroutine KSToCoord(ixI^L,ixO^L,xKS,xCoord)

    ! Transforms the KS coordinate to the modified KS-coordinates.
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    integer                                                :: ix^D
    !-----------------------------------------------------------------------------

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    xCoord(ix^D,1) = s(cmplx(xKS(ix^D,1), 0.0d0, kind(1.d0)))
    {^IFZIN
    xCoord(ix^D,^Z) = vartheta(cmplx(xKS(ix^D,^Z), 0.0d0, kind(1.d0)))
    }
    {^IFPHIIN
    xCoord(ix^D,^PHI) = xKS(ix^D,^PHI)
    }
    {end do\}

  end subroutine KSToCoord
  !=============================================================================
  subroutine CoordToKS(ixI^L,ixO^L,xCoord,xKS)

    ! Transforms the modified KS coordinate to KS-coordinates.
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xKS
    ! .. local ..
    integer                                                :: ix^D
    !-----------------------------------------------------------------------------

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    xKS(ix^D,1) = rks(cmplx(xCoord(ix^D,1), 0.0d0, kind(1.d0)))
    {^IFZIN
    xKS(ix^D,^Z) = thetaks(cmplx(xCoord(ix^D,^Z), 0.0d0, kind(1.d0)))
    }
    {^IFPHIIN
    xKS(ix^D,^PHI) = xCoord(ix^D,^PHI)
    }
    {end do\}

  end subroutine CoordToKS
  !=============================================================================
  subroutine CoordToCart(ixI^L,ixO^L,x,xCart)

    ! Transforms the coordinates to "Kerr-Schild Cartesian coordinates"
    ! Two transformations are needed:  
    ! rKS     = R0 + exp(s)
    ! thetaKS = vartheta + 2h/pi**2 vartheta (pi-2vartheta) (pi-vartheta)
    ! phiKS   = phiMKS
    ! 
    ! x = rKS Sin(thetaKS) Cos(phiKS)
    ! y = rKS Sin(thetaKS) Sin(phiKS)
    ! z = rKS Cos(thetaKS)
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCart
    ! .. local ..
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, ra
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    integer                                                :: zcart_, ycart_
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

    call CoordToKS(ixI^L,ixO^L,x,xKS)

    {^IFPHIIN
    sph(ixO^S) = sin(xKS(ixO^S,phi_))
    cph(ixO^S) = cos(xKS(ixO^S,phi_))
    }{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    }

    {^IFZIN
    sth(ixO^S) = sin(xKS(ixO^S,z_))
    cth(ixO^S) = cos(xKS(ixO^S,z_))
    }{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    }

    ra(ixO^S) = xKS(ixO^S,1)!sqrt(xKS(ixO^S,1)**2+eqpar(a_)**2)

    xCart(ixO^S,1) = ra(ixO^S) * sth(ixO^S)*cph(ixO^S)
    if (ycart_.ne.1) &
         xCart(ixO^S,ycart_) = ra(ixO^S) * sth(ixO^S)*sph(ixO^S)
    if (zcart_.ne.1) &
         xCart(ixO^S,zcart_) = ra(ixO^S) * cth(ixO^S)

  end subroutine CoordToCart
  !=============================================================================
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)

    ! Transforms any contravariant three vector u3Coord in the current coordinates
    ! to a vector in the basis of "Kerr-Schild Cartesian coordinates".
    !
    ! Two transformations are needed:  
    ! rKS     = R0 + exp(s)
    ! thetaKS = vartheta + 2h/pi**2 vartheta (pi-2vartheta) (pi-vartheta)
    ! phiKS   = phiMKS
    ! 
    ! x = rKS Sin(thetaKS) Cos(phiKS)
    ! y = rKS Sin(thetaKS) Sin(phiKS)
    ! z = rKS Cos(thetaKS)
    ! 
    ! The transformation matrix J=del(x,y,z)/del(s,vartheta,phiMKS)
    ! Since the coordinates don't depend on time, this is the same transformation
    ! as for the four-vector u4CoordToCart()
    !
    ! u3Cart = J u3Coord
    !
    ! see modKS.nb for details
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^NC), intent(in)   :: u3Coord
    double precision, dimension(ixI^S,1:^NC), intent(out)  :: u3Cart
    double precision, dimension(ixI^S,1:^NC,1:^NC), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixI^S,1:^NC,1:^NC)         :: Jac
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, tmp
    integer                                                :: zcart_, ycart_, ix, jx
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

    call CoordToKS(ixI^L,ixO^L,x,xKS)

    {^IFPHIIN
    sph(ixO^S) = sin(xKS(ixO^S,phi_))
    cph(ixO^S) = cos(xKS(ixO^S,phi_))
    }{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    }

    {^IFZIN
    sth(ixO^S) = sin(xKS(ixO^S,z_))
    cth(ixO^S) = cos(xKS(ixO^S,z_))
    }{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    }


    tmp(ixO^S) = (dpi**2*(1.0d0+2.0d0*eqpar(h_))  &
         - 1.20d1*dpi*eqpar(h_)*x(ixO^S,^Z)+1.20d1*eqpar(h_)*x(ixO^S,^Z)**2)/dpi**2

    Jac(ixO^S,1,1) = exp(x(ixO^S,1)) * sth(ixO^S)*cph(ixO^S)
    {^IFZ
    Jac(ixO^S,1,^Z) = (cph(ixO^S)*cth(ixO^S)*tmp(ixO^S)*xKS(ixO^S,1))
    }
    {^IFPHI
    Jac(ixO^S,1,^PHI) = - sph(ixO^S)*sth(ixO^S)*xKS(ixO^S,1)
    }

    if (ycart_.ne.1) then
       Jac(ixO^S,ycart_,1) = exp(x(ixO^S,1))*sph(ixO^S)*sth(ixO^S)
       {^IFZ
       Jac(ixO^S,ycart_,^ZZ) = cth(ixO^S)*sph(ixO^S)*tmp(ixO^S)*xKS(ixO^S,1)
       }
       {^IFPHI
       Jac(ixO^S,ycart_,^PPHI) = cph(ixO^S)*sth(ixO^S)*xKS(ixO^S,1)
       }
    end if

    if(zcart_.ne.1) then
       Jac(ixO^S,zcart_,1) = cth(ixO^S)*exp(x(ixO^S,1))
       {^IFZ
       Jac(ixO^S,zcart_,^ZZ) = - sth(ixO^S)*tmp(ixO^S)*xKS(ixO^S,1)
       }
       {^IFPHI
       Jac(ixO^S,zcart_,^PPHI) = 0.0d0
       }
    end if

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
  ! Dummies
  !=============================================================================
  subroutine get_dgammainvdk(x^D,dgammainvdk,k)
    ! Obtains the derivatives of the inverse _spatial_ metric
    ! \partial_k gamma^{ij} for a given k. 
    ! This is not required when initialising gamma, lapse and shift
    ! e.g. when set init_from_g4 = .false.

    double precision, intent(in)           :: x^D
    double precision, intent(out)          :: dgammainvdk(1:^NC,1:^NC)
    integer, intent(in)                    :: k
    ! .. local ..
    !-----------------------------------------------------------------------------

    call mpistop("get_dgammainvdk: Not required and not implemented.")

    dgammainvdk = 0.0d0
    
  end subroutine get_dgammainvdk
  !=============================================================================
  ! End of coordinate-specific definitions.
  !=============================================================================
