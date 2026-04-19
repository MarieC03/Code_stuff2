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
  ! Metric components for Kerr-Schild coordinates -coord=ks
  !
  ! Check KS-coordinates.nb for the definitions.
  ! Oliver Porth
  ! 07.12.2015
  !=============================================================================
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------


  character*20, parameter              :: coord="ks"
  logical, save                        :: init_from_g4 = .false.  ! We don't provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: ncoordpar = 0
  double precision, save               :: coordpar(ncoordpar) 

contains
!=============================================================================
  subroutine init_coord

    LRNFu => KS_LRNFu

    u4CoordToKS => KS_u4CoordToKS

  end subroutine init_coord
  !=============================================================================
  subroutine get_sqrtgamma_analytic(x^D,sqrtgamma,is_analytic)

    include 'amrvacdef.f'

    ! You can specify sqrtgamma here.
    ! If you don't want to, set is_analytic = .false.
    ! In that case, get_sqrtgamma() will calculate from the metric,
    ! which can be a major slowdown.
    ! init_metric() will set sqrtgamma_is_analytic from the return-value
    ! of "is_analytic".  

    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    !-----------------------------------------------------------------------------

    if(present(is_analytic)) is_analytic = .true.


    {^IFZIN
    sqrtgamma = sqrt((abs(x1)*(2.0d0*eqpar(m_)+abs(x1))  &
         + eqpar(a_)**2*Cos(x^Z)**2)*Sin(x^Z)**2*(eqpar(a_)**2  &
         + abs(x1)**2-1.0d0*eqpar(a_)**2*Sin(x^Z)**2))
    }{^IFZOUT
    sqrtgamma = sqrt(abs(x1)**3*(2.0d0*eqpar(m_)+abs(x1)))
    }

  end subroutine get_sqrtgamma_analytic
  !=============================================================================
  subroutine get_alpha(x^D,alpha,iszero,dalphadj_iszero,dalphadj,jdir)

    include 'amrvacdef.f'

    ! get the lapse.  Optional parameter is true if lapse is
    ! identically zero (does not really make sense)
    ! Optional parameters jdir and dalphadj request derivatives
    ! \partial_j \alpha ; j=jdir
    double precision, intent(in)                     :: x^D
    integer, optional, intent(in)                    :: jdir
    double precision, intent(out)                    :: alpha
    logical, optional, intent(out)                   :: iszero, dalphadj_iszero
    double precision, optional, intent(out)          :: dalphadj
    ! .. local ..
    !-----------------------------------------------------------------------------
    if(present(dalphadj) .and. .not. present(jdir) .or. &
         present(dalphadj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_alpha: derivatives requested without direction or output-slot given.")

    if(present(iszero)) iszero = .false.

    {^IFZIN  
    alpha = 1.0d0/Sqrt(1.0d0+(2.0d0*eqpar(m_)*x1)/(x1**2+  &
         eqpar(a_)**2*Cos(x^Z)**2))
    }{^IFZOUT
    alpha = 1.0d0/Sqrt(1.0d0+(2.0d0*eqpar(m_))/x1)
    }

    if (present(jdir)) then

       select case(jdir)

       case(1)
          ! Radial derivative:
          if (present(dalphadj)) then
             {^IFZIN
             dalphadj = (2.0d0*eqpar(m_)*(-eqpar(a_)**2+2.0d0*x1**2-  &
                  eqpar(a_)**2*Cos(2*x^Z)))/((1+(2*eqpar(m_)*x1)/(x1**2+  &
                  eqpar(a_)**2*Cos(x^Z)**2))**1.50d0*(eqpar(a_)**2+  &
                  2.0d0*x1**2+eqpar(a_)**2*Cos(2*x^Z))**2)
             }{^IFZOUT
             dalphadj = eqpar(m_)/((1.0d0+(2.0d0*eqpar(m_))/x1)**1.50d0*x1**2)
             }
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .false.
          end if

       case(^Z)
          ! Theta-derivative:
          {^IFZIN
          if (present(dalphadj)) then
             dalphadj = (-2.0d0*eqpar(a_)**2*eqpar(m_)*x1*Cos(x^Z)*Sin(x^Z))/((x1**2  &
                  +eqpar(a_)**2*Cos(x^Z)**2)**2*(1.0d0+  &
                  (2.0d0*eqpar(m_)*x1)/(x1**2+  &
                  eqpar(a_)**2*Cos(x^Z)**2))**1.50d0)
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .false.
          end if
          }{^IFZOUT
          if (present(dalphadj)) then
             dalphadj = 0.0d0
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .true.
          end if
          }

       case default
          ! Other cases (phi):
          if (present(dalphadj)) then
             dalphadj = 0.0d0
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .true.
          end if
       end select

    end if ! present(jdir)

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
       {^IFZIN
       beta = (2.0d0*eqpar(m_)*x1)/(2.0d0*eqpar(m_)*x1+x1**2+  &
            eqpar(a_)**2*Cos(x^Z)**2)
       }{^IFZOUT
       beta = (2.0d0*eqpar(m_)*x1)/(2.0d0*eqpar(m_)*x1+x1**2)
       }

       if(present(iszero)) then
          iszero = .false.
       end if

       if (present(jdir)) then
          select case(jdir)

          case(1)
             ! dbetardr:
             if (present(dbetaidj)) then
                {^IFZIN
                dbetaidj = (-2.0d0*eqpar(m_)*(x1**2-  &
                     eqpar(a_)**2*Cos(x^Z)**2))/(x1*(2.0d0*eqpar(m_)+x1)+  &
                     eqpar(a_)**2*Cos(x^Z)**2)**2
                }{^IFZOUT
                dbetaidj = (-2.0d0*eqpar(m_))/(2.0d0*eqpar(m_)+x1)**2
                }
             end if ! present(dbetaidj)

             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .false.
             end if
          case(^Z)
             ! dbetardtheta:
             {^IFZIN
             if (present(dbetaidj)) then
                dbetaidj = (4.0d0*eqpar(a_)**2*eqpar(m_)*x1*Cos(x^Z)*Sin(x^Z))/(x1*(2*eqpar(m_) &
                     +x1)+eqpar(a_)**2*Cos(x^Z)**2)**2
             end if
             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .false.
             end if
             }{^IFZOUT
             if (present(dbetaidj)) then
                dbetaidj = 0.0d0
             end if
             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .true.
             end if
             }


          case default
             ! dbbetardphi is zero
             if (present(dbetaidj)) then
                dbetaidj = 0.0d0
             end if ! present(dbetaidj)           

             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .true.
             end if

          end select
       end if ! present(jdir)


    case default
       ! All zeroes except r:
       beta = 0.0d0
       if(present(iszero))           iszero          = .true.
       if (present(dbetaidj))         dbetaidj       = 0.0d0
       if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.

    end select

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
       {^IFZIN
       g = 1.0d0+(2.0d0*eqpar(m_)*x1)/(x1**2+  &
            eqpar(a_)**2*Cos(x^Z)**2)
       }{^IFZOUT
       g = 1.0d0+(2.0d0*eqpar(m_))/x1
       }
       if (present(kdir)) then
          ! Derivative info requested

          select case(kdir)

          case (1)
             ! dgrrdr:
             if (present(dgdk)) then
                {^IFZIN
                dgdk = (-4.0d0*eqpar(m_)*x1**2)/(x1**2+  &
                     eqpar(a_)**2*Cos(x^Z)**2)**2+(2.0d0*eqpar(m_))/(x1**2+  &
                     eqpar(a_)**2*Cos(x^Z)**2)
                }{^IFZOUT
                dgdk = (-2.0d0*eqpar(m_))/x1**2
                }
             end if
             if (present(dgdk_iszero)) dgdk_iszero = .false.

          case(^Z)
             !dgrrdtheta:
             {^IFZIN
             if (present(dgdk)) then
                dgdk = (4.0d0*eqpar(a_)**2*eqpar(m_)*x1*Cos(x^Z)*Sin(x^Z))/(x1**2  &
                     +eqpar(a_)**2*Cos(x^Z)**2)**2
             end if
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) then
                dgdk = 0.0d0
             end if
             if (present(dgdk_iszero)) dgdk_iszero = .true.
             }

          case default
             ! phi derivatives zero:
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.

          end select
       end if

    else if (i .eq. ^PHI .and. j .eq. ^PHI) then
       ! gphiphi
       {^IFZIN
       g = Sin(x^Z)**2*(eqpar(a_)**2+x1**2+  &
            (2.0d0*eqpar(a_)**2*eqpar(m_)*x1*Sin(x^Z)**2)/(x1**2+  &
            eqpar(a_)**2*Cos(x^Z)**2))
       }{^IFZOUT
       ! equatorial plane
       g = eqpar(a_)**2+(2.0d0*eqpar(a_)**2*eqpar(m_))/x1+x1**2
       }

       if (present(kdir)) then
          ! Derivative info requested
          select case(kdir)
          case(1)
             !DgphiphiDr
             {^IFZIN
             if (present(dgdk)) dgdk = 2.0d0*Sin(x^Z)**2*(x1+(eqpar(a_)**2*eqpar(m_)*(-x1**2+  &
                  eqpar(a_)**2*Cos(x^Z)**2)*Sin(x^Z)**2)/(x1**2+  &
                  eqpar(a_)**2*Cos(x^Z)**2)**2)
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) dgdk = 2.0d0*(-((eqpar(a_)**2*eqpar(m_))/x1**2)+x1)
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }
          case(^Z)
             !DgphiphiDtheta
             {^IFZIN
             if (present(dgdk)) dgdk = (eqpar(a_)**2+x1*(-2.0d0*eqpar(m_)+x1)+  &
                  (8.0d0*eqpar(m_)*x1*(eqpar(a_)**2+  &
                  x1**2)**2)/(eqpar(a_)**2+2.0d0*x1**2+  &
                  eqpar(a_)**2*Cos(2.0d0*x^Z))**2)*Sin(2.0d0*x^Z)
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.
             }
          case default
             !DphiphiDphi
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.
          end select
       end if

    else if (i .eq. ^Z .and. j .eq. ^Z) then
       ! gthetatheta
       {^IFZIN
       g = x1**2+eqpar(a_)**2*Cos(x^Z)**2
       }{^IFZOUT
       g = x1**2
       }

       if (present(kdir)) then
          ! Derivative info requested
          select case(kdir)

          case(1)
             ! DgthetathetaDr:
             if (present(dgdk)) dgdk = 2.0d0*x1
             if (present(dgdk_iszero)) dgdk_iszero = .false.

          case(^Z)
             ! DthetathetaDtheta:
             {^IFZIN
             if (present(dgdk)) dgdk = -2.0d0*eqpar(a_)**2*Cos(x^Z)*Sin(x^Z)
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.              
             }

          case default
             ! DthetathetaDphi
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.

          end select
       end if

    else if (i .eq. 1 .and. j .eq. ^PHI) then
       ! grphi
       {^IFZIN
       g = -1.0d0*eqpar(a_)*(1.0d0+(2.0d0*eqpar(m_)*x1)/(x1**2+  &
            eqpar(a_)**2*Cos(x^Z)**2))*Sin(x^Z)**2
       }{^IFZOUT
       g = -1.0d0*eqpar(a_)*(1.0d0+(2.0d0*eqpar(m_))/x1)
       }

       if (present(kdir)) then
          ! Derivative info requested
          select case(kdir)

          case(1)
             ! DgrphiDr:
             {^IFZIN
             if (present(dgdk)) dgdk = (2.0d0*eqpar(a_)*eqpar(m_)*(x1**2-  &
                  eqpar(a_)**2*Cos(x^Z)**2)*Sin(x^Z)**2)/(x1**2+  &
                  eqpar(a_)**2*Cos(x^Z)**2)**2
             }{^IFZOUT
             if (present(dgdk)) dgdk = (2.0d0*eqpar(a_)*eqpar(m_))/x1**2
             }
             if (present(dgdk_iszero)) dgdk_iszero = .false.

          case(^Z)
             ! DgrphiDtheta:
             {^IFZIN
             if (present(dgdk)) dgdk = eqpar(a_)*(-1.0d0-(8.0d0*eqpar(m_)*x1*(eqpar(a_)**2+  &
                  x1**2))/(eqpar(a_)**2+2.0d0*x1**2+  &
                  eqpar(a_)**2*Cos(2*x^Z))**2)*Sin(2.0d0*x^Z)
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.              
             }

          case default
             ! DrphiDphi
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.

          end select
       end if

    else

       g = 0.0d0

       if(present(iszero))      iszero       = .true.
       if (present(dgdk))       dgdk         = 0.0d0
       if(present(dgdk_iszero)) dgdk_iszero  = .true.

    end if

  end subroutine get_g_component
  !=============================================================================
    double precision function outerhorizon()
    
    include 'amrvacdef.f'
    !-----------------------------------------------------------------------------

    outerhorizon = eqpar(m_) + sqrt(eqpar(m_)**2 - eqpar(a_)**2)
    
  end function outerhorizon
  !=============================================================================
  subroutine u4BLtoCoord(ixI^L,ixO^L,xBL,u4BL,u4Coord)

    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to the current (KS) coordinates u4Coord.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    !
    ! See also Font et al., Mon. Not. R. Astron. Soc. 305, 920±936 (1999) for the
    ! transformations.
    !
    ! utKS       = utBL + 2 M rBL / Delta * urBL
    ! uphiKS     = uphiBL + a/Delta * urBL
    ! urKS       = urBL
    ! uthetaKS   = uthetaBL
    !
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    double precision, dimension(ixI^S)                     :: Delta, twoMr
    !-----------------------------------------------------------------------------

    twoMr(ixO^S) = 2.0d0*eqpar(m_)*abs(xBL(ixO^S,1))
    Delta(ixO^S) = xBL(ixO^S,1)**2 - twoMr(ixO^S) + eqpar(a_)**2


    u4Coord(ixO^S,0) = u4BL(ixO^S,0) + twoMr(ixO^S)/Delta(ixO^S) * u4BL(ixO^S,1)
    u4Coord(ixO^S,1) = u4BL(ixO^S,1)
    {^IFPHI
    u4Coord(ixO^S,^PHI) = u4BL(ixO^S,^PHI) + eqpar(a_)/Delta(ixO^S) * u4BL(ixO^S,1)
    }
    {^IFZ
    u4Coord(ixO^S,^Z) = u4BL(ixO^S,^Z)
    }

  end subroutine u4BLtoCoord
  !=============================================================================
  subroutine KSToCOORD(ixI^L,ixO^L,xKS,xCoord)

    ! Identity
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    !-----------------------------------------------------------------------------

    xCoord(ixO^S,:) = xKS(ixO^S,:)

  end subroutine KSToCoord
  !=============================================================================
  subroutine CoordToKS(ixI^L,ixO^L,xCoord,xKS)

    ! Identity
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xKS
    ! .. local ..
    !-----------------------------------------------------------------------------

    xKS(ixO^S,:) = xCoord(ixO^S,:)

  end subroutine CoordToKS
  !=============================================================================
  subroutine u4KStoCoord(ixI^L,ixO^L,xKS,u4KS,u4Coord)

    ! Identity
    !
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4KS
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    !-----------------------------------------------------------------------------

    u4Coord(ixO^S,:) = u4KS(ixO^S,:)

  end subroutine u4KStoCoord
  !=============================================================================
  subroutine KS_u4CoordToKS(ixI^L,ixO^L,x,u4Coord,u4KS,J)

    ! Identity
    !
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4Coord
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4KS
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    !-----------------------------------------------------------------------------

    if (present(J)) then
       call mpistop('KS_u4CoordToKS: Jacobian (triv.) not implemented.')
    else
       u4KS(ixO^S,:) = u4Coord(ixO^S,:)
    end if

  end subroutine KS_u4CoordToKS
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
    !-----------------------------------------------------------------------------

    call BLToKS(ixI^L,ixO^L,xBL,xCoord)

  end subroutine BLToCoord
  !=============================================================================
  subroutine CoordToBL(ixI^L,ixO^L,xCoord,xBL)

    ! Transforms the (spatial) KS coordinate to Boyer-Lindquist coordinates.
    ! We don't return the new time-component although it changes.
    !
    use mod_transform, only: KSToBL
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xBL
    ! .. local ..
    !-----------------------------------------------------------------------------

    call KSToBL(ixI^L,ixO^L,xCoord,xBL)

  end subroutine CoordToBL
  !=============================================================================
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)

    ! Transforms any contravariant three vector u3Coord in the current coordinates
    ! to a vector in the basis of "Kerr-Schild Cartesian coordinates".
    ! x = ra Sin(thetaKS) Cos(phiKS)
    ! y = ra Sin(thetaKS) Sin(phiKS)
    ! z = rKS Cos(thetaKS)
    ! where ra=sqrt(rKS**2+a**2) ! WARNING: I decided to put ra=rKS for now...
    !
    ! We are already in KS-cordinates, so we only have to do one transformation.  
    ! The transformation matrix J=del(x,y,z)/del(rKS,thetaKS,phiKS)
    ! Since the coordinates don't depend on time, this is the same transformation as for the
    ! four-vector u4CoordToCart()
    !
    !     | rKS/ra Sin(thetaKS) Cos(phiBl)    ra Cos(thetaKS) Cos(phiKS)    - ra Sin(thetaKS) Sin(phiKS)  |
    ! J = | rKS/ra Sin(thetaKS) Sin(phiKS)    ra Cos(thetaKS) Sin(phiKS)      ra Sin(thetaKS) Cos(phiKS)  |
    !     | Cos(thetaKS)                    - rKS Sin(thetaKS)                 0                          |
    !
    ! where ra=sqrt(rKS**2+a**2) ! WARNING: I decided to put ra=rKS for now...
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
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, ra
    integer                                                :: zcart_, ycart_, ix,jx
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


    {^IFPHIIN
    sph(ixO^S) = sin(x(ixO^S,phi_))
    cph(ixO^S) = cos(x(ixO^S,phi_))
    }{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    }

    {^IFZIN
    sth(ixO^S) = sin(x(ixO^S,z_))
    cth(ixO^S) = cos(x(ixO^S,z_))
    }{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    }

    ra(ixO^S) = x(ixO^S,1)!sqrt(x(ixO^S,1)**2+eqpar(a_)**2)

    Jac(ixO^S,1,1) = x(ixO^S,1)/ra(ixO^S) * sth(ixO^S)*cph(ixO^S)
    {^IFZ
    Jac(ixO^S,1,^Z) = ra(ixO^S) * cth(ixO^S)*cph(ixO^S)
    }
    {^IFPHI
    Jac(ixO^S,1,^PHI) = - ra(ixO^S)*sth(ixO^S)*sph(ixO^S)
    }

    if (ycart_.ne.1) then
       Jac(ixO^S,ycart_,1) = x(ixO^S,1)/ra(ixO^S) * sth(ixO^S)*sph(ixO^S)
       {^IFZ
       Jac(ixO^S,ycart_,^ZZ) = ra(ixO^S) * cth(ixO^S)*sph(ixO^S)
       }
       {^IFPHI
       Jac(ixO^S,ycart_,^PPHI) = ra(ixO^S)*sth(ixO^S)*cph(ixO^S)
       }
    end if

    if(zcart_.ne.1) then
       Jac(ixO^S,zcart_,1) = cth(ixO^S)
       {^IFZ
       Jac(ixO^S,zcart_,^ZZ) = - x(ixO^S,1)*sth(ixO^S)
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
  subroutine CoordToCart(ixI^L,ixO^L,x,xCart)

    ! Transforms the coordinates to "Kerr-Schild Cartesian coordinates"
    ! x = ra Sin(thetaKS) Cos(phiKS)
    ! y = ra Sin(thetaKS) Sin(phiKS)
    ! z = rKS Cos(thetaKS)
    ! where ra=sqrt(rKS**2+a**2) ! WARNING: I decided to put ra=rKS for now...
    ! We are already in KS-cordinates, so we only have to do one transformation.
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCart
    ! .. local ..
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, ra
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

    {^IFPHIIN
    sph(ixO^S) = sin(x(ixO^S,phi_))
    cph(ixO^S) = cos(x(ixO^S,phi_))
    }{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    }

    {^IFZIN
    sth(ixO^S) = sin(x(ixO^S,z_))
    cth(ixO^S) = cos(x(ixO^S,z_))
    }{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    }

    ra(ixO^S) = x(ixO^S,1)!sqrt(x(ixO^S,1)**2+eqpar(a_)**2)


    xCart(ixO^S,1) = ra(ixO^S) * sth(ixO^S)*cph(ixO^S)
    if (ycart_.ne.1) &
         xCart(ixO^S,ycart_) = ra(ixO^S) * sth(ixO^S)*sph(ixO^S)
    if (zcart_.ne.1) &
         xCart(ixO^S,zcart_) = x(ixO^S,1)*cth(ixO^S)

  end subroutine CoordToCart
  !=============================================================================
  subroutine KS_LRNFu(ixI^L,ixO^L,m,ehatu)

    ! Obtains orthonormal tetrad vectors in matrix form.
    ! Follows Rohta Takahashi, MNRAS, 383, 1155-116, 2008
    ! His equations (11) - (14)
    ! Checked with Mathematica.
    ! This takes in already the KS metric datastructure 'm', I might regret
    ! this choice later on. 
    
    include 'amrvacdef.f'
    integer, intent(in)                                         :: ixI^L, ixO^L
    type(metric)                                                :: m
    double precision, dimension(ixI^S,0:^NC,0:^NC), intent(out) :: ehatu
    ! .. local ..
    !-----------------------------------------------------------------------------

    ehatu = zero
    
    ehatu(ixO^S,0,0) = m%alpha(ixO^S)

    ehatu(ixO^S,r_,0)  = m%beta(r_)%elem(ixO^S)/sqrt(m%gammainv(r_,r_)%elem(ixO^S))
    ehatu(ixO^S,r_,r_) = one/sqrt(m%gammainv(r_,r_)%elem(ixO^S))

    {^IFZ
    ehatu(ixO^S,z_,z_) = sqrt(m%g(z_,z_)%elem(ixO^S))
    }

    {^IFPHI
    ehatu(ixO^S,phi_,0)    = m%beta(r_)%elem(ixO^S) * m%g(r_,phi_)%elem(ixO^S) / sqrt(m%g(phi_,phi_)%elem(ixO^S))
    ehatu(ixO^S,phi_,r_)   = m%g(r_,phi_)%elem(ixO^S) / sqrt(m%g(phi_,phi_)%elem(ixO^S))
    ehatu(ixO^S,phi_,phi_) = sqrt(m%g(phi_,phi_)%elem(ixO^S))
    }
    
  end subroutine KS_LRNFu
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
