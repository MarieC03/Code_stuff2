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
  ! Metric components for Schwarzschild coordinates -coord=ss
  ! This file is included by the preprocessor in mod_metric.t
  ! It will be part of the module mod_metric.
  !  
  ! Oliver Porth
  ! 04.12.2015
  !=============================================================================
  
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------
  
  character*20, parameter              :: coord="ss"
  logical, save                        :: init_from_g4 = .false.  ! We don't provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: ncoordpar = 0
  double precision, save               :: coordpar(ncoordpar) 

contains
!=============================================================================
  subroutine init_coord

   get_gammainv_component_analytic => ss_get_gammainv_component_analytic
   get_gammainv_component => get_gammainv_component_analytic
    
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
    sqrtgamma = (x1**2 * abs(sin(x^Z)) / sqrt(1.0d0 - 2.0d0*eqpar(m_)/x1))
    }{^IFZOUT
    sqrtgamma = (x1**2 / sqrt(1.0d0 - 2.0d0*eqpar(m_)/x1))
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

    alpha = sqrt(1.0d0-2.0d0*eqpar(m_)/x1)

    if (present(jdir)) then
       if (jdir .eq. 1) then
          if (present(dalphadj)) then
             dalphadj = eqpar(m_)/alpha/x1**2
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

    beta = 0.0d0
    if(present(iszero)) &
         iszero = .true.

    ! \partial_j \beta^i
    if (present(dbetaidj)) &
         dbetaidj = 0.0d0
    if (present(dbetaidj_iszero)) &
         dbetaidj_iszero = .true.

  end subroutine get_beta
  !=============================================================================
  subroutine ss_get_gammainv_component_analytic(iin,jin,x^D,ginv,iszero,dginvdk_iszero,dginvdk,kdir)

    include 'amrvacdef.f'

    ! Indices of the 3-metric are up (contravariant) gammainv^{ij}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dginvdk and kdir request derivatives of the metric
    ! \partial_k gammainv^{ij} ; i=iin, j=jin, k=kdir
    ! The optional argument dginvdk_iszero

    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: ginv
    logical, optional, intent(out)           :: iszero, dginvdk_iszero
    double precision, optional, intent(out)  :: dginvdk
    integer                                  :: i,j
    !-----------------------------------------------------------------------------

    if(present(dginvdk) .and. .not. present(kdir) .or. &
         present(dginvdk_iszero) .and. .not. present(kdir)) &
         call mpistop("get_gammainv_component: derivatives requested without &
         &direction or output-slot given.")

    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    ! Diagonal for now
    if (i==j) then
       if(present(iszero)) &
            iszero = .false.

       if (i .eq. 1) then ! r direction
          ginv = (1.0d0-2.0d0*eqpar(m_)/x1)
          if (present(kdir)) then
             if (kdir .eq.1) then
                if (present(dginvdk)) &
                     dginvdk = 2.0d0*eqpar(m_)/x1**2
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             else
                if (present(dginvdk)) dginvdk = 0.0d0
                if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             end if
          end if

       else if (i .eq. ^PHI) then ! phi direction
          {^IFZIN
          ginv = 1.0d0/(x1*sin(x^Z))**2
          }{^IFZOUT
          ! equatorial region
          ginv = 1.0d0/x1**2      
          }

          if (present(kdir)) then
             if (kdir .eq.1) then
                {^IFZIN
                if (present(dginvdk)) dginvdk = - 2.0d0/(x1**3*sin(x^Z)**2)
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
                }{^IFZOUT
                if (present(dginvdk)) dginvdk = - 2.0d0/x1**3
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
                }
             else if (kdir .eq. ^Z) then
                {^IFZIN
                if (present(dginvdk)) dginvdk = - 2.0d0*cos(x^Z)/(x1**2*sin(x^Z)**3)
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
                }{^IFZOUT
                if (present(dginvdk)) dginvdk = 0.0d0
                if (present(dginvdk_iszero)) dginvdk_iszero = .true.
                }
             else
                if (present(dginvdk)) dginvdk = 0.0d0
                if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             end if
          end if

       else if (i .eq. ^Z) then ! theta direction
          ginv = 1.d0/x1**2

          if (present(kdir)) then
             if (kdir .eq.1) then
                if (present(dginvdk)) dginvdk = - 2.0d0/x1**3
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             else
                if (present(dginvdk)) dginvdk = 0.0d0
                if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             end if
          end if

       end if
    else
       ginv = 0.0d0

       if(present(iszero)) &
            iszero = .true.

       if (present(dginvdk)) &
            dginvdk = 0.0d0

       if(present(dginvdk_iszero)) &
            dginvdk_iszero = .true.

    end if
    
  end subroutine ss_get_gammainv_component_analytic
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
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    ! Diagonal for now
    if (i==j) then
       if(present(iszero)) &
            iszero = .false.

       if (i .eq. 1) then ! r direction
          g = 1.0d0/(1.0d0-2.0d0*eqpar(m_)/x1)
          if (present(kdir)) then
             if (kdir .eq.1) then
                if (present(dgdk)) &
                     dgdk = - 2.0d0*eqpar(m_)*(g/x1)**2
                if(present(dgdk_iszero)) dgdk_iszero = .false.
             else
                if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.
             end if
          end if

       else if (i .eq. ^PHI) then ! phi direction
          {^IFZIN
          g = (x1*sin(x^Z))**2
          }{^IFZOUT
          ! equatorial region
          g = x1**2      
          }

          if (present(kdir)) then
             if (kdir .eq.1) then
                {^IFZIN
                if (present(dgdk)) dgdk = 2.0d0*x1*sin(x^Z)**2
                if (present(dgdk_iszero)) dgdk_iszero = .false.
                }{^IFZOUT
                if (present(dgdk)) dgdk = 2.0d0*x1
                if (present(dgdk_iszero)) dgdk_iszero = .false.
                }
             else if (kdir .eq. ^Z) then
                {^IFZIN
                if (present(dgdk)) dgdk = 2.0d0*x1**2*cos(x^Z)*sin(x^Z)
                if (present(dgdk_iszero)) dgdk_iszero = .false.
                }{^IFZOUT
                if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.
                }
             else
                if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.
             end if
          end if

       else if (i .eq. ^Z) then ! theta direction
          g = x1**2

          if (present(kdir)) then
             if (kdir .eq.1) then
                if (present(dgdk)) dgdk = 2.0d0*x1
                if (present(dgdk_iszero)) dgdk_iszero = .false.
             else
                if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.
             end if
          end if

       end if
    else
       g = 0.0d0

       if(present(iszero)) &
            iszero = .true.

       if (present(dgdk)) &
            dgdk = 0.0d0

       if(present(dgdk_iszero)) &
            dgdk_iszero = .true.

    end if

  end subroutine get_g_component
  !=============================================================================
    double precision function outerhorizon()
    
    include 'amrvacdef.f'
    !-----------------------------------------------------------------------------

    outerhorizon = 2.0d0 * eqpar(m_)
    
  end function outerhorizon
  !=============================================================================
  subroutine BLToCoord(ixI^L,ixO^L,xBL,xCoord)

    ! identity
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    !-----------------------------------------------------------------------------

    xCoord = xBL

  end subroutine BLToCoord
  !=============================================================================
  subroutine CoordToBL(ixI^L,ixO^L,xCoord,xBL)

    ! identity
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xBL
    ! .. local ..
    !-----------------------------------------------------------------------------

    xBL = xCoord

  end subroutine CoordToBL
  !=============================================================================
  subroutine u4BLtoCoord(ixI^L,ixO^L,x,u4BL,u4Coord)

    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to the current (SS) coordinates u4Coord.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    !
    ! The transformation is trivial since SS coordinates are just BL for a=0
    ! If a ne zero this makes no sense.
    !
    include 'amrvacdef.f'

    integer                                                :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    !-----------------------------------------------------------------------------
    if (eqpar(a_) .ne. zero) &
         call mpistop("u4BLtoCoord: Spin not zero, transformation to Schwarzschild makes no sense.")

    u4Coord(ixO^S,0:^NC) = u4BL(ixO^S,0:^NC)

  end subroutine u4BLtoCoord
  !============================================================================= 
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)

    ! Transforms any contravariant three vector u3Coord in the current coordinates
    ! to a vector in the basis of "Boyer-Lindquist Cartesian coordinates".
    ! x = sqrt(rBL**2+a**2) Sin(thetaBL) Cos(phiBL)
    ! y = sqrt(rBL**2+a**2) Sin(thetaBL) Sin(phiBL)
    ! z = rBL Cos(thetaBL)
    ! We are already in BL-cordinates, so we only have to do one transformation.  
    ! The transformation matrix J=del(x,y,z)/del(rBL,thetaBL,phiBL)
    ! Since the coordinates don't depend on time, this is the same transformation as for the
    ! four-vector u4CoordToCart()
    !
    !     | rBL/ra Sin(thetaBL) Cos(phiBl)    ra Cos(thetaBL) Cos(phiBL)    - ra Sin(thetaBL) Sin(phiBL)  |
    ! J = | rBL/ra Sin(thetaBL) Sin(phiBL)    ra Cos(thetaBL) Sin(phiBL)      ra Sin(thetaBL) Cos(phiBL)  |
    !     | Cos(thetaBL)                    - rBL Sin(thetaBL)                 0                          |
    !
    ! where ra=sqrt(rBL**2+a**2) ! WARNING: I decided to put ra=rBL for now...
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

    ! Transforms the coordinates to "Boyer-Lindquist Cartesian coordinates"
    ! x = ra Sin(thetaBL) Cos(phiBL)
    ! y = ra Sin(thetaBL) Sin(phiBL)
    ! z = rBL Cos(thetaBL)
    ! where ra=sqrt(rBL**2+a**2) ! WARNING: I decided to put ra=rBL for now...
    ! We are already in BL-cordinates, so we only have to do one transformation.
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
