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
  ! Metric components for spherical coordinates -coord=splog
  ! This file is included by the preprocessor in mod_metric.t
  ! It will be part of the module mod_metric.
  !  
  ! Pushpita Das
  ! 23.01.2020
  !=============================================================================
  
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------
  
  character*20, parameter              :: coord="splog"
  logical, save                        :: init_from_g4 = .false.  ! We don't provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: ncoordpar = 0
  double precision, save               :: coordpar(ncoordpar) 

contains
  !=============================================================================
  subroutine init_coord

  end subroutine init_coord
  !=============================================================================
  double precision function s(rsp)
    include 'amrvacdef.f'

    double precision :: rsp

    s = log(rsp)

  end function s
  !================================================================================
  double precision function rsp(s)

    double precision :: s

    rsp = exp(s)

  end function rsp
  !===============================================================================
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
    sqrtgamma = (exp(x1)) * (rsp(x1))**2 * abs(sin(x^Z))
    }{^IFZOUT
    sqrtgamma = (exp(x1)) * (rsp(x1))**2
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

    alpha = 1.0d0

    if (present(jdir)) then
       if (jdir .eq. 1) then
          if (present(dalphadj)) then
             dalphadj = 0.0d0
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .true.
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
          g = exp(2*x1)
          if (present(kdir)) then
             if (kdir .eq.1) then
                if (present(dgdk)) dgdk = 2.0d0*rsp(x1)*exp(x1)
                if(present(dgdk_iszero)) dgdk_iszero = .false.
             else
                if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.
             end if
          end if

       else if (i .eq. ^PHI) then ! phi direction
          {^IFZIN
          g = (rsp(x1)*sin(x^Z))**2
          }{^IFZOUT
          ! equatorial region
          g = (rsp(x1))**2      
          }

          if (present(kdir)) then
             if (kdir .eq.1) then
                {^IFZIN
                if (present(dgdk)) dgdk = 2.0d0*rsp(x1)*exp(x1)*(sin(x^Z))**2
                if (present(dgdk_iszero)) dgdk_iszero = .false.
                }{^IFZOUT
                if (present(dgdk)) dgdk = 2.0d0*rsp(x1)*exp(x1)
                if (present(dgdk_iszero)) dgdk_iszero = .false.
                }
             else if (kdir .eq. ^Z) then
                {^IFZIN
                if (present(dgdk)) dgdk = (rsp(x1))**2*(2.0d0*cos(x^Z)*sin(x^Z))
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
          g = (rsp(x1))**2

          if (present(kdir)) then
             if (kdir .eq.1) then
                if (present(dgdk)) dgdk = 2.0d0*rsp(x1)*exp(x1)
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

    call mpistop('outerhorizon: not possible for spherical coordinates')

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
    integer                                                :: ix^D
    !-----------------------------------------------------------------------------
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    xCoord(ix^D,1) = s(xBL(ix^D,1))
    {^IFZIN
    xCoord(ix^D,^Z) = xBL(ix^D,^Z)
    }
    {^IFPHIIN
    xCoord(ix^D,^PHI) = xBL(ix^D,^PHI)
    }
    {end do\}
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
    integer                                                :: ix^D
    !-----------------------------------------------------------------------------
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    xBL(ix^D,1) = rsp(xCoord(ix^D,1))
    {^IFZIN
    xBL(ix^D,^Z) = xCoord(ix^D,^Z)
    }
    {^IFPHIIN
    xBL(ix^D,^PHI) = xCoord(ix^D,^PHI)
    }
    {end do\}


  end subroutine CoordToBL
  !=============================================================================
  subroutine u4BLtoCoord(ixI^L,ixO^L,xBL,u4BL,u4Coord)

    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to the current coordinates u4Coord.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    !
    include 'amrvacdef.f'

    integer                                                :: ixI^L, ixO^L
    double precision, dimension(ixI^S,0:^ND),intent(in)    :: xBL
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    !double precision,dimension(ixI^S)                      :: ems
    !-----------------------------------------------------------------------------
    if (eqpar(a_) .ne. zero) &
         call mpistop("u4BLtoCoord: Spin not zero, transformation to Schwarzschild makes no sense.")


    u4Coord(ixO^S,0) = u4BL(ixO^S,0)
    u4Coord(ixO^S,r_) = u4BL(ixO^S,r_)/xBL(ixO^S,1)
    {^IFZIN
    u4Coord(ixO^S,z_) = u4BL(ixO^S,z_)
    }
    {^IFPHI
    u4Coord(ixO^S,phi_) = u4BL(ixO^S,phi_)
    }
  end subroutine u4BLtoCoord
  !============================================================================= 
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)

    ! Transforms any contravariant three vector u3Coord in the current coordinates
    ! to a vector in the basis of "Boyer-Lindquist Cartesian coordinates".
    ! rsp = exp(s)
    ! x = rsp Sin(thetaBL) Cos(phiBL)
    ! y = rsp Sin(thetaBL) Sin(phiBL)
    ! z = rsp Cos(thetaBL)  
    ! The transformation matrix J=del(x,y,z)/del(s,thetaBL,phiBL)
    ! Since the coordinates don't depend on time, this is the same transformation as for the
    ! four-vector u4CoordToCart()
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
    double precision, dimension(ixI^S,1:^ND)               :: xBL
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

    call CoordToBL(ixI^L,ixO^L,x,xBL)

    {^IFPHIIN
    sph(ixO^S) = sin(xBL(ixO^S,phi_))
    cph(ixO^S) = cos(xBL(ixO^S,phi_))
    }{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    }

    {^IFZIN
    sth(ixO^S) = sin(xBL(ixO^S,z_))
    cth(ixO^S) = cos(xBL(ixO^S,z_))
    }{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    }


    Jac(ixO^S,1,1) = exp(x(ixO^S,1)) * sth(ixO^S)*cph(ixO^S)
    {^IFZ
    Jac(ixO^S,1,^Z) = xBL(ixO^S,1) * cth(ixO^S)*cph(ixO^S)
    }
    {^IFPHI
    Jac(ixO^S,1,^PHI) = - xBL(ixO^S,1)*sth(ixO^S)*sph(ixO^S)
    }

    if (ycart_.ne.1) then
       Jac(ixO^S,ycart_,1) = exp(x(ixO^S,1)) * sth(ixO^S)*sph(ixO^S)
       {^IFZ
       Jac(ixO^S,ycart_,^ZZ) = xBL(ixO^S,1) * cth(ixO^S)*sph(ixO^S)
       }
       {^IFPHI
       Jac(ixO^S,ycart_,^PPHI) = xBL(ixO^S,1)*sth(ixO^S)*cph(ixO^S)
       }
    end if

    if(zcart_.ne.1) then
       Jac(ixO^S,zcart_,1) = exp(x(ixO^S,1)) * cth(ixO^S)
       {^IFZ
       Jac(ixO^S,zcart_,^ZZ) = - xBL(ixO^S,1)*sth(ixO^S)
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
    ! rsp = exp(s)
    !
    ! x = rsp Sin(thetaBL) Cos(phiBL)
    ! y = rsp Sin(thetaBL) Sin(phiBL)
    ! z = rsp Cos(thetaBL)
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCart
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xBL
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

    call CoordToBL(ixI^L,ixO^L,x,xBL) 

    {^IFPHIIN
    sph(ixO^S) = sin(xBL(ixO^S,phi_))
    cph(ixO^S) = cos(xBL(ixO^S,phi_))
    }{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    }

    {^IFZIN
    sth(ixO^S) = sin(xBL(ixO^S,z_))
    cth(ixO^S) = cos(xBL(ixO^S,z_))
    }{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    }

    ra(ixO^S) = xBL(ixO^S,1)!sqrt(x(ixO^S,1)**2+eqpar(a_)**2)


    xCart(ixO^S,1) = ra(ixO^S) * sth(ixO^S)*cph(ixO^S)
    if (ycart_.ne.1) &
         xCart(ixO^S,ycart_) = ra(ixO^S) * sth(ixO^S)*sph(ixO^S)
    if (zcart_.ne.1) &
         xCart(ixO^S,zcart_) = ra(ixO^S)*cth(ixO^S)

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
