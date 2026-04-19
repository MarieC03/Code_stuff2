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
  ! Metric components for braneworld black hole in Randall-Sundrum model
  ! This file is included by the preprocessor in mod_metric.t
  ! It will be part of the module mod_metric.
  !  
  ! Hector Olivares and Bart Ripperda
  ! 25.01.2019
  !=============================================================================
  
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------
  
  character*20, parameter              :: coord="brane"
  logical, save                        :: init_from_g4 = .false.  ! We don't provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: R0_=1,h_=2,q_=3,Lambda_=4
  integer, parameter                   :: ncoordpar = 4
  double precision, save               :: coordpar(ncoordpar) 

contains
!=============================================================================
  subroutine init_coord

    SPToCoord=>BLToCoord
    
  end subroutine init_coord
  !=============================================================================
  double precision function rH()
    ! Calculates the radial position of the outer event horizon in BL-type
    ! coordinates
    include 'amrvacdef.f'
    ! .. local ..
    !-----------------------------------------------------------------------------

    rH = eqpar(m_)+sqrt(eqpar(m_)**2 - eqpar(a_)**2 - coordpar(q_)**2)

  end function rH 
  !=============================================================================
  double precision function s(rks)

    include 'amrvacdef.f'
    double precision                               :: rks
    ! .. local ..
    !-----------------------------------------------------------------------------

    s = log( rks - coordpar(R0_) )

  end function s
  !=============================================================================
  double precision function vartheta(theta)

    use mod_rtsafe2D, only: MAXIT, rt1D
    include 'amrvacdef.f'
    double precision                               :: theta
    ! .. local ..
    double precision                               :: rtmp
    integer                                        :: ierror
    !-----------------------------------------------------------------------------

    vartheta = theta
    if (coordpar(h_) .eq. zero) then 
       return
    else

       MAXIT = 30
       call rt1D(func_vartheta,vartheta,absaccnr,ierror)

       if (ierror .ne. 0) call mpistop('vartheta: Does not converge')

    end if
  contains
    !=============================================================================
    subroutine func_vartheta(myvartheta,f,df)

      double precision, intent(in)      :: myvartheta
      double precision, intent(out)     :: f, df
      ! .. local ..
      !-----------------------------------------------------------------------------

      f  = myvartheta+(coordpar(h_)*Sin(2.0d0*myvartheta))/2.0d0 - theta
      df = 1.0d0 + coordpar(h_) * Cos(2.0d0*myvartheta)

    end subroutine func_vartheta
    !=============================================================================
  end function vartheta
  !=============================================================================
  double precision function rks(s)

    double precision                               :: s
    ! .. local ..
    !-----------------------------------------------------------------------------

    rks = coordpar(R0_) + exp(s)

  end function rks
  !=============================================================================
  double precision function thetaks(vartheta)

    double precision                               :: vartheta

    ! .. local ..
    !-----------------------------------------------------------------------------

    thetaks = vartheta+0.5d0*coordpar(h_)*Sin(2.0d0*vartheta)

  end function thetaks
  !=============================================================================
  subroutine get_coords(x^D,r,th,phi)

    include 'amrvacdef.f'
    double complex, intent(in)  :: x^D
    double complex              :: r, th, phi
    ! .. local ..
  !-----------------------------------------------------------------------------

    r =  coordpar(R0_) + exp(x1)
    th = 0.5d0*dpi
    phi = 0.0d0

    {^IFZIN
    th = x^Z + 0.5d0*coordpar(h_)*sin(2.0d0*x^Z)

    }
    {^PHIIN
    phi = x^PHI

    }

  end subroutine get_coords
  !=============================================================================
  subroutine gethelpers(r,th,phi,Sigma,Deltar,Deltath,Zeta)
    include 'amrvacdef.f'

    double complex, intent(in)          :: r, th, phi
    double complex, intent(out)         :: Sigma, Deltar, Deltath, Zeta
    ! .. local ..

  !-----------------------------------------------------------------------------

    Sigma = r**2 + (eqpar(a_)*cos(th))**2


!    Deltar = (r**2 + eqpar(a_)**2) * (one - coordpar(Lambda_)*r**2/3.d0) &
!              - two*eqpar(m_)*r + coordpar(q_)

    ! !!! Not standard notation... !!!
    Deltar = two*eqpar(m_)*r - coordpar(q_)**2

    Deltath = one + coordpar(Lambda_)*(eqpar(a_)*cos(th))**2/3.d0

    Zeta = one + coordpar(Lambda_)*eqpar(a_)**2/3.d0


  end subroutine gethelpers
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


    sqrtgamma = one

    if(present(is_analytic)) is_analytic = .false.

  end subroutine get_sqrtgamma_analytic
  !=============================================================================
  double complex function cplx_alpha(x^D) 
    include 'amrvacdef.f'
    double complex, intent(in)  :: x^D
    ! .. local ..
    double complex              :: r, th, phi, Sigma, Deltar, Deltath, Zeta
    !-----------------------------------------------------------------------------

    call get_coords(x^D,r,th,phi)

    call gethelpers(r,th,phi,Sigma,Deltar,Deltath,Zeta)

    cplx_alpha = sqrt(Sigma/(Sigma + Deltar))

  end function cplx_alpha
  !=============================================================================
  double complex function cplx_betar(x^D) 
    include 'amrvacdef.f'
    double complex, intent(in)  :: x^D
    ! .. local ..
    double complex              :: r, th, phi, Sigma, Deltar, Deltath, Zeta
    !-----------------------------------------------------------------------------

    call get_coords(x^D,r,th,phi)

    call gethelpers(r,th,phi,Sigma,Deltar,Deltath,Zeta)

    cplx_betar = Deltar/(Deltar + Sigma)

    ! Coordinate transformation to logarithmic

    cplx_betar = cplx_betar/r

  end function cplx_betar
  !=============================================================================
  double complex function cplx_grr(x^D) 
    include 'amrvacdef.f'
    double complex, intent(in)  :: x^D
    ! .. local ..
    double complex              :: r, th, phi, Sigma, Deltar, Deltath, Zeta
    !-----------------------------------------------------------------------------

    call get_coords(x^D,r,th,phi)

    call gethelpers(r,th,phi,Sigma,Deltar,Deltath,Zeta)

    cplx_grr = 1.0d0 + Deltar/Sigma

    ! Coordinate transformation to logarithmic

    cplx_grr = cplx_grr*r**2

  end function cplx_grr
  !=============================================================================
  double complex function cplx_gphph(x^D) 
    include 'amrvacdef.f'
    double complex, intent(in)  :: x^D
    ! .. local ..
    double complex              :: r, th, phi, Sigma, Deltar, Deltath, Zeta
    !-----------------------------------------------------------------------------

    call get_coords(x^D,r,th,phi)

    call gethelpers(r,th,phi,Sigma,Deltar,Deltath,Zeta)

    cplx_gphph = (Deltar*(eqpar(a_)*sin(th))**2/Sigma + r**2 + eqpar(a_)**2)*sin(th)**2

  end function cplx_gphph
  !=============================================================================
  double complex function cplx_gthth(x^D) 
    include 'amrvacdef.f'
    double complex, intent(in)  :: x^D
    ! .. local ..
    double complex              :: r, th, phi, Sigma, Deltar, Deltath, Zeta
    !-----------------------------------------------------------------------------

    call get_coords(x^D,r,th,phi)

    call gethelpers(r,th,phi,Sigma,Deltar,Deltath,Zeta)

    cplx_gthth = Sigma

    ! Coordinate transformation to logarithmic

{^IFZIN
    cplx_gthth = cplx_gthth*(1.0d0+coordpar(h_)*cos(2.0d0*x^Z))**2
}

  end function cplx_gthth
  !=============================================================================
  double complex function cplx_grph(x^D) 
    include 'amrvacdef.f'
    double complex, intent(in)  :: x^D
    ! .. local ..
    double complex              :: r, th, phi, Sigma, Deltar, Deltath, Zeta
    !-----------------------------------------------------------------------------

    call get_coords(x^D,r,th,phi)

    call gethelpers(r,th,phi,Sigma,Deltar,Deltath,Zeta)

    cplx_grph = -(1.0d0 + Deltar/Sigma)*eqpar(a_)*(sin(th))**2

    ! Coordinate transformation to logarithmic

    cplx_grph = cplx_grph*r
 
  end function cplx_grph
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


    alpha = cplx_alpha({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})

    if(present(iszero)) iszero = .false.

    if (present(dalphadj_iszero)) dalphadj_iszero = .false.

    if (present(dalphadj)) call complex_derivative(x^D,cplx_alpha,jdir,dalphadj)

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

    beta = zero

    if(present(iszero)) iszero = .true.
    ! \partial_j \beta^i
    if (present(dbetaidj)) dbetaidj = 0.0d0
    if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.

    if (idir.eq.1) then

       beta = cplx_betar({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
       if(present(iszero)) iszero = .false.
       if (present(dbetaidj)) call complex_derivative(x^D,cplx_betar,jdir,dbetaidj)
       if (present(dbetaidj_iszero)) dbetaidj_iszero = .false.

    end if

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

    ! Initialize everything to zero to fill later only the non-zero components
    if(present(iszero)) iszero = .true.
    g = zero
    if (present(dgdk)) dgdk = zero
    if(present(dgdk_iszero)) dgdk_iszero = .true.

    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    ! Diagonal unity

    select case(i)
    case(r_)
       select case(j)
       case(r_)
          g = cplx_grr({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if(present(iszero)) iszero = .false.
          if (present(dgdk_iszero)) dgdk_iszero = .false.
          if (present(dgdk)) call complex_derivative(x^D,cplx_grr ,kdir,dgdk)

       case(phi_)
          g = cplx_grph({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if(present(iszero)) iszero = .false.
          if (present(dgdk_iszero)) dgdk_iszero = .false.
          if (present(dgdk)) call complex_derivative(x^D,cplx_grph ,kdir,dgdk)


       end select
    case(z_)
       if (j .eq. z_) then
          g = cplx_gthth({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if(present(iszero)) iszero = .false.
          if (present(dgdk_iszero)) dgdk_iszero = .false.
          if (present(dgdk)) call complex_derivative(x^D,cplx_gthth ,kdir,dgdk)
       end if

    case(phi_)
       if (j .eq. phi_) then
          g = cplx_gphph({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if(present(iszero)) iszero = .false.
          if (present(dgdk_iszero)) dgdk_iszero = .false.
          if (present(dgdk)) call complex_derivative(x^D,cplx_gphph ,kdir,dgdk)
       end if

    end select

  end subroutine get_g_component
  !=============================================================================
    double precision function outerhorizon()
    
    !-----------------------------------------------------------------------------

    print *, 'Function other horizon still not implemented'
    outerhorizon = 0.0d0
    
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
    ! uthetaKS   = uthetaBL / (1+h*Cos(2.0d0*vartheta))
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

    ems(ixO^S)   = one/(xBL(ixO^S,1) - coordpar(R0_))

    {^IFZIN
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    varth(ix^D) = vartheta(xBL(ix^D,^Z))
    {end do \}
    }{^IFZOUT
    varth(ixO^S) = dpi/2.0d0
    }

    twoMr(ixO^S) = 2.0d0*eqpar(m_)*abs(xBL(ixO^S,1)) - eqpar(q_)**2
    Delta(ixO^S) = xBL(ixO^S,1)**2 - twoMr(ixO^S) + eqpar(a_)**2

    u4Coord(ixO^S,0) = u4BL(ixO^S,0) + twoMr(ixO^S)/Delta(ixO^S) * u4BL(ixO^S,1)
    u4Coord(ixO^S,1) = ems(ixO^S) * u4BL(ixO^S,1)

    {^IFZ
    u4Coord(ixO^S,^Z) = 1/(1+coordpar(h_)*Cos(2.0d0*varth(ixO^S))) * u4BL(ixO^S,^Z)
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
    xCoord(ix^D,1) = s(xKS(ix^D,1))
    {^IFZIN
    xCoord(ix^D,^Z) = vartheta(xKS(ix^D,^Z))
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
    xKS(ix^D,1) = rks(xCoord(ix^D,1))
    {^IFZIN
    xKS(ix^D,^Z) = thetaks(xCoord(ix^D,^Z))
    }
    {^IFPHIIN
    xKS(ix^D,^PHI) = xCoord(ix^D,^PHI)
    }
    {end do\}

  end subroutine CoordToKS
  !=============================================================================
  subroutine MKS_u4CoordToKS(ixI^L,ixO^L,xCoord,u4Coord,u4KS,J)

    ! Transforms the modified KS four vector u4Coord to
    ! the ordinary KS four vector u4KS.
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4Coord
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4KS
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    integer                                                :: ix^D
    !-----------------------------------------------------------------------------

    u4KS(ixO^S,0)    = u4Coord(ixO^S,0)
    u4KS(ixO^S,r_)   = exp(xCoord(ixO^S,r_)) * u4Coord(ixO^S,r_)
    {^IFZIN
    u4KS(ixO^S,z_)   = (1.0d0 + coordpar(h_)*cos(2.0d0*xCoord(ixO^S,z_))) * u4Coord(ixO^S,z_)
    }
    {^IFPHI
    u4KS(ixO^S,phi_) = u4Coord(ixO^S,phi_)
    }

    if (present(J)) then
       J(ixO^S,:,:) = 0.0d0
       J(ixO^S,r_,r_) = exp(xCoord(ixO^S,r_))
       {^IFZIN
       J(ixO^S,z_,z_) = (1.0d0 + coordpar(h_)*cos(2.0d0*xCoord(ixO^S,z_)))
       }
    end if

  end subroutine MKS_u4CoordToKS
  !=============================================================================
  subroutine CoordToCart(ixI^L,ixO^L,x,xCart)

    ! Transforms the coordinates to "Kerr-Schild Cartesian coordinates"
    ! Two transformations are needed:  
    ! rKS     = R0 + exp(s)
    ! thetaKS = vartheta + h/2 Sin(2 vartheta)
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
    ! thetaKS = vartheta + h/2 Sin(2 vartheta)
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

    {^IFZIN
    tmp(ixO^S) = 1+coordpar(h_)*Cos(2.0d0*x(ixO^S,^Z))
    }{^IFZOUT
    tmp(ixO^S) = 1.0d0 - coordpar(h_)
    }

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
  subroutine MKS_LRNFu(ixI^L,ixO^L,m,ehatu)

    ! Obtains orthonormal tetrad vectors in matrix form.
    ! This is the locally non-rotating reference frame (LRNF)
    ! which to my understanding is just the Eulerian one, but its normalized.  
    ! Follows Rohta Takahashi, MNRAS, 383, 1155-116, 2008
    ! His equations (11) - (14)
    ! Should obey the same structure as KS case.
    ! This takes in already the MKS metric datastructure 'm', I might regret
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
    
  end subroutine MKS_LRNFu
  !=============================================================================
  subroutine KNKSToLRNF(ixI^L,ixO^L,xKS,ehatu)

    ! Hector Olivares 26.04.2019
    ! Based in KSToLRNF from Oliver Porth
    ! Obtains the tetrad transform to the locally non-rotating frame (LRNF)
    ! in the Kerr-Newman spacetime in Kerr-Schild coordinates
    ! Needs KS coordinates on input.
    ! Needs to be applied to tensor in (spherical) KS coordinates.
    ! Returns the transformation matrix e^\hat{\mu}_\nu
    ! 
    ! Follows Rohta Takahashi, MNRAS, 383, 1155-1165, 2008
    ! His equations (11) - (14), replacing 2Mr --> 2Mr - Q**2
    ! and Rosquist 2009 Gen Relativ Gravit, 41:2619¿2632, 2009
    ! DOI 10.1007/s10714-009-0789-7
    ! Checked with Sage.

    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS ! KS-coordinate
    double precision, dimension(ixI^S,0:^NC,0:^NC), intent(out) :: ehatu
    ! .. local ..
    double precision, dimension(ixI^S)        :: sth, cth, rho2, TwoMrQ
    !-----------------------------------------------------------------------------

    {^IFZIN
    sth(ixO^S) = sin(xKS(ixO^S,^Z))
    cth(ixO^S) = cos(xKS(ixO^S,^Z))
    }{^IFZOUT
    sth(ixO^S) = 1.0d0   ! xy-Plane
    cth(ixO^S) = 0.0d0   ! xy-Plane
    }
    rho2(ixO^S) = xKS(ixO^S,r_)**2 + eqpar(a_)**2 * cth(ixO^S)**2

    ehatu(ixO^S,:,:) = zero

    TwoMrQ(ixO^S) = 2.0d0*eqpar(m_)*xKS(ixO^S,r_) - coordpar(q_)**2

    ehatu(ixO^S,0,0) = 1/Sqrt(1+TwoMrQ(ixO^S)/rho2(ixO^S))

    ehatu(ixO^S,r_,0) = TwoMrQ(ixO^S)/((rho2(ixO^S)  &
         +TwoMrQ(ixO^S))*Sqrt(-TwoMrQ(ixO^S)/(rho2(ixO^S)  &
         +TwoMrQ(ixO^S))+1/(1  &
         - (eqpar(a_)**2*sth(ixO^S)**2)/(eqpar(a_)**2  &
         + xKS(ixO^S,r_)**2))))

    ehatu(ixO^S,r_,r_) = 1/Sqrt(-TwoMrQ(ixO^S)/(rho2(ixO^S)  &
         + TwoMrQ(ixO^S))+1/(1  &
         - (1.0d0*eqpar(a_)**2*sth(ixO^S)**2)/(eqpar(a_)**2  &
         + xKS(ixO^S,r_)**2)))

    {^IFZ
    ehatu(ixO^S,z_,z_) = sqrt(rho2(ixO^S))
    }

    {^IFPHI
    ehatu(ixO^S,phi_,0) = (-TwoMrQ(ixO^S)*eqpar(a_)*sth(ixO^S)**2*(1  &
         + TwoMrQ(ixO^S)/rho2(ixO^S)))/((rho2(ixO^S)  &
         + TwoMrQ(ixO^S))*Sqrt(sth(ixO^S)**2*(eqpar(a_)**2  &
         + (TwoMrQ(ixO^S)*eqpar(a_)**2*sth(ixO^S)**2)/rho2(ixO^S)  &
         + xKS(ixO^S,r_)**2)))

    ehatu(ixO^S,phi_,r_) = (-1.0d0*eqpar(a_)*sth(ixO^S)**2*(1+  &
         TwoMrQ(ixO^S)/rho2(ixO^S)))/Sqrt(sth(ixO^S)**2*(eqpar(a_)**2  &
         + (TwoMrQ(ixO^S)*eqpar(a_)**2*sth(ixO^S)**2)/rho2(ixO^S)  &
         + xKS(ixO^S,r_)**2))

    ehatu(ixO^S,phi_,phi_) = sqrt(sth(ixO^S)**2*(eqpar(a_)**2  &
         + (TwoMrQ(ixO^S)*eqpar(a_)**2*sth(ixO^S)**2)/rho2(ixO^S) &
         + xKS(ixO^S,r_)**2))
    }
    
  end subroutine KNKSToLRNF
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
