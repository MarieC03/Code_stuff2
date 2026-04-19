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
  ! Metric components for Cartesian Kerr-Schild coordinates -coord=cylks
  !
  ! Cylindrical Kerr-Schild coordinates, related to Cartesian Kerr-Schild
  ! coordinates by R = x^2 + y^2 , phi = acrtan(y/x).
  !
  ! Hector Olivares
  ! 26.05.2019
  !=============================================================================
  
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------
  
  character*20, parameter              :: coord="cylks"
  logical, save                        :: init_from_g4 = .true.  ! We provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put e.g. the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: rexc_= 1, rcut_=2             ! Excision radius
  integer, parameter                   :: ncoordpar = 2
  double precision, save               :: coordpar(ncoordpar) 


contains
!=============================================================================
  subroutine init_coord

    SPToCoord => KSToCoord
    u4CoordToKS => CYLKS_u4CoordToKS
    u4CoordToBL => CYLKS_u4CoordToBL
    d4CoordToKS => CYLKS_d4CoordToKS
    u4CoordToBL => CYLKS_u4CoordToBL
   
  end subroutine init_coord
  !=============================================================================
  subroutine get_helpers(rcoord,z,r,r2,tmr,idel)
    include 'amrvacdef.f'

    double complex, intent(in)          :: rcoord, z
    double complex, intent(out)         :: r, r2, tmr, idel
    ! .. local ..
    double complex                      :: a2, rho
    double complex                      :: sig, del, isig
    !-----------------------------------------------------------------------------

    a2  = eqpar(a_)*eqpar(a_)
    rho = rcoord*rcoord + z*z - a2
    r   = SQRT(0.5D0*(rho + SQRT(rho*rho + 4.D0*a2*z*z)))

    r2  = r*r
    tmr = 2.D0*eqpar(m_)*r

    sig  = r2*r2 + a2*z*z
    del  = r2 + a2

    ! A small number is added in order to remove the singularity,
    ! but only well inside the horizon (same as for Cartesian Kerr-Schild)
    if (abs(r).lt.1.1d0*eqpar(m_)) then
        isig = 1.D0/(sig + 1.0d-3)
    else
        isig = 1.D0/sig
    end if

    tmr  = tmr*isig

    idel = 1.D0/(del + 1.0d-7)

  end subroutine get_helpers
  !=============================================================================
  double complex function g00(x^D)
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g00 = -(1.D0 - tmr*r2)

  end function g00
  !=============================================================================
  double complex function g01(x^D)
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g01 = tmr*rcoord*r2*r*idel

  end function g01
  !=============================================================================
  double complex function g0^PHI(x^D)
    ! 0-phi component
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g0^PHI = -tmr*r2*r*rcoord*rcoord*eqpar(a_)*idel

  end function g0^PHI
  !=============================================================================
  double complex function g0^Z(x^D)
    ! 0-z component
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g0^Z = tmr*r2*z

  end function g0^Z
  !=============================================================================
  double complex function g11(x^D)
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g11 = 1.D0 + tmr*r2*r2*rcoord*rcoord*idel*idel

  end function g11
  !=============================================================================
  double complex function g1^PHI(x^D)
    ! r-phi component
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g1^PHI =-tmr*r2*rcoord*rcoord*rcoord*r*eqpar(a_)*idel*idel

  end function g1^PHI
  !=============================================================================
  double complex function g1^Z(x^D)
    ! r-z component
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g1^Z = tmr*r2*rcoord*z*idel
    
  end function g1^Z
  !=============================================================================
  double complex function g^PHI^PHI(x^D)
    ! phi-phi component
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED0
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g^PHI^PHI = 1.D0 + tmr*r2*rcoord*rcoord*r2*eqpar(a_)*eqpar(a_)*idel*idel

  end function g^PHI^PHI
  !=============================================================================
  double complex function g23(x^D)
    ! phi-z component
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g23 =-tmr*r*rcoord*rcoord*eqpar(a_)*z*idel

  end function g23
  !=============================================================================
  double complex function g^Z^Z(x^D)
    ! phi-phi component
    include 'amrvacdef.f'
    double complex, intent(in)          :: x^D
    ! .. local ..
    double complex :: rcoord, z
    double complex :: r, r2, tmr, idel
    !-----------------------------------------------------------------------------

    rcoord = x1
    {^NOONED
    z = x^Z
    }{^IFONED
    z = 0.0d0
    }

    call get_helpers(rcoord,z,r,r2,tmr,idel)

    g^Z^Z = 1.D0 + tmr*z*z

  end function g^Z^Z
  !=============================================================================
  double complex function g(iin,jin,x^D)
    include 'amrvacdef.f'
    integer, intent(in)                 :: iin, jin
    double complex, intent(in)          :: x^D
    ! .. local ..
    integer                             :: i, j
    !-----------------------------------------------------------------------------

    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    select case(i)
    case(0)
       select case(j)
       case(0)
          g = g00(x^D)
       case(1)
          g = g01(x^D)
       case(2)
          g = g02(x^D)
       case(3)
          g = g03(x^D)
       end select
    case(1)
       select case(j)
       case(1)
          g = g11(x^D)
       case(2)
          g = g12(x^D)
       case(3)
          g = g13(x^D)
       end select
    case(2)
       select case(j)
       case(2)
          g = g22(x^D)
       case(3)
          g = g23(x^D)
       end select
    case(3)
       select case(j)
       case(3)
          g = g33(x^D)
       end select
    end select


  end function g
  !=============================================================================
  subroutine get_dgammainvdk(x^D,dgammainvdk,k)
    ! Obtains the derivatives of the inverse _spatial_ metric
    ! \partial_k gamma^{ij} for a given k. 
    ! This is required to obtain derivatives e.g. of the contravariant shift
    ! and the lapse.  

    use mod_lu
    include 'amrvacdef.f'

    double precision, intent(in)           :: x^D
    double precision, intent(out)          :: dgammainvdk(1:^NC,1:^NC)
    integer, intent(in)                    :: k
    ! .. local ..
    integer                                :: i, j
    integer, dimension(^NC)                :: indx
    integer                                :: code, d
    double complex                         :: xloc^D
    double complex, dimension(1:^NC,1:^NC) :: a
    double complex, dimension(1:^NC)       :: b
    !-----------------------------------------------------------------------------

    xloc^D=cmplx(x^D , kr(k,^D)*smalldouble , kind(1.d0));

    do i=1,^NC
       do j=i,^NC
          a(i,j) = g(i,j,xloc^D)
       end do
    end do
    do i=1,^NC
       do j=i+1,^NC
          a(j,i) = a(i,j)
       end do
    end do

    call LUDCMP(a,^NC,indx,d,code)

    do j=1,^NC
       ! Directly calculate inverse metric:
       b    = cmplx(0.0d0,0.0d0,kind(1.0d0))
       b(j) = cmplx(1.0d0,0.0d0,kind(1.0d0))
       call LUBKSB(a,^NC,indx,b)
       do i=1,j
          dgammainvdk(i,j) = aimag(b(i))/smalldouble
       end do
    end do
    do j=1,^NC
       do i=j+1,^NC
          dgammainvdk(i,j) = dgammainvdk(j,i)
       end do
    end do

  end subroutine get_dgammainvdk
  !=============================================================================
  subroutine get_sqrtgamma_analytic(x^D,sqrtgamma,is_analytic)

    include 'amrvacdef.f'

    ! At the moment sqrtgamma is calculated numerically. Perhaps later
    ! I will obtain the analytic expression to speed-up the initialization. 

    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    !-----------------------------------------------------------------------------

    if(present(is_analytic)) is_analytic = .false.

    sqrtgamma = 0.0

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
    double precision                                 :: myalpha2  
    !-----------------------------------------------------------------------------
    if(present(dalphadj) .and. .not. present(jdir) .or. &
         present(dalphadj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_alpha: derivatives requested without direction or output-slot given.")

    if(present(iszero)) iszero = .false.

    if (present(jdir)) then
       if (present(dalphadj)) then
          call mpistop("get_alpha: derivatives not implemented here.")
       end if
       if (present(dalphadj_iszero)) then
          dalphadj_iszero = .false.
       end if
    end if

    if(.not.present(iszero) .and. .not. present(dalphadj_iszero)) &
         call mpistop("get_alpha: should only be used to check whether coefficient is zero.")

    ! return dummy:
    myalpha = 0.0d0

  end subroutine get_alpha
  !=============================================================================
  subroutine get_beta(idir,x^D,mybeta,iszero,dbetaidj_iszero,dbetaidj,jdir)

    include 'amrvacdef.f'

    ! get the (contravariant!!) shift vector.
    ! The optional argument iszero is true if shift-component is 
    ! identically zero.
    ! if requested, dbetaidj is the derivative of the contravariant shift.
    ! \partial_j \beta^i ; i=idir, j=jdir
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

    if(present(iszero)) iszero = .false.

    if (present(jdir)) then
       if (present(dbetaidj)) then
          call mpistop("get_beta: derivatives not implemented here.")
       end if
       if (present(dbetaidj_iszero)) then
          dbetaidj_iszero = .false.
       end if
    end if


    {^IFZOUT
    ! no z-shift or derivatives in the xy-plane:
    if (idir .eq. ^Z) then 
       if(present(iszero)) iszero = .true.
       mybeta = 0.0d0
       if (present(dbetaidj)) dbetaidj = 0.0d0
       if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.
    end if
    }

    if(.not.present(iszero) .and. .not. present(dbetaidj_iszero)) &
         call mpistop("get_beta: should only be used to check whether coefficient is zero.")

    ! return dummy:
    mybeta = 0.0d0

  end subroutine get_beta
  !=============================================================================
  subroutine get_g_component(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir)

    include 'amrvacdef.f'

    ! This is at the heart of the scheme: Since init_from_g4 = .true.,
    ! we set the full four-metric components here and only here...
    ! Indices of the metric are down (covariant) g_{\mu \nu}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dgdk and kdir request derivatives of the metric
    ! \partial_k g_{ij} ; i=iin, j=jin, k=kdir
    ! The optional argument dgdk_iszero is true if the element is identically zero.
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
    case(0) ! t-something coordinate
       select case(j)
       case(0) ! tt component
          g = g00({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g00,kdir,dgdk)          
       case(1) ! tx component
          g = g01({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g01,kdir,dgdk)
       case(2) ! ty component
          g = g02({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g02,kdir,dgdk)
       case(3) ! tz component
          g = g03({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g03,kdir,dgdk)
       end select
    case(1) ! x-something coordinate
       select case(j)
       case(1) ! xx component
          g = g11({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g11,kdir,dgdk)
       case(2) ! xy component
          g = g12({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g12,kdir,dgdk)
       case(3) ! xz component
          g = g13({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g13,kdir,dgdk)
       end select
    case(2) ! y-something coordinate
       select case(j)
       case(2) ! yy coordinate
          g = g22({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g22,kdir,dgdk)
       case(3) ! yz coordinate
          g = g23({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g23,kdir,dgdk)          
       end select
    case(3) ! z-something coordinate
       select case(j)
       case(3) ! zz coordinate
          g = g33({^D& cmplx(x^D, 0.0d0, kind(1.d0))|,})
          if (present(dgdk)) call complex_derivative(x^D,g33,kdir,dgdk)
       end select
    end select

    if(present(iszero))       iszero      = .false.
    if(present(dgdk_iszero))  dgdk_iszero = .false.
    {^IFZOUT ! We are actually in the xy plane and dz is zero:
    if (present(kdir)) then
       if (kdir .eq. 3) then
          if(present(dgdk_iszero))  dgdk_iszero = .true.
       end if
    end if
    }

  end subroutine get_g_component
  !=============================================================================
  double precision function outerhorizon()

    include 'amrvacdef.f'
    !-----------------------------------------------------------------------------

    outerhorizon = eqpar(m_) + sqrt(eqpar(m_)**2 - eqpar(a_)**2)

  end function outerhorizon
  !=============================================================================
  double precision function innerhorizon()

    include 'amrvacdef.f'
    !-----------------------------------------------------------------------------

    innerhorizon = eqpar(m_) - sqrt(eqpar(m_)**2 - eqpar(a_)**2)

  end function innerhorizon

  !=============================================================================
  subroutine CoordToBL(ixI^L,ixO^L,xCoord,xBL)

    use mod_transform, only: CylKSToKS, KSToBL
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xBL
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    !-----------------------------------------------------------------------------

    call CylKSToKS(ixI^L,ixO^L,xCoord,xKS)
    call KSToBL(ixI^L,ixO^L,xKS,xBL)

  end subroutine CoordToBL
  !=============================================================================
  subroutine CoordToKS(ixI^L,ixO^L,xCoord,xKS)

    use mod_transform, only: CylKSToKS
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xKS
    ! .. local ..
    !-----------------------------------------------------------------------------

    call CylKSToKS(ixI^L,ixO^L,xCoord,xKS)

  end subroutine CoordToKS
  !=============================================================================
  subroutine u4KStoCoord(ixI^L,ixO^L,xKS,u4KS,u4CylKS,J)

    use mod_transform, only: u4KStoCylKS
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4KS
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4CylKS
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    !-----------------------------------------------------------------------------

    if (present(J)) then 
       call u4KSToCylKS(ixI^L,ixO^L,xKS,u4KS,u4CylKS,J=J)
    else
       call u4KSToCylKS(ixI^L,ixO^L,xKS,u4KS,u4CylKS)
    end if
    
  end subroutine u4KStoCoord
  !=============================================================================  
  subroutine CYLKS_u4CoordToKS(ixI^L,ixO^L,xCylKS,u4CylKS,u4KS,J)

    use mod_transform, only: u4CylKSToKS
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCylKS
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4CylKS
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4KS
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    !-----------------------------------------------------------------------------
    
    if (present(J)) then 
       call u4CylKSToKS(ixI^L,ixO^L,xCylKS,u4CylKS,u4KS,J=J)
    else
       call u4CylKSToKS(ixI^L,ixO^L,xCylKS,u4CylKS,u4KS)
    end if
    
  end subroutine CylKS_u4CoordToKS
  !=============================================================================  
  subroutine CYLKS_d4CoordToKS(ixI^L,ixO^L,xCylKS,d4CylKS,d4KS)

    use mod_transform, only: d4CylKSToKS
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCylKS
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: d4CylKS
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: d4KS
    ! .. local ..
    !-----------------------------------------------------------------------------
    
    call d4CylKSToKS(ixI^L,ixO^L,xCylKS,d4CylKS,d4KS)
    
  end subroutine CylKS_d4CoordToKS
  !=============================================================================
  subroutine CYLKS_u4CoordToBL(ixI^L,ixO^L,xCylKS,u4CylKS,u4BL,J)

    use mod_transform, only: u4CylKSToBL
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCylKS
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4CylKS
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4BL
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    !-----------------------------------------------------------------------------
    
    if (present(J)) then 
       call u4CylKSToBL(ixI^L,ixO^L,xCylKS,u4CylKS,u4BL,J=J)
    else
       call u4CylKSToBL(ixI^L,ixO^L,xCylKS,u4CylKS,u4BL)
    end if
    
  end subroutine CYLKS_u4CoordToBL
  !=============================================================================
  subroutine u4BLtoCoord(ixI^L,ixO^L,xBL,u4BL,u4Coord)

    use mod_transform, only: BLToKS, u4BLToKS, u4KStoCylKS
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    double precision, dimension(ixI^S,0:^NC)               :: u4KS
    !-----------------------------------------------------------------------------

    ! BL to KS in spherical coordinates
    call u4BLtoKS(ixI^L,ixO^L,xBL,u4BL,u4KS)

    ! Coordinate transformation from spherical BL to spherical KS
    call BLToKS(ixI^L,ixO^L,xBL,xKS)

    ! tranformation of u4 from spherical to cartesian (KS)
    call u4KStoCylKS(ixI^L,ixO^L,xKS,u4KS,u4Coord)

  end subroutine u4BLToCoord
  !============================================================================= 
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)

    ! Transforms any contravariant three vector u3Coord in the current coordinates
    ! to a vector in the Cartesian basis.
    ! We just return the same Cartesian KS coordinates.
    ! Better would be to go the route CylKS->BL->Cartesian. 
    !
    !     | 1   0   0 |
    ! J = | 0   1   0 |
    !     | 0   0   1 |
    !
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
    do ix=1,^NC
       Jac(ixO^S,ix,ix) = one
    end do

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

    ! We just return the same Cartesian KS coordinates.
    ! Better would be to go the route CylKS->BL->Cartesian
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCart
    ! .. local ..
    !-----------------------------------------------------------------------------

    {^IFPHIIN
    xCart(ixO^S,1) = x(ixO^S,1)*cos(x(ixO^S,^PHI))
    xCart(ixO^S,2) = x(ixO^S,1)*sin(x(ixO^S,^PHI))
    }{^IFPHIOUT
    xCart(ixO^S,1) = x(ixO^S,1)
    }

    {^IFZIN
    xCart(ixO^S,^Z) = x(ixO^S,^Z)
    }

  end subroutine CoordToCart
  !=============================================================================
  subroutine KSToCoord(ixI^L,ixO^L,xKS,xCylKS)

    ! Alejandro Cruz 12.07.2017
    use mod_transform, only: KSToCylKS
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCylKS
    !-----------------------------------------------------------------------------
    
    call KSToCylKS(ixI^L,ixO^L,xKS,xCylKS)

  end subroutine KSToCoord
  !=============================================================================
  subroutine BLToCoord(ixI^L,ixO^L,xBL,xCoord)

    ! Alejandro Cruz 12.07.2017
    use mod_transform, only: BLToKS

    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord

    double precision, dimension(ixI^S,1:^ND)              :: xKS

    call BLToKS(ixI^L,ixO^L,xBL,xKS)
    call KSToCoord(ixI^L,ixO^L,xKS, xCoord)

  end subroutine BLToCoord
  !=============================================================================
  ! End of coordinate-specific definitions.
  !=============================================================================
