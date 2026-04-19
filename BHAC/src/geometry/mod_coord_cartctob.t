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
  ! Metric components for Cartesian coordinates -coord=cartctob
  ! This file is included by the preprocessor in mod_metric.t
  ! It will be part of the module mod_metric.
  !  
  ! Oliver Porth
  ! 28.06.2016
  !=============================================================================
  
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------
  
  character*20, parameter              :: coord="cartctob"
  logical, save                        :: init_from_g4 = .false.  ! We don't provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: ncoordpar = 0
  double precision, save               :: coordpar(ncoordpar) 

contains
!=============================================================================
  subroutine init_coord

    SPToCoord=>BLToCoord
    
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


    sqrtgamma = one

    if(present(is_analytic)) is_analytic = .false.

  end subroutine get_sqrtgamma_analytic
  !=============================================================================
  subroutine get_alpha(x^D,alpha,iszero,dalphadj_iszero,dalphadj,jdir)

    
    use mod_carpetgrid
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
    double precision                                 :: x_interpol1, x_interpol2, x_interpol3
    !-----------------------------------------------------------------------------
    if(present(dalphadj) .and. .not. present(jdir) .or. &
         present(dalphadj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_alpha: derivatives requested without direction or output-slot given.")

    ! TODO: ifnot3D mpistop!

    !assumes slicing at 0 for lower dimensional problems. TODO: obtain slice coordinate from amrvac.par?
    x_interpol1 = 0.d0
    x_interpol2 = 0.d0
    x_interpol3 = 0.d0

    {x_interpol^D = x^D;}
    alpha = finterpolate("alp", x_interpol1, x_interpol2, x_interpol3)

!      if (mype==0 .and. alpha>1.000001) then
!         print*,'-------------------------------------------------------------------------------'
!         write(*,'(a,f17.3, f17.3, f17.3, a, f17.3)')' alpha at point ',x1, x2, x3, 'is', alpha
!         print*,'-------------------------------------------------------------------------------'
!      end if

    if(present(iszero)) iszero = .false.
    if (present(dalphadj_iszero)) dalphadj_iszero = .false.
    if (present(dalphadj)) dalphadj = fderive("alp", x_interpol1, x_interpol2, x_interpol3, jdir)

  end subroutine get_alpha
  !=============================================================================
  subroutine get_beta(idir,x^D,beta,iszero,dbetaidj_iszero,dbetaidj,jdir)

    use mod_carpetgrid
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
    double precision                         :: x_interpol1, x_interpol2, x_interpol3
    !-----------------------------------------------------------------------------
    if(present(dbetaidj) .and. .not. present(jdir) .or. &
         present(dbetaidj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_beta: derivatives requested &
         &without direction or output-slot given.")

    x_interpol1 = 0.d0
    x_interpol2 = 0.d0
    x_interpol3 = 0.d0

    {x_interpol^D = x^D;}
 
    select case(idir)
    case(1)
        beta = finterpolate("betax", x_interpol1, x_interpol2, x_interpol3)
        if (present(dbetaidj)) dbetaidj = fderive("betax", x_interpol1, x_interpol2, x_interpol3, jdir)
    case(2) 
        beta = finterpolate("betay", x_interpol1, x_interpol2, x_interpol3)
        if (present(dbetaidj)) dbetaidj = fderive("betay", x_interpol1, x_interpol2, x_interpol3, jdir)
    case(3) 
        beta = finterpolate("betaz", x_interpol1, x_interpol2, x_interpol3)
        if (present(dbetaidj)) dbetaidj = fderive("betaz", x_interpol1, x_interpol2, x_interpol3, jdir)
   end select 

    if(present(iszero)) iszero = .false.
    ! \partial_j \beta^i
    if (present(dbetaidj_iszero)) dbetaidj_iszero = .false.

  end subroutine get_beta
  !=============================================================================
  subroutine get_g_component(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir)

    use mod_carpetgrid
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
    double precision                         :: x_interpol1, x_interpol2, x_interpol3
    !-----------------------------------------------------------------------------
    if(present(dgdk) .and. .not. present(kdir) .or. &
         present(dgdk_iszero) .and. .not. present(kdir)) &
         call mpistop("get_g_component: derivatives requested without &
         &direction or output-slot given.")


    x_interpol1 = 0.d0
    x_interpol2 = 0.d0
    x_interpol3 = 0.d0

    {x_interpol^D = x^D;}
 
    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    if(present(iszero)) iszero = .false.
    if(present(dgdk_iszero)) dgdk_iszero = .false.

  select case(i)
      case(1)

        select case(j)
            case(1)
              g = finterpolate("gxx", x_interpol1, x_interpol2, x_interpol3)
              if (present(dgdk)) then
                dgdk = fderive("gxx", x_interpol1, x_interpol2, x_interpol3, kdir)
              end if
            case(2)
              g = finterpolate("gxy", x_interpol1, x_interpol2, x_interpol3)
              if (present(dgdk)) then
                dgdk = fderive("gxy", x_interpol1, x_interpol2, x_interpol3, kdir)
              end if
            case(3)
              g = finterpolate("gxz", x_interpol1, x_interpol2, x_interpol3)
              if (present(dgdk)) then
                dgdk = fderive("gxz", x_interpol1, x_interpol2, x_interpol3, kdir)
              end if
        end select

      case(2)
        select case(j)
            case(2)
              g = finterpolate("gyy", x_interpol1, x_interpol2, x_interpol3)
              if (present(dgdk)) then
                dgdk = fderive("gyy", x_interpol1, x_interpol2, x_interpol3, kdir)
              end if


            case(3)
              g = finterpolate("gyz", x_interpol1, x_interpol2, x_interpol3)
              if (present(dgdk)) then
                dgdk = fderive("gyz", x_interpol1, x_interpol2, x_interpol3, kdir)
              end if

       end select

      case(3)
              g = finterpolate("gzz", x_interpol1, x_interpol2, x_interpol3)
              if (present(dgdk)) then
                dgdk = fderive("gzz", x_interpol1, x_interpol2, x_interpol3, kdir)
              end if
    end select


  end subroutine get_g_component
  !=============================================================================
    double precision function outerhorizon()
    
    !-----------------------------------------------------------------------------

    outerhorizon = 0.0d0
    
  end function outerhorizon
  !=============================================================================
  subroutine BLToCoord(ixI^L,ixO^L,xBL,xCoord)

    ! Assuming zero mass and do the transformation for spherical coordinates
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    !-----------------------------------------------------------------------------

    xCoord(ixO^S,1) = xBL(ixO^S,1) {^IFPHIIN * cos(xBL(ixO^S,^PHI))} {^IFZIN * sin(xBL(ixO^S,^Z))}

    if (^ND .lt. 3) then
       ! Simulation in a plane
       if (^PHI .le. ^ND .and. ^PHI .gt. 0) then
          {^NOONED
          !Simulation in the r-phi plane:
          xCoord(ixO^S,2) = xBL(ixO^S,1) {^IFPHIIN * sin(xBL(ixO^S,^PHI))} {^IFZIN * sin(xBL(ixO^S,^Z))}
          }
       else if (^Z .le. ^ND .and. ^Z .gt. 0) then
          !Simulation in the r-z plane:
          {^NOONED
          xCoord(ixO^S,2) = xBL(ixO^S,1) {^IFZIN * cos(xBL(ixO^S,^Z))}
          }
       end if
    else
       ! Simulation in 3D
       {#IFDEF D3
       xCoord(ixO^S,2) = xBL(ixO^S,1) {^IFPHIIN * sin(xBL(ixO^S,^PHI))} {^IFZIN * sin(xBL(ixO^S,^Z))}
       xCoord(ixO^S,3) = xBL(ixO^S,1) {^IFZIN * cos(xBL(ixO^S,^Z))}
       }
    end if
    
  end subroutine BLToCoord
  !=============================================================================
  subroutine CoordToBL(ixI^L,ixO^L,xCoord,xBL)

    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xBL
    ! .. local ..
    !-----------------------------------------------------------------------------

    xBL(ixO^S,1) = sqrt({^D& xCoord(ixO^S,^D)**2 |+})

    if (^ND .lt. 3) then
       ! Simulation in a plane
       if (^PHI .le. ^ND .and. ^PHI .gt. 0) then 
          !Simulation in the r-phi plane:
          {^IFPHIIN xBL(ixO^S,^PPHI) = atan2(xCoord(ixO^S,2),xCoord(ixO^S,1)) }
       else if (^Z .le. ^ND .and. ^Z .gt. 0) then
          !Simulation in the r-z plane:
          {^IFZIN xBL(ixO^S,^ZZ) = acos(xCoord(ixO^S,^ZZ)) / xBL(ixO^S,1) }
       end if
    else
       ! Simulation in 3D
       {^IFPHIIN xBL(ixO^S,^PPHI) = atan2(xCoord(ixO^S,2),xCoord(ixO^S,1)) }
       {^IFZIN xBL(ixO^S,^ZZ)   = acos(xCoord(ixO^S,^ZZ)) / xBL(ixO^S,1) }
    end if
    
  end subroutine CoordToBL
  !=============================================================================
  subroutine u4BLtoCoord(ixI^L,ixO^L,x,u4BL,u4Coord)

    ! Assuming zero mass and do the transformation for spherical coordinates
    ! Most likely we will not use this ever.  
    !
    include 'amrvacdef.f'

    integer                                                :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    !-----------------------------------------------------------------------------

    if (eqpar(m_) .gt. zero) &
         call mpistop("BLToCoord: Cannot transform to Cartesian for non-zero mass")

    call mpistop("u4BLtoCoord: Transformation not implemented")

    u4Coord = u4BL
    
  end subroutine u4BLtoCoord
  !============================================================================= 
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)

    ! u3Cart = u3Coord
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^NC), intent(in)   :: u3Coord
    double precision, dimension(ixI^S,1:^NC), intent(out)  :: u3Cart
    double precision, dimension(ixI^S,1:^NC,1:^NC), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixI^S,1:^NC,1:^NC)         :: Jac
    integer                                                :: ix, jx
    !-----------------------------------------------------------------------------

    u3Cart = u3Coord
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Assemble (diagonal) Jacobian:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (present(J)) then
       Jac(ixO^S,:,:) = zero
       do ix=1,ndir
          do jx=1,ndir
             Jac(ixO^S,ix,jx) = one
          end do
       end do
       J = Jac
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end subroutine u3CoordToCart
  !=============================================================================
  subroutine CoordToCart(ixI^L,ixO^L,x,xCart)

    ! identity
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCart
    ! .. local ..
    !-----------------------------------------------------------------------------

    xCart = x

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
