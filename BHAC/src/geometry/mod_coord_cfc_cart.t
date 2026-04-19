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
! Routine to fill the components from a numerical metric, read from a file
! The numerical metric file's name must be 'numerical.met'
! Its header must contain:
! Metric name
! Metric variable names (metric and derivatives)
! Total number of points	Points in dim. 1	Points in dim. 2   ...
! Values of the metric variables
! The order in which the metric components and their
! derivatives should appear in the file is the following:
! Metric components:
! alpha beta1 beta2 beta3 g11 g12 g13 g22 g23 g33
! ...
!
! Hector Olivares
! 28.05.2017
!=============================================================================
module mod_coord_cfc_cart
use mod_metric_aux
use mod_cfc_parameters

character*20, parameter     :: coord="cfc_cart"

logical, save               :: init_from_g4=.false.
integer, parameter          :: ncoordpar=1
double precision, save      :: coordpar(ncoordpar)

procedure(sub_get_sqrtgamma_analytic), pointer        :: get_sqrtgamma_analytic       => null()
procedure(sub_get_alpha), pointer                     :: get_alpha                    => null()
procedure(sub_get_beta), pointer                      :: get_beta                     => null()
procedure(sub_get_g_component), pointer               :: get_g_component              => null()

  abstract interface

     subroutine sub_get_sqrtgamma_analytic(x^D,sqrtgamma,is_analytic,w_pt)
       include 'amrvacdef.f'
       double precision, intent(in)                     :: x^D
       double precision, intent(in), optional           :: w_pt(1:nw)
       double precision, intent(out)                    :: sqrtgamma
       logical, optional                                :: is_analytic
     end subroutine sub_get_sqrtgamma_analytic

     subroutine sub_get_alpha(x^D,alpha,iszero,dalphadj_iszero,dalphadj,jdir,w_pt)
       include 'amrvacdef.f'
       double precision, intent(in)                     :: x^D
       double precision, intent(in), optional           :: w_pt(1:nw)
       double precision, intent(out)                    :: alpha
       logical, optional, intent(out)                   :: iszero, dalphadj_iszero
       integer, optional, intent(in)                    :: jdir
       double precision, optional, intent(out)          :: dalphadj
     end subroutine sub_get_alpha

     subroutine sub_get_beta(idir,x^D,beta,iszero,dbetaidj_iszero,dbetaidj,jdir,w_pt)
       include 'amrvacdef.f'
       integer, intent(in)                      :: idir
       double precision, intent(in), optional   :: w_pt(1:nw)
       double precision, intent(in)             :: x^D
       double precision, intent(out)            :: beta
       logical, optional, intent(out)           :: iszero, dbetaidj_iszero
       double precision, optional, intent(out)  :: dbetaidj
       integer, optional, intent(in)            :: jdir
     end subroutine sub_get_beta
     
     subroutine sub_get_g_component(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir,w_pt)
       include 'amrvacdef.f'
       integer, intent(in)                      :: iin,jin
       double precision, intent(in), optional   :: w_pt(1:nw)
       double precision, intent(in)             :: x^D
       double precision, intent(out)            :: g
       logical, optional, intent(out)           :: iszero, dgdk_iszero
       double precision, optional, intent(out)  :: dgdk
       integer, optional, intent(in)            :: kdir
     end subroutine sub_get_g_component

  end interface

contains
!=============================================================================
  subroutine init_coord

    coordinate = Cartesian
    !only need cfc_cart is okay, no need cart
    get_gammainv_component_analytic => get_gammainv_component_analytic_cfc_cart
    get_gammainv_component => get_gammainv_component_analytic
    call switch_coord_to_cfc_cart
  end subroutine init_coord

  subroutine switch_coord_to_cfc_cart
      get_sqrtgamma_analytic  => get_sqrtgamma_analytic_cfc_cart
      get_alpha               => get_alpha_cfc_cart
      get_beta                => get_beta_cfc_cart
      get_g_component         => get_g_component_cfc_cart
  end subroutine switch_coord_to_cfc_cart
  
  subroutine switch_coord_to_cart
      get_sqrtgamma_analytic  => get_sqrtgamma_analytic_cart
      get_alpha               => get_alpha_cart
      get_beta                => get_beta_cart
      get_g_component         => get_g_component_cart
  end subroutine switch_coord_to_cart

  subroutine get_sqrtgamma_analytic_cfc_cart(x^D,sqrtgamma,is_analytic,w_pt)

    include 'amrvacdef.f'

    ! Since the metric is numeric, sqrtgamma is not calculated analytically.

    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    double precision, intent(in), optional           :: w_pt(1:nw)
    double precision                                 :: w_pt_local(1:nw)

    if(present(is_analytic)) is_analytic = .true.
    if(present(w_pt)) then
       w_pt_local = w_pt
    else
       w_pt_local(psi_metric_) = 1.0d0
    endif

    if (w_pt_local(psi_metric_) .lt. 1.0d0 &
            .or. w_pt_local(psi_metric_) .ge. 2.0d0) then
       write(*,*) w_pt_local(psi_metric_), 'psi_metric not filled in sqrtgamma'
       w_pt_local(psi_metric_) = 1.0d0
    endif

    sqrtgamma = 1.0d0 * w_pt_local(psi_metric_) **6.0d0

  end subroutine get_sqrtgamma_analytic_cfc_cart


!=============================================================================
  subroutine get_alpha_cfc_cart(x^D,alpha,iszero,dalphadj_iszero,dalphadj,jdir,w_pt)

    include 'amrvacdef.f'

    ! get the lapse.  Optional parameter is true if lapse is
    ! identically zero (does not really make sense)
    ! Optional parameters jdir and dalphadj request derivatives
    ! \partial_j \alpha ; j=jdir
    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: alpha
    logical, optional, intent(out)                   :: iszero, dalphadj_iszero
    integer, optional, intent(in)                    :: jdir
    double precision, optional, intent(out)          :: dalphadj
    double precision, intent(in), optional           :: w_pt(1:nw)
    ! .. local ..
    double precision, dimension(^ND)                 :: x
    double precision                                 :: d^Calpha
    double precision                                 :: dx^D
    double precision                                 :: alpha_local
    double precision                                 :: alpha_lower, alpha_upper
    double precision                                 :: w_pt_local(1:nw)


    if(present(w_pt)) then
           w_pt_local = w_pt
    else
           w_pt_local(alp_metric_) = 1.0d0
    endif

    if (w_pt_local(alp_metric_) .gt. 1.0d0 &
            .or. w_pt_local(alp_metric_) .le. 0.4d0) then   ! < 0.4 ?
       write(*,*) w_pt_local(alp_metric_), 'alp not well filled in get alp'
       w_pt_local(alp_metric_) = 1.0d0
    endif

    if(present(dalphadj) .and. .not. present(jdir) .or. &
         present(dalphadj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_alpha: derivatives requested without direction or output-slot given.")

    if(present(iszero)) iszero = .false.
 
      alpha = w_pt_local(alp_metric_) 
      if (.not.present(jdir)) then
      else  ! present jdir
          if (present(dalphadj)) then
             call mpistop('cfc not calucating metric dervatives here! here is get_alpha_cfc')
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .true. !.false.
          end if
      end if 

    if(alpha .ne. alpha) then
           write(*,*) 'alpha, w_pt_local(alp_metric_)'
           write(*,*) alpha, w_pt_local(alp_metric_)
         stop  "and yes, it's also a NaN alpha sometimes"
    endif

  end subroutine get_alpha_cfc_cart
  !=============================================================================
  subroutine get_beta_cfc_cart(idir,x^D,beta,iszero,dbetaidj_iszero,dbetaidj,jdir,w_pt)

    include 'amrvacdef.f'

    ! get the (contravariant!!) shift vector.
    ! The optional argument iszero is true if shift-component is 
    ! identically zero.
    ! if requested, dbetaidj is the derivative of the contravariant shift.
    ! \partial_j \beta^i ; i=idir, j=jdir
    integer, intent(in)                      :: idir
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: beta

    logical, optional, intent(out)           :: iszero, dbetaidj_iszero
    double precision, optional, intent(out)  :: dbetaidj
    integer, optional, intent(in)            :: jdir
    double precision, intent(in), optional   :: w_pt(1:nw)
    ! .. local ..
    double precision, dimension(^ND)         :: x
    double precision                         :: d^Cbeta
    double precision                         :: w_pt_local(1:nw)

    if(present(w_pt)) w_pt_local = w_pt

    if(present(dbetaidj) .and. .not. present(jdir) .or. &
         present(dbetaidj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_beta: derivatives requested &
         &without direction or output-slot given.")


    select case(idir)
     case(1)
        if (present(iszero)) iszero = .false.
          beta = w_pt_local(beta_metric1_)
          
          if (present(jdir)) then
             if (present(dbetaidj)) then
                call mpistop('cfc not calucating metric dervatives here! here is get_beta_cfc')
             end if
          endif

          if (present(dbetaidj_iszero)) then
             dbetaidj_iszero = .true. !.false.
          end if

     case(2)
        if (present(iszero)) iszero = .false.
           beta = w_pt_local(beta_metric2_)

          if (present(jdir)) then
             if (present(dbetaidj)) then
                call mpistop('cfc not calucating metric dervatives here! here is get_beta_cfc')
             end if
          endif

          if (present(dbetaidj_iszero)) then
             dbetaidj_iszero = .true. !.false.
          end if
     case(3)
        if (present(iszero)) iszero = .false.
          beta = w_pt_local(beta_metric3_)

          if (present(jdir)) then
             if (present(dbetaidj)) then
                call mpistop('cfc not calucating metric dervatives here! here is get_beta_cfc')
             end if
          endif

          if (present(dbetaidj_iszero)) then
             dbetaidj_iszero = .true. !.false.
          end if
    end select

    if(beta .ne. beta) then
           write(*,*) 'beta, w_pt_local(beta_metric1_),&
                       w_pt_local(beta_metric2_), w_pt_local(beta_metric3_)'
           write(*,*) beta, w_pt_local(beta_metric1_),&
                       w_pt_local(beta_metric2_), w_pt_local(beta_metric3_)
         stop  "and yes, it's also a NaN beta sometimes"
    endif


  end subroutine get_beta_cfc_cart
  !=============================================================================
  subroutine get_g_component_cfc_cart(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir,w_pt)

    include 'amrvacdef.f'

    ! This is at the heart of the scheme: Set the (spatial) metric components here
    ! and only here...
    ! Indices of the metric are down (covariant) g_{ij\}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dgdk and kdir request derivatives of the metric
    ! \partial_k g_{ij\} ; i=iin, j=jin, k=kdir
    ! The optional argument dgdk_iszero
    integer, intent(in)                      :: iin,jin
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: g

    logical, optional, intent(out)           :: iszero, dgdk_iszero
    double precision, optional, intent(out)  :: dgdk
    integer, optional, intent(in)            :: kdir
    double precision, intent(in), optional   :: w_pt(1:nw)
    ! .. local ..
    integer                                  :: i,j, component_
    double precision, dimension(^ND)         :: x
    double precision                         :: d^Cg
    double precision                         :: g_local
    double precision                         :: w_pt_local(1:nw)


    if(present(w_pt)) then 
         w_pt_local = w_pt
    else
         w_pt_local(psi_metric_) = 1.0d0
    endif

    if (w_pt_local(psi_metric_) .lt. 1.0d0 &
            .or. w_pt_local(psi_metric_) .ge. 2.0d0) then
       write(*,*) w_pt_local(psi_metric_), 'psi_metric not well filled in g_com'
       w_pt_local(psi_metric_) = 1.0d0
    endif

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

    if (present(iszero)) iszero = .false.
    if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.


    !initialize g
    g = 0.0d0
    if (i==j) then
       if(present(iszero)) iszero = .false.
       g = one * w_pt_local(psi_metric_)**4.0d0
    else
       if(present(iszero)) iszero = .true.
       g = zero
    end if

    if(g .ne. g) then 
           write(*,*) 'g, w_pt_local(psi_metric_)'
           write(*,*) g, w_pt_local(psi_metric_)
         stop  "and yes, it's also a NaN g sometimes" 
    endif

  end subroutine get_g_component_cfc_cart

  subroutine get_gammainv_component_analytic_cfc_cart(iin,jin,x^D,ginv,iszero,dginvdk_iszero,dginvdk,kdir,w_pt)

    include 'amrvacdef.f'

    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: ginv
    logical, optional, intent(out)           :: iszero, dginvdk_iszero
    double precision, optional, intent(out)  :: dginvdk
    double precision, optional, intent(in)   :: w_pt(1:nw)
    ! .. local ..
    double precision                         :: myrks, mythetaks, rho2
    double precision                         :: d1rks, d2thks, d1rho2, d2rho2
    integer                                  :: i, j
    double precision                         :: w_pt_local(1:nw)

    if(present(w_pt)) then 
         w_pt_local = w_pt
    else
         w_pt_local(psi_metric_) = 1.0d0
    endif
    if (w_pt_local(psi_metric_) .lt. 1.0d0 &
            .or. w_pt_local(psi_metric_) .ge. 2.0d0) then
       write(*,*) w_pt_local(psi_metric_), 'psi_metric not well filled in get gammainv'
       w_pt_local(psi_metric_) = 1.0d0
    endif



    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    if (present(iszero)) iszero = .false.
    if (present(dginvdk_iszero)) dginvdk_iszero = .true. !.false.

    !initialize g
    ginv = 0.0d0
    if (i==j) then
       if(present(iszero)) iszero = .false.
       ginv = one / (w_pt_local(psi_metric_)**4.0d0)
    else
       if(present(iszero)) iszero = .true.
       ginv = zero
    end if


    if(ginv .ne. ginv) then 
           write(*,*) 'ginv, w_pt_local(psi_metric_)'
           write(*,*) ginv, w_pt_local(psi_metric_)
         stop  "and yes, it's also a NaN ginv sometimes" 
    endif


   call mpistop('passed gammainv_componet analytic')

  end subroutine get_gammainv_component_analytic_cfc_cart



!=============================================================================
!  Below are flat metric
!=============================================================================

  subroutine get_sqrtgamma_analytic_cart(x^D,sqrtgamma,is_analytic,w_pt)

    include 'amrvacdef.f'

    ! Since the metric is numeric, sqrtgamma is not calculated analytically.

    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    double precision, intent(in), optional           :: w_pt(1:nw)


    if(present(is_analytic)) is_analytic = .true.

    sqrtgamma = one

  end subroutine get_sqrtgamma_analytic_cart

!=============================================================================
  subroutine get_alpha_cart(x^D,alpha,iszero,dalphadj_iszero,dalphadj,jdir,w_pt)

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
    double precision, intent(in), optional           :: w_pt(1:nw)
    ! .. local ..
    double precision, dimension(^ND)                 :: x
    double precision                                 :: d^Calpha
    double precision                                 :: dx^D
    double precision                                 :: alpha_lower, alpha_upper
    if(present(dalphadj) .and. .not. present(jdir) .or. &
         present(dalphadj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_alpha: derivatives requested without direction or output-slot given.")

    if(present(iszero)) iszero = .false.
    if (present(dalphadj_iszero)) then
       dalphadj_iszero = .true. !.false.
    end if
    if (present(dalphadj)) then
         call mpistop('cfc not calucating metric dervatives here! here is get_alpha_cart')
    endif

    alpha = 1.0d0

  end subroutine get_alpha_cart
  !=============================================================================
  subroutine get_beta_cart(idir,x^D,beta,iszero,dbetaidj_iszero,dbetaidj,jdir,w_pt)

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
    double precision, intent(in), optional   :: w_pt(1:nw)
    ! .. local ..
    double precision, dimension(^ND)         :: x
    double precision                         :: d^Cbeta
    if(present(dbetaidj) .and. .not. present(jdir) .or. &
         present(dbetaidj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_beta: derivatives requested &
         &without direction or output-slot given.")

    beta = 0.0d0
    if(present(iszero)) &
         iszero = .false.

    ! \partial_j \beta^i
    if (present(dbetaidj)) then
         call mpistop('cfc not calucating metric dervatives here! here is get_beta_cart')
    endif
    if (present(dbetaidj_iszero)) &
         dbetaidj_iszero = .true. !.false.

  end subroutine get_beta_cart
  !=============================================================================
  subroutine get_g_component_cart(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir,w_pt)

    include 'amrvacdef.f'

    ! This is at the heart of the scheme: Set the (spatial) metric components here
    ! and only here...
    ! Indices of the metric are down (covariant) g_{ij\}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dgdk and kdir request derivatives of the metric
    ! \partial_k g_{ij\} ; i=iin, j=jin, k=kdir
    ! The optional argument dgdk_iszero
    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: g
    logical, optional, intent(out)           :: iszero, dgdk_iszero
    double precision, optional, intent(out)  :: dgdk
    double precision, intent(in), optional   :: w_pt(1:nw)
    ! .. local ..
    integer                                  :: i,j, component_
    double precision, dimension(^ND)         :: x
    double precision                         :: d^Cg

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

       g = 0.0d0
    ! Diagonal unity
    if (i==j) then
       if(present(iszero)) iszero = .false.

       g = one

    else
       if(present(iszero)) iszero = .true.

       g = zero

    end if

    if (present(dgdk)) dgdk = zero
    if(present(dgdk_iszero)) dgdk_iszero = .true.


  end subroutine get_g_component_cart

  subroutine BLToCoord(ixI^L,ixO^L,xBL,xCoord)

    ! identity
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    integer                                                :: ix^D
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

    ! Trivial case
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xBL
    ! .. local ..
    integer                                                :: ix^D

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

    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to the current (BL) coordinates u4Coord.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    !
    !
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..


    if (eqpar(m_) .gt. zero) &
         call mpistop("BLToCoord: Cannot transform to Cartesian for non-zero mass")

    call mpistop("u4BLtoCoord: Transformation not implemented")

    u4Coord = u4BL



  end subroutine u4BLtoCoord
  !=============================================================================
  subroutine CoordToCart(ixI^L,ixO^L,x,xCart)

    ! Transforms the coordinates to "Boyer-Lindquist Cartesian coordinates"
    ! x = ra Sin(thetaBL) Cos(phiBL)
    ! y = ra Sin(thetaBL) Sin(phiBL)
    ! z = rBL Cos(thetaBL)
    ! where ra=sqrt(rBL**2+a**2) ! WARNING: I decided to put ra=rBL for now...
    ! We are already in BL-cordinates, so we only have to do one transformation.

    include 'amrvacdef.f'
    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCart
    ! .. local ..
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, r
    double precision, dimension(ixI^S,1:^ND)               :: xBL
    integer                                                :: zcart_, ycart_

    xCart = x

  end subroutine CoordToCart
  !=============================================================================
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)


    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^NC), intent(in)   :: u3Coord
    double precision, dimension(ixI^S,1:^NC), intent(out)  :: u3Cart
    double precision, dimension(ixI^S,1:^NC,1:^NC), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xBL
    double precision, dimension(ixI^S,1:^NC,1:^NC)         :: Jac
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, ra
    integer                                                :: zcart_, ycart_, ix,jx

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
  subroutine get_dgammainvdk(x^D,dgammainvdk,k)
    ! Obtains the derivatives of the inverse _spatial_ metric
    ! \partial_k gamma^{ij\} for a given k. 
    ! This is not required when initialising gamma, lapse and shift
    ! e.g. when set init_from_g4 = .false.

    double precision, intent(in)           :: x^D
    double precision, intent(out)          :: dgammainvdk(1:^NC,1:^NC)
    integer, intent(in)                    :: k
    ! .. local ..

    call mpistop("get_dgammainvdk: Not required and not implemented.")

    dgammainvdk = 0.0d0
    
  end subroutine get_dgammainvdk
  !=============================================================================
  ! End of coordinate-specific definitions.
  !=============================================================================

end module mod_coord_cfc_cart
