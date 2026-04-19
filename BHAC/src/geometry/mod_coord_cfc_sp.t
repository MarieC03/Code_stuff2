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
module mod_coord_cfc_sp
use mod_metric_aux
use mod_cfc_parameters
character*20, parameter     :: coord="cfc_sp"

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

    coordinate = spherical
    !only need cfc_sp is okay, no need sp
    get_gammainv_component_analytic => get_gammainv_component_analytic_cfc_sp
    get_gammainv_component => get_gammainv_component_analytic

    call switch_coord_to_cfc_sp
  end subroutine init_coord

  subroutine switch_coord_to_cfc_sp
      get_sqrtgamma_analytic  => get_sqrtgamma_analytic_cfc_sp
      get_alpha               => get_alpha_cfc_sp
      get_beta                => get_beta_cfc_sp
      get_g_component         => get_g_component_cfc_sp
  end subroutine switch_coord_to_cfc_sp
  
  subroutine switch_coord_to_sp
      get_sqrtgamma_analytic  => get_sqrtgamma_analytic_sp
      get_alpha               => get_alpha_sp
      get_beta                => get_beta_sp
      get_g_component         => get_g_component_sp
  end subroutine switch_coord_to_sp

  subroutine get_sqrtgamma_analytic_cfc_sp(x^D,sqrtgamma,is_analytic,w_pt)

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


    !code test, includes a psi**6
    {^IFZIN
    sqrtgamma = x1**2 * abs(sin(x^Z)) * w_pt_local(psi_metric_) **6.0d0
    }{^IFZOUT
    sqrtgamma = x1**2 * w_pt_local(psi_metric_) **6.0d0
    }

  end subroutine get_sqrtgamma_analytic_cfc_sp


!=============================================================================
  subroutine get_alpha_cfc_sp(x^D,alpha,iszero,dalphadj_iszero,dalphadj,jdir,w_pt)

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

  end subroutine get_alpha_cfc_sp
  !=============================================================================
  subroutine get_beta_cfc_sp(idir,x^D,beta,iszero,dbetaidj_iszero,dbetaidj,jdir,w_pt)

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


  end subroutine get_beta_cfc_sp


  subroutine get_g_component_cfc_sp(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir,w_pt)

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

    ! Diagonal for now
    if (i==j) then
       if(present(iszero)) &
            iszero = .false.

       if (i .eq. 1) then ! r direction
          g = 1.0d0 * w_pt_local(psi_metric_) ** 4.0d0
          if (present(kdir)) then
                if (present(dgdk_iszero)) dgdk_iszero = .true.
          end if

       else if (i .eq. ^PHI) then ! phi direction
          {^IFZIN
          g = (x1*sin(x^Z))**2 * w_pt_local(psi_metric_) ** 4.0d0
          }{^IFZOUT
          ! equatorial region
          g = x1**2 * w_pt_local(psi_metric_) ** 4.0d0
          }

          if (present(kdir)) then
                if (present(dgdk_iszero)) dgdk_iszero = .true.
          end if

       else if (i .eq. ^Z) then ! theta direction
          g = x1**2 * w_pt_local(psi_metric_) ** 4.0d0

          if (present(kdir)) then
                if (present(dgdk_iszero)) dgdk_iszero = .true.
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

    if(g .ne. g) then 
           write(*,*) 'g, w_pt_local(psi_metric_)'
           write(*,*) g, w_pt_local(psi_metric_)
         stop  "and yes, it's also a NaN g sometimes" 
    endif

  end subroutine get_g_component_cfc_sp


!  !=============================================================================
!  subroutine get_g_component_cfc_sp(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir,w_pt)
!
!    include 'amrvacdef.f'
!
!    ! This is at the heart of the scheme: Set the (spatial) metric components here
!    ! and only here...
!    ! Indices of the metric are down (covariant) g_{ij\}
!    ! The optional argument iszero is true if the element is identically zero
!    ! The optional arguments dgdk and kdir request derivatives of the metric
!    ! \partial_k g_{ij\} ; i=iin, j=jin, k=kdir
!    ! The optional argument dgdk_iszero
!    integer, intent(in)                      :: iin,jin
!    double precision, intent(in)             :: x^D
!    double precision, intent(out)            :: g
!
!    logical, optional, intent(out)           :: iszero, dgdk_iszero
!    double precision, optional, intent(out)  :: dgdk
!    integer, optional, intent(in)            :: kdir
!    double precision, intent(in), optional   :: w_pt(1:nw)
!    ! .. local ..
!    integer                                  :: i,j, component_
!    double precision, dimension(^ND)         :: x
!    double precision                         :: d^Cg
!    double precision                         :: g_local
!    double precision                         :: w_pt_local(1:nw)
!
!
!    if(present(w_pt)) then 
!         w_pt_local = w_pt
!    else
!         w_pt_local(psi_metric_) = 1.0d0
!    endif
!
!    if (w_pt_local(psi_metric_) .lt. 1.0d0 &
!            .or. w_pt_local(psi_metric_) .ge. 2.0d0) then
!       write(*,*) w_pt_local(psi_metric_), 'psi_metric not well filled in g_com'
!       w_pt_local(psi_metric_) = 1.0d0
!    endif
!
!    if(present(dgdk) .and. .not. present(kdir) .or. &
!         present(dgdk_iszero) .and. .not. present(kdir)) &
!         call mpistop("get_g_component: derivatives requested without &
!         &direction or output-slot given.")
!
!    ! metric is symmetric: swap indices if needed:
!    ! User needs only to provide values for i<=j (upper triangle).  
!    if (iin>jin) then
!       i=jin; j=iin
!    else
!       i=iin; j=jin
!    end if
!
!    if (present(iszero)) iszero = .false.
!    if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
!
!
!    !initialize g
!    g = 0.0d0
!
!    ! set the component index 
!    select case(i)
!      case(1)
!        select case(j)
!        case(1)
!         g = 1.0d0 * (w_pt_local(psi_metric_) ** 4.0d0)
!         if (present(kdir)) then
!             if (present(dgdk)) then
!                call mpistop('cfc not calucating metric dervatives here! here is get_g_com_cfc')
!             end if
!             if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
!         endif
!        case(2)
!          g = 0.0d0
!          if(present(iszero))       iszero      = .true.
!          if (present(dgdk))        dgdk        = 0.0d0
!          if(present(dgdk_iszero))  dgdk_iszero = .true.
!        case(3) 
!          g = 0.0d0
!          if(present(iszero))       iszero      = .true.
!          if (present(dgdk))        dgdk        = 0.0d0
!          if(present(dgdk_iszero))  dgdk_iszero = .true.
!        end select 
!
!      case(2)
!        select case(j)
!        case(2)
!          g = x1**2 * (w_pt_local(psi_metric_)** 4.0d0)
!         if (present(kdir)) then
!             if (present(dgdk)) then
!                call mpistop('cfc not calucating metric dervatives here! here is get_g_com_cfc')
!             end if
!             if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
!         endif
!        case(3)
!          g = 0.0d0
!          if(present(iszero))       iszero      = .true.
!          if (present(dgdk))        dgdk        = 0.0d0
!          if(present(dgdk_iszero))  dgdk_iszero = .true.
!        end select
!
!      case(3)
!       {^IFZIN
!         g = (x1*dsin(x2))**2 * (w_pt_local(psi_metric_) ** 4.0d0)  ! phi-direc: r**2 * theta**2 * psi**4
!       }{^IFZOUT
!       ! equatorial region
!         g = x1**2 * (w_pt_local(psi_metric_) ** 4.0d0)
!       }
!
!         if (present(kdir)) then
!             if (present(dgdk)) then
!                call mpistop('cfc not calucating metric dervatives here! here is get_g_com_cfc')
!             end if
!             if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
!         endif
!    end select
!
!    if(g .ne. g) then 
!           write(*,*) 'g, w_pt_local(psi_metric_)'
!           write(*,*) g, w_pt_local(psi_metric_)
!         stop  "and yes, it's also a NaN g sometimes" 
!    endif
!
!  end subroutine get_g_component_cfc_sp

  subroutine get_gammainv_component_analytic_cfc_sp(iin,jin,x^D,ginv,iszero,dginvdk_iszero,dginvdk,kdir,w_pt)

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

    ! set the component index 
    select case(i)
      case(1)
        select case(j)
        case(1)
         ginv = 1.0d0 / (1.0d0 * (w_pt_local(psi_metric_) ** 4.0d0))
         if (present(kdir)) then
             if (present(dginvdk)) then
                call mpistop('cfc not calucating metric dervatives here! here is get_ginv_cfc')
             end if
             if (present(dginvdk_iszero)) dginvdk_iszero = .true. !.false.
         endif
        case(2)
          ginv = 0.0d0
          if(present(iszero))       iszero            = .true.
          if (present(dginvdk))     dginvdk           = 0.0d0
          if(present(dginvdk_iszero))  dginvdk_iszero = .true.
        case(3) 
          ginv = 0.0d0  
          if(present(iszero))       iszero            = .true.
          if (present(dginvdk))     dginvdk           = 0.0d0
          if(present(dginvdk_iszero))  dginvdk_iszero = .true.
        end select 

      case(2)
        select case(j)
        case(2)
          ginv = 1.0d0 / (x1**2 * (w_pt_local(psi_metric_)** 4.0d0))
         if (present(kdir)) then
             if (present(dginvdk)) then
                call mpistop('cfc not calucating metric dervatives here! here is get_ginv_cfc')
             end if
             if (present(dginvdk_iszero)) dginvdk_iszero = .true. !.false.
         endif
        case(3)
          ginv = 0.0d0
          if(present(iszero))          iszero         = .true.
          if (present(dginvdk))        dginvdk        = 0.0d0
          if(present(dginvdk_iszero))  dginvdk_iszero = .true.
        end select

      case(3)
       {^IFZIN
         ginv = 1.0d0 / ((x1*dsin(x2))**2 * (w_pt_local(psi_metric_) ** 4.0d0))  ! phi-direc: r**2 * theta**2 * psi**4
       }{^IFZOUT
       ! equatorial region
         ginv = 1.0d0 / (x1**2 * (w_pt_local(psi_metric_) ** 4.0d0))
       }

         if (present(kdir)) then
             if (present(dginvdk)) then
                call mpistop('cfc not calucating metric dervatives here! here is get_ginv_cfc')
             end if
             if (present(dginvdk_iszero)) dginvdk_iszero = .true. !.false.
         endif
    end select

    if(ginv .ne. ginv) then 
           write(*,*) 'ginv, w_pt_local(psi_metric_)'
           write(*,*) ginv, w_pt_local(psi_metric_)
         stop  "and yes, it's also a NaN ginv sometimes" 
    endif


   call mpistop('passed gammainv_componet analytic')

  end subroutine get_gammainv_component_analytic_cfc_sp



!=============================================================================
!  Below are flat metric
!=============================================================================

  subroutine get_sqrtgamma_analytic_sp(x^D,sqrtgamma,is_analytic,w_pt)

    include 'amrvacdef.f'

    ! Since the metric is numeric, sqrtgamma is not calculated analytically.

    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    double precision, intent(in), optional           :: w_pt(1:nw)


    if(present(is_analytic)) is_analytic = .true.

    {^IFZIN
    sqrtgamma = x1**2 * abs(sin(x^Z))
    }{^IFZOUT
    sqrtgamma = x1**2
    }

  end subroutine get_sqrtgamma_analytic_sp

!=============================================================================
  subroutine get_alpha_sp(x^D,alpha,iszero,dalphadj_iszero,dalphadj,jdir,w_pt)

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
         call mpistop('cfc not calucating metric dervatives here! here is get_alpha_sp')
    endif

    alpha = 1.0d0

  end subroutine get_alpha_sp
  !=============================================================================
  subroutine get_beta_sp(idir,x^D,beta,iszero,dbetaidj_iszero,dbetaidj,jdir,w_pt)

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
         call mpistop('cfc not calucating metric dervatives here! here is get_beta_sp')
    endif
    if (present(dbetaidj_iszero)) &
         dbetaidj_iszero = .true. !.false.

  end subroutine get_beta_sp
  !=============================================================================
  subroutine get_g_component_sp(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir,w_pt)

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

    ! Diagonal for now
    if (i==j) then
       if(present(iszero)) &
            iszero = .false.

       if (i .eq. 1) then ! r direction
          g = 1.0d0
          if (present(kdir)) then
             if (kdir .eq.1) then
                if(present(dgdk_iszero)) dgdk_iszero = .true. !.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
             else
                if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
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
                !if (present(dgdk)) dgdk = 2.0d0*x1*sin(x^Z)**2
                if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
                }{^IFZOUT
                !if (present(dgdk)) dgdk = 2.0d0*x1
                if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
                }
             else if (kdir .eq. ^Z) then
                {^IFZIN
                !if (present(dgdk)) dgdk = 2.0d0*x1**2*cos(x^Z)*sin(x^Z)
                if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
                }{^IFZOUT
                !if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
                }
             else
                !if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true. !.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
             end if
          end if

       else if (i .eq. ^Z) then ! theta direction
          g = x1**2

          if (present(kdir)) then
             if (kdir .eq.1) then
                !if (present(dgdk)) dgdk = 2.0d0*x1
                if (present(dgdk_iszero)) dgdk_iszero = .true.!.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
             else
                !if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.!.false.
                if (present(dgdk)) then
                    call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
                endif
             end if
          end if

       end if
    else
       g = 0.0d0

       if(present(iszero)) &
            iszero = .true.

       !dgdk = 0.0d0
       if (present(dgdk)) then
           call mpistop('cfc not calucating metric dervatives here! here is get_g_com_sp')
       endif

       if(present(dgdk_iszero)) &
            dgdk_iszero = .true.

    end if


  end subroutine get_g_component_sp



  subroutine BLToCoord(ixI^L,ixO^L,xBL,xCoord)

    ! identity
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    integer                                                :: ix^D
!    call mpistop("called BLtoCoord, cfc do not need coordinate trans")
    xCoord = xBL


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

!    call mpistop("called CoordtoBL, cfc do not need coordinate trans")
    xBL = xCoord


  end subroutine CoordToBL
  !=============================================================================
  subroutine CoordToKS(ixI^L,ixO^L,xCoord,xKS)

    ! Of course these are not really KS coordinates,
    ! this routine only converts to non-logarithmic coordinates
    ! but is useful when the analogous routine is called for black hole
    ! spacetimes
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xKS
    ! .. local ..
    integer                                                :: ix^D

!    call mpistop("called CoordToKS, cfc do not need coordinate trans")
    xKS = xCoord


  end subroutine CoordToKS
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

!    call mpistop("called u4BLtoCoord, cfc do not need coordinate trans")
!    if (eqpar(a_) .ne. zero) &
!         call mpistop("u4BLtoCoord: Spin not zero, transformation to Schwarzschild makes no sense.")

    u4Coord(ixO^S,0:^NC) = u4BL(ixO^S,0:^NC)


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
! call mpistop('called coordtocart')

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

     call CoordToBL(ixI^L,ixO^L,x,xBL)

    {^IFPHIIN
    sph(ixO^S) = sin(xBL(ixO^S,phi_))
    cph(ixO^S) = cos(xBL(ixO^S,phi_))
    \}{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    \}

    {^IFZIN
    sth(ixO^S) = sin(xBL(ixO^S,z_))
    cth(ixO^S) = cos(xBL(ixO^S,z_))
    \}{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    \}

    r(ixO^S) = xBL(ixO^S,1)

    xCart(ixO^S,1) = r(ixO^S) * sth(ixO^S)*cph(ixO^S)
    if (ycart_.ne.1) &
         xCart(ixO^S,ycart_) = r(ixO^S) * sth(ixO^S)*sph(ixO^S)
    if (zcart_.ne.1) &
         xCart(ixO^S,zcart_) = r(ixO^S)*cth(ixO^S)



!    call mpistop("called CoordTocart, cfc do not need coordinate trans")
!   if (ndim.eq.3) then
!       ycart_=2
!       zcart_=3
!    else if (ndim.eq.2 .and. ^Z.eq.2) then
!       ycart_=1 ! must catch ycart_=1 case and don't use
!       zcart_=2
!    else if (ndim.eq.2 .and. ^PHI.eq.2) then
!       ycart_=2
!       zcart_=1 ! must catch zcart_=1 case and don't use
!    else if (ndim.eq.1) then
!       ycart_=1
!       zcart_=1
!    else
!       call mpistop("CoordToCart: unknown parameter combination of -d, -phi and -z!")
!    end if
!
!    {^IFPHIIN
!    sph(ixO^S) = sin(x(ixO^S,phi_))
!    cph(ixO^S) = cos(x(ixO^S,phi_))
!    }{^IFPHIOUT
!    sph(ixO^S) = zero
!    cph(ixO^S) = one
!    }
!
!    {^IFZIN
!    sth(ixO^S) = sin(x(ixO^S,z_))
!    cth(ixO^S) = cos(x(ixO^S,z_))
!    }{^IFZOUT
!    sth(ixO^S) = one
!    cth(ixO^S) = zero
!    }
!
!    ra(ixO^S) = x(ixO^S,1)!sqrt(x(ixO^S,1)**2+eqpar(a_)**2)
!
!
!    xCart(ixO^S,1) = ra(ixO^S) * sth(ixO^S)*cph(ixO^S)
!    if (ycart_.ne.1) &
!         xCart(ixO^S,ycart_) = ra(ixO^S) * sth(ixO^S)*sph(ixO^S)
!    if (zcart_.ne.1) &
!         xCart(ixO^S,zcart_) = x(ixO^S,1)*cth(ixO^S)

  end subroutine CoordToCart


  subroutine u4CoordToCart(ixI^L,ixO^L,x,u4Coord,u4Cart,J)

    ! Transforms any contravariant four vector u4Coord in the current coordinates
    ! to a vector in the basis of "Boyer-Lindquist Cartesian coordinates".
    ! t = tBL
    ! x = sqrt(rBL**2+a**2) Sin(thetaBL) Cos(phiBL)
    ! y = sqrt(rBL**2+a**2) Sin(thetaBL) Sin(phiBL)
    ! z = rBL Cos(thetaBL)
    ! We are already in BL-cordinates, so we only have to do one transformation.
    ! The transformation matrix J=del(t,x,y,z)/del(tBL,rBL,thetaBL,phiBL)
    !
    !     | 1               0                             0                        0                          |
    !     | 0   rBL/ra Sin(thetaBL) Cos(phiBl)    ra Cos(thetaBL) Cos(phiBL)    - ra Sin(thetaBL) Sin(phiBL)  |
    ! J = | 0   rBL/ra Sin(thetaBL) Sin(phiBL)    ra Cos(thetaBL) Sin(phiBL)      ra Sin(thetaBL) Cos(phiBL)  |
    !     | 0   Cos(thetaBL)                    - rBL Sin(thetaBL)                 0                          |
    !
    ! where ra=sqrt(rBL**2+a**2) ! WARNING: I decided to put ra=rBL for now...
    !
    ! u4Cart = J u4Coord
    !
    ! 3D simulation, three vector components, (-d=33 -z=2, phi=3 , the standard):
    ! u4Coord = (utBL,urBL,uthetaBL,uphiBL)
    ! u4Cart  = (ut,ux,uy,uz)
    !
    ! r-theta plane (phi=0), three vector components (-d=23 -z=2 -phi=3):
    ! u4Coord = (utBL,urBL,uthetaBL,uphiBL)
    ! u4Cart = (ut,ux,uy,uz)
    !
    ! r-theta plane (phi = 0), two vector components (-d=22 -z=2 -phi=0):
    ! u4Coord = (utBL,urBL,uthetaBL)
    ! u4Cart = (ut,ux,uz)
    !
    ! r-phi plane (theta = pi/2), three vector components (-d=23 -phi=2 -z=3)
    ! u4Coord = (utBL,urBL,uphiBL,uthetaBL)
    ! u4Cart = (ut,ux,uy,uz)
    !
    ! r-phi plane (theta = pi/2), two vector components (-d=22 -phi=2 -z=0)
    ! u4Coord = (utBL,urBL,uphiBL)
    ! uCart = (ut,ux,uy)

    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4Coord
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Cart
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixI^S,0:^NC,0:^NC)         :: Jac
    double precision, dimension(ixI^S,1:^ND)               :: xBL
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, ra
    integer                                                :: zcart_, ycart_, ix,jx

!    call mpistop("called u4CoordTocart, cfc do not need coordinate trans")
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
       ycart_=0 ! must catch ycart_=0 case and don't use
       zcart_=2
    else if (ndir.eq.2 .and. ^PHI.eq.2) then
       ycart_=2
       zcart_=0 ! must catch zcart_=0 case and don't use
    else if (ndir.eq.1) then
       ycart_=0
       zcart_=0
    else
       call mpistop("u4CoordToCart: unknown parameter combination of -d, -phi and -z!")
    end if

    call CoordToBL(ixI^L,ixO^L,x,xBL)

    {^IFPHIIN
    sph(ixO^S) = sin(xBL(ixO^S,phi_))
    cph(ixO^S) = cos(xBL(ixO^S,phi_))
    \}{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    \}

    {^IFZIN
    sth(ixO^S) = sin(xBL(ixO^S,z_))
    cth(ixO^S) = cos(xBL(ixO^S,z_))
    \}{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    \}


    Jac(ixO^S,:,:) = zero
    Jac(ixO^S,0,0) = one

    Jac(ixO^S,1,1) = sth(ixO^S)*cph(ixO^S)

    {^IFZ
    Jac(ixO^S,1,^Z) = xBL(ixO^S,1) * cth(ixO^S)*cph(ixO^S)
    \}
    {^IFPHI
    Jac(ixO^S,1,^PHI) = - xBL(ixO^S,1)*sth(ixO^S)*sph(ixO^S)
    \}

    if (ycart_.ne.1) then

       Jac(ixO^S,ycart_,1) = sth(ixO^S)*sph(ixO^S)

       {^IFZ
       Jac(ixO^S,ycart_,^ZZ) = xBL(ixO^S,1)*cth(ixO^S)*sph(ixO^S)
       \}
       {^IFPHI
       Jac(ixO^S,ycart_,^PPHI) = xBL(ixO^S,1)*sth(ixO^S)*cph(ixO^S)
       \}
    end if

    if(zcart_.ne.1) then
       Jac(ixO^S,zcart_,1) = cth(ixO^S)

       {^IFZ
       Jac(ixO^S,zcart_,^ZZ) = - xBL(ixO^S,1)*sth(ixO^S)
       \}
       {^IFPHI
       Jac(ixO^S,zcart_,^PPHI) = 0.0d0
       \}
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled, now
    ! transform contravariant four-vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    u4Cart(ixO^S,:) = zero
    do ix=0,ndir
       do jx=0,ndir
          u4Cart(ixO^S,ix) = u4Cart(ixO^S,ix) + Jac(ixO^S,ix,jx) * u4Coord(ixO^S,jx)
       end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (present(J)) J = Jac



  end subroutine u4CoordToCart
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

! call mpistop('u3Coordtocart')

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
    \}{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    \}

    {^IFZIN
    sth(ixO^S) = sin(xBL(ixO^S,z_))
    cth(ixO^S) = cos(xBL(ixO^S,z_))
    \}{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    \}

    Jac(ixO^S,1,1) = sth(ixO^S)*cph(ixO^S)

    {^IFZ
    Jac(ixO^S,1,^Z) = xBL(ixO^S,1) * cth(ixO^S)*cph(ixO^S)
    \}
    {^IFPHI
    Jac(ixO^S,1,^PHI) = - xBL(ixO^S,1)*sth(ixO^S)*sph(ixO^S)
    \}

    if (ycart_.ne.1) then

       Jac(ixO^S,ycart_,1) = sth(ixO^S)*sph(ixO^S)

       {^IFZ
       Jac(ixO^S,ycart_,^ZZ) = xBL(ixO^S,1)*cth(ixO^S)*sph(ixO^S)
       \}
       {^IFPHI
       Jac(ixO^S,ycart_,^PPHI) = xBL(ixO^S,1)*sth(ixO^S)*cph(ixO^S)
       \}
    end if

    if(zcart_.ne.1) then
       Jac(ixO^S,zcart_,1) = cth(ixO^S)

       {^IFZ
       Jac(ixO^S,zcart_,^ZZ) = - xBL(ixO^S,1)*sth(ixO^S)
       \}
       {^IFPHI
       Jac(ixO^S,zcart_,^PPHI) = 0.0d0
       \}
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



!    if (ndim.eq.3.and.ndir.eq.3) then
!       ycart_=2
!       zcart_=3
!    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. ^PHI .eq. 2) then
!       ycart_=2
!       zcart_=3
!    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. ^Z .eq. 2) then
!       ycart_=3
!       zcart_=2
!    else if (ndir.eq.2 .and. ^Z.eq.2) then
!       ycart_=1 ! must catch ycart_=1 case and don't use
!       zcart_=2
!    else if (ndir.eq.2 .and. ^PHI.eq.2) then
!       ycart_=2
!       zcart_=1 ! must catch zcart_=1 case and don't use
!    else if (ndir.eq.1) then
!       ycart_=1
!       zcart_=1
!    else
!       call mpistop("u3CoordToCart: unknown parameter combination of -d, -phi and -z!")
!    end if
!
!
!    {^IFPHIIN
!    sph(ixO^S) = sin(x(ixO^S,phi_))
!    cph(ixO^S) = cos(x(ixO^S,phi_))
!    }{^IFPHIOUT
!    sph(ixO^S) = zero
!    cph(ixO^S) = one
!    }
!
!    {^IFZIN
!    sth(ixO^S) = sin(x(ixO^S,z_))
!    cth(ixO^S) = cos(x(ixO^S,z_))
!    }{^IFZOUT
!    sth(ixO^S) = one
!    cth(ixO^S) = zero
!    }
!
!    ra(ixO^S) = x(ixO^S,1)
!
!    Jac(ixO^S,1,1) = x(ixO^S,1)/ra(ixO^S) * sth(ixO^S)*cph(ixO^S)
!    {^IFZ
!    Jac(ixO^S,1,^Z) = ra(ixO^S) * cth(ixO^S)*cph(ixO^S)
!    }
!    {^IFPHI
!    Jac(ixO^S,1,^PHI) = - ra(ixO^S)*sth(ixO^S)*sph(ixO^S)
!    }
!
!    if (ycart_.ne.1) then
!       Jac(ixO^S,ycart_,1) = x(ixO^S,1)/ra(ixO^S) * sth(ixO^S)*sph(ixO^S)
!       {^IFZ
!       Jac(ixO^S,ycart_,^ZZ) = ra(ixO^S) * cth(ixO^S)*sph(ixO^S)
!       }
!       {^IFPHI
!       Jac(ixO^S,ycart_,^PPHI) = ra(ixO^S)*sth(ixO^S)*cph(ixO^S)
!       }
!    end if
!
!    if(zcart_.ne.1) then
!       Jac(ixO^S,zcart_,1) = cth(ixO^S)
!       {^IFZ
!       Jac(ixO^S,zcart_,^ZZ) = - x(ixO^S,1)*sth(ixO^S)
!       }
!       {^IFPHI
!       Jac(ixO^S,zcart_,^PPHI) = 0.0d0
!       }
!    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! Jacobian fully assembled, now
!    ! transform contravariant three-vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    u3Cart(ixO^S,:) = zero
!    do ix=1,ndir
!       do jx=1,ndir
!          u3Cart(ixO^S,ix) = u3Cart(ixO^S,ix) + Jac(ixO^S,ix,jx) * u3Coord(ixO^S,jx)
!       end do
!    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    if (present(J)) J = Jac

!    call mpistop("called u3CoordTocart, cfc do not need coordinate trans")

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

end module mod_coord_cfc_sp
