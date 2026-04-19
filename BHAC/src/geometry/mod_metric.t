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
!
! Metric and geometry for general covariant coordinates.
! To use, set the 'typecoord' keyword to a string containting 'covariant'.  
! The metric is initialised by providing values for alpha, shift-beta and 
! the spatial part of the metric.
! 
! For the dynamics, the metric is assumed to be of the form
! ds^2 = (-\alpha^2+\beta_i \beta^i)dt^2 + 2 \beta_i dt dx^i + \gamma_{ij} dx^i dx^j
! where \beta_i = \gamma_{ij}beta^j.  
! See e.g. the book of Alcubierre, Introduction to 3+1 Numerical Relativity (2008)
! Chapter 2.2
!
! The Valencia forumation is used for the evolution with first class citizens:
! Spatial metric \gamma_{ij}
! Lapse function \alpha
! Contravariant shift vector \beta^i
! and their spatial partial derivatives.  
! 
! It is also possible to initialise the four-metric first and the code will
! automatically obtain lapse alpha and shift beta which are needed for the
! dynamical evolution.
! To do so, set init_from_g4 = .true. in your corresponding mod_coord_ file.
! 
! In the rest of the code, the metric is available via the pointer myM
! or via the geometry datastructure.
! However, if you define DY_SP for BSSN or xCFC dynamical spacetime, 
! we should not use myM
! 
! For the definition of the datastructure see 'mod_physicaldata.t'.
! 2014-03-23 by Oliver Porth
! 
!=============================================================================
module mod_metric
  use mod_metric_aux
  use mod_cfc_parameters
{#IFDEF DY_SP
!  use mod_coord_cfc_sp
!  use mod_coord_cfc_cart
  use mod_coord_cfc
}
  implicit none

  integer, save                        :: nnonzero_metric ! number of non-zero elements in four-metric
  integer, save                        :: nnonzero_beta   ! number of non-zero elements in contravariant shift
  integer, save                        :: nnonzero_dalphadj   ! number of non-zero elements in lapse-derivatives  
  integer, save                        :: nnonzero_dbetaidj   ! number of non-zero elements in shift-derivatives  
  integer, save                        :: nnonzero_dgdk   ! number of non-zero elements in metric-derivatives
  logical, save,dimension(0:^NC,0:^NC) :: g_is_zero       ! four-metric
  logical, save,dimension(1:^NC)       :: beta_is_zero    ! contravariant part, not in metric
  logical, save,dimension(1:^NC,1:^NC,1:^NC) :: dgdk_is_zero       ! partial derivatives of three-metric
  logical, save,dimension(1:^NC)       :: dalphadj_is_zero    ! partial derivatives of lapse
  logical, save,dimension(1:^NC,1:^NC) :: dbetaidj_is_zero    ! partial derivatives of shift 
  logical, save                        :: space_is_diagonal, sqrtgamma_is_analytic


{#IFNDEF DY_SP
INCLUDE:geometry/mod_coord_^COORD.t
}

{#IFDEF DY_SP
contains
}
!=============================================================================
subroutine init_coord_check

  if (.not. associated(get_gammainv_component_analytic)) &
       get_gammainv_component_analytic => dummy_get_gammainv_component_analytic
  if (.not. associated(get_gammainv_component)) &
       get_gammainv_component => get_gammainv_component_numerical
  if (.not. associated(SPToCoord)) &
       SPToCoord => dummy_SPToCoord
  if (.not. associated(LRNFu)) &
       LRNFu => dummy_LRNFu ! locally non-rotation orthonormal tetrad transform
  if (.not. associated(u4CoordToKS)) &
       u4CoordToKS => dummy_u4CoordToKS ! Transform four-vector to KS coordinates
  if (.not. associated(u4CoordToBL)) &
       u4CoordToBL => dummy_u4CoordToBL ! Transform four-vector to BL coordinates
  if (.not. associated(d4CoordToKS)) &
       d4CoordToKS => dummy_d4CoordToKS ! Transform four-covariant vector to KS coordinates

end subroutine init_coord_check
!=============================================================================
subroutine get_g4(x^D,g,dgdk)
  ! Return the four-metric at point x^D
  ! This is just a convenience routine and fairly slow.
  ! Can also return metric derivatives when dgdk given.  
  ! Avoid extensive use, expecially when init_from_g4=.false.
  double precision, intent(in)                           :: x^D
  double precision, dimension(0:^NC,0:^NC),intent(out)   :: g
  double precision, dimension(0:^NC,0:^NC,1:^NC), optional, intent(out) :: dgdk
  ! .. local ..
  integer                                                :: i, j, k
  double precision, dimension(1:^NC)                     :: betaU, betaD, dalphadj
  double precision, dimension(1:^NC,1:^NC)               :: dbetaUdj, dbetaDdj
  double precision                                       :: beta2, alpha, dummy
  !-----------------------------------------------------------------------------

  if (.not. init_from_g4) then

     do j=1,^NC
        do i=1,^NC
           if (present(dgdk)) then
              do k = 1, ^NC
                 call get_g_component(i,j,x^D,g(i,j),dgdk=dgdk(i,j,k),kdir=k)
              end do
           else
              call get_g_component(i,j,x^D,g(i,j))
           end if
        end do
     end do

     do i=1,^NC
        call get_alpha(x^D,alpha,dalphadj=dalphadj(i),jdir=i)
        if (present(dgdk)) then
           do j=1,^NC
              call get_beta(i,x^D,betaU(i),dbetaidj=dbetaUdj(i,j),jdir=j)
           end do
        else
           call get_beta(i,x^D,betaU(i))   
        end if
     end do
     ! lower the beta:
     betaD(:) = 0.0d0
     do j=1,^NC
        do i=1,^NC
           if (g_is_zero(i,j) .or. beta_is_zero(j)) cycle
           betaD(i) = betaD(i) + g(i,j)*betaU(j)
        end do
     end do

     if (present(dgdk)) then
        do j=1,^NC
           do i=1,^NC
              dbetaDdj(i,j) = {^C& betaU(^C) * dgdk(i,^C,j) |+} & 
                   + {^C& g(i,^C)*dbetaUdj(^C,j) |+}
           end do
        end do
     end if

     beta2 = {^C& betaU(^C)*betaD(^C) |+}

     g(0,0) = -alpha**2 + beta2

     do i=1,^NC
        g(i,0) = betaD(i)
        g(0,i) = betaD(i)
     end do

     if (present(dgdk)) then
        do i=1,^NC
	   dgdk(0,0,i) = -2.0d0*alpha*dalphadj(i) + {^C& betaD(^C)*dbetaUdj(^C,i) |+} & 
                + {^C& betaU(^C)*dbetaDdj(^C,i) |+}
	end do
        do i=1,^NC
           do j=1,^NC
              dgdk(0,i,j) = dbetaDdj(i,j)
              dgdk(i,0,j) = dbetaDdj(i,j)
           end do
        end do
     end if

  else ! init_from_g4

     do j=0,^NC
        do i=0,^NC
           if (present(dgdk)) then
              do k = 1, ^NC
                 call get_g_component(i,j,x^D,g(i,j),dgdk=dgdk(i,j,k),kdir=k)
              end do
           else
              call get_g_component(i,j,x^D,g(i,j))
           end if
        end do
     end do

  end if ! init_from_g4

end subroutine get_g4
!=============================================================================
subroutine g4inv(ixI^L,ixO^L,m,ginv)
  include 'amrvacdef.f'
  ! Calculates the full inverse 4-metric using 3+1 quantities in m
  integer, intent(in)                      :: ixI^L, ixO^L
  type(metric), pointer                    :: m
  double precision, intent(out), dimension(ixI^S,0:ndir,0:ndir)  :: ginv
  ! .. local ..
  integer                                  :: i,j
  !-----------------------------------------------------------------------------

  ginv(ixO^S,0,0) = -1.0d0/m%alpha(ixO^S)**2
  do j=1,^NC
     do i=0,j
        if (i .eq. 0) then
           ginv(ixO^S,0,j) = m%beta(j)%elem(ixO^S)/m%alpha(ixO^S)**2
        else
           ginv(ixO^S,i,j) = m%gammainv(i,j)%elem(ixO^S) &
                - m%beta(i)%elem(ixO^S)*m%beta(j)%elem(ixO^S) &
                / m%alpha(ixO^S)**2
        end if
        if (i .ne. j) ginv(ixO^S,j,i) = ginv(ixO^S,i,j)
     end do
  end do
  
  
end subroutine g4inv
!=============================================================================
subroutine get_sqrtgamma(x^D,sqrtgamma,w_pt)
    include 'amrvacdef.f'

    ! Calculate the determinant of the spatial metric at point x
    ! Takes analytic values if provided, otherwise bruteforcing...
    double precision, intent(in)                            :: x^D
    double precision, intent(out)                           :: sqrtgamma
    double precision, dimension(1:^NC,1:^NC)                :: g
    double precision, intent(in), optional                  :: w_pt(1:nw)
    double precision                                        :: w_pt_local(1:nw)

    ! .. local ..
    integer                                                 :: i, j
    !-----------------------------------------------------------------------------
   
    if (present(w_pt)) w_pt_local = w_pt
 
    if (sqrtgamma_is_analytic) then
       {#IFDEF DY_SP
       call get_sqrtgamma_analytic(x^D,sqrtgamma,w_pt=w_pt_local(1:nw))
       }
       {#IFNDEF DY_SP
       call get_sqrtgamma_analytic(x^D,sqrtgamma)
       }
       
    else

       {#IFDEF DY_SP
           call mpistop('sqrtgamma must be analytic in cfc')
       }

       do j=1,^NC
          do i=1,j
             call get_g_component(i,j,x^D,g(i,j))
             if (i .ne. j) g(j,i) = g(i,j)
          end do
       end do

       if (space_is_diagonal) then 
          ! If space is diagonal, just multiply diagonal entries
          sqrtgamma = 1.0d0
          do j=1,^NC
              sqrtgamma = sqrtgamma*g(j,j)
          end do
       else
          ! Else, compute the determinant by LU decomposition
          sqrtgamma = determinant(g,^NC)
       end if

       ! Check that det(gamma) is positive and
       ! send warning message if not.
       ! HO: I commented this since it spams the standard output
!       if (sqrtgamma.lt.0.0d0) &
!          print *, 'WARNING: metric determinant is negative at x^D =',&
!                    x^D,'taking absolute value'

       sqrtgamma = sqrt(abs(sqrtgamma))
    end if

  end subroutine get_sqrtgamma
  !=============================================================================
  subroutine get_gammainv_component_numerical(iin,jin,x^D,ginv,iszero,dginvdk_iszero,dginvdk,kdir,w_pt)

    ! Obtains a component of the inverse spatial metric gamma
    ! Not particularly efficient.
    
    use mod_lu
    include 'amrvacdef.f'

    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: ginv
    logical, optional, intent(out)           :: iszero, dginvdk_iszero
    double precision, optional, intent(out)  :: dginvdk
    double precision, optional, intent(in)   :: w_pt(1:nw)

    ! .. local ..
    double precision                         :: w_pt_local(1:nw)
    integer, dimension(1:^NC)                :: indx
    integer                                  :: i, j, code, d
    double precision, dimension(1:^NC,1:^NC) :: g
    double precision, dimension(1:^NC)       :: b
    !-----------------------------------------------------------------------------

    if (present(w_pt)) w_pt_local = w_pt


    if (space_is_diagonal) then 

       if (iin .eq. jin) then
          {#IFDEF DY_SP
          call get_g_component(iin,iin,x^D,g(iin,iin),w_pt=w_pt_local(1:nw))
          }
          {#IFNDEF DY_SP
          call get_g_component(iin,iin,x^D,g(iin,iin))
          }
          ginv = 1.0d0/g(iin,iin)
       else
          ginv = 0.0d0
       end if

    else ! space is non diagonal

       {#IFDEF DY_SP
         call mpistop('cfc space is diagonal')
       }

       do j=1,^NC
          do i=1,j
             call get_g_component(i,j,x^D,g(i,j))
             if ( i .ne. j) g(j,i) = g(i,j)
          end do
       end do

       call LUDCMP(g,^NC+1,indx,d,code)

       {b(^C)=kr(jin,^C)|;}

       call LUBKSB(g,^NC,indx,b)

       ginv = b(iin)

    end if 

   ! !!!! NaN catcher
   ! if (ginv.ne.ginv) then
   !    ginv = 0.0d0
   !    if (iin.eq.jin) ginv = 1.0d0
   ! end if

    if (present(iszero)) iszero = .false. ! At least we cant tell with certainty
    if (present(dginvdk_iszero)) dginvdk_iszero = .false. ! At least we cant tell with certainty
    if (present(dginvdk)) call mpistop('get_gammainv_component_numerical: dginvdk not yet implemented.')
    
  end subroutine get_gammainv_component_numerical
  !=============================================================================
  double precision function determinant(g,n)
    ! determinant for up to nxn matrices
    use mod_lu
    integer, intent(in) :: n
    double precision, dimension(1:n,1:n), intent(in)   :: g
    ! .. local ..
    double precision :: det
    integer          :: idx
    integer          :: code,indx(3),d
    !-----------------------------------------------------------------------------

    ! Compute the determinant as the product of the diagonal
    ! entries of U in the LU decompostion
    ! (the diagonal of L just contains ones).
    ! The integer d=+-1 keeps track of the permutations and
    ! changes the sign accortingly.

    call ludcmp_d(g,n,indx,d,code)

    determinant = d
    do idx=1,n
       determinant = determinant*g(idx,idx)
    end do

  end function determinant
  !=============================================================================
  subroutine get_sqrtgammai(idims,x^D,s,w_pt)
    include 'amrvacdef.f'
    ! Calculate the determinant of the induced metric in direction idims
    
    integer, intent(in)                                 :: idims
    double precision,intent(in)                         :: x^D
    double precision, intent(out)                       :: s
    double precision, optional, intent(in)              :: w_pt(1:nw)

    ! .. local ..
    integer                                             :: ii,jj,i,j
    double precision                                    :: w_pt_local(1:nw)
{#IFNDEF C1
    double precision, dimension(1:^NC-1,1:^NC-1)        :: gi
}
    !-----------------------------------------------------------------------------

    if (present(w_pt)) w_pt_local = w_pt

{#IFDEF C1
s = 1.0d0
}

{#IFNDEF C1
    i=0
    do ii=1,^NC
       if (ii.eq.idims) cycle
       i = i+1
       j = 0
       do jj=1,^NC
          if (jj.eq.idims) cycle
          j = j+1
          {#IFDEF DY_SP
          call get_g_component(ii,jj,x^D,gi(i,j),w_pt=w_pt_local(1:nw))
          }
          {#IFNDEF DY_SP
          call get_g_component(ii,jj,x^D,gi(i,j))
          }
       end do
    end do

{#IFDEF C2
    s = sqrt(gi(1,1))}
{#IFDEF C3
    s = sqrt(gi(1,1)*gi(2,2) - gi(1,2)*gi(2,1))}
}
  end subroutine get_sqrtgammai
  !=============================================================================
  subroutine init_metric()
    include 'amrvacdef.f'
    
    ! .. local ..
    integer                            :: i,j,k,inonzero
    double precision                   :: dummy
    !-----------------------------------------------------------------------------
    {^IFCOORDNUM
    call read_metric('numerical.met')}
    {^IFCOORDNUMCTOB 
    call read_metric('numerical.met')} 


    ! Initialize the auxiliary subroutines (better strategy):
    call init_coord
    call init_coord_check
    
    ! Check if we want to use analytic expression for determinant:
    call get_sqrtgamma_analytic(xprobmin^D,dummy,is_analytic=sqrtgamma_is_analytic)
    
    ! Check which components from the spatial part are identical zero:
    do j=1,^NC
       do i=1,^NC
          call get_g_component(i,j,xprobmin^D,dummy,iszero=g_is_zero(i,j))
       end do
    end do

    ! Time component 00 cannot be zero:
    g_is_zero(0,0) = .false.

    ! Contravariant shifts:
    inonzero = 0
    do i=1,^NC
       call get_beta(i,xprobmin^D,dummy,beta_is_zero(i))
       if (.not.beta_is_zero(i)) inonzero = inonzero+1
    end do
    nnonzero_beta = inonzero

    ! Time-space components elements of metric
    ! Check if lowering beta would still result in zero:
    do i=1,^NC
       g_is_zero(0,i) = .true.
       do j=1,^NC
          if ( .not.(g_is_zero(i,j) .or. beta_is_zero(j)) ) then
             g_is_zero(0,i) = .false.
          end if
       end do
       g_is_zero(i,0) = g_is_zero(0,i)
    end do
       

    ! Check if space is diagonal and set the flag
    ! also count the number of non-zero elements
    space_is_diagonal = .true.
    inonzero = 0
    do j=0,^NC
       do i=0,^NC
          if (.not.g_is_zero(i,j)) inonzero = inonzero + 1
          if (i.eq.j .or. i.eq.0 .or. j.eq.0) cycle
          space_is_diagonal = space_is_diagonal .and. g_is_zero(i,j)
       end do
    end do
    nnonzero_metric = inonzero
    
    ! Derivatives of lapse:
    inonzero = 0
    do j=1,^NC
       call get_alpha(xprobmin^D,dummy,dalphadj_iszero=dalphadj_is_zero(j),jdir=j)
       if (.not. dalphadj_is_zero(j)) inonzero = inonzero + 1
    end do
    nnonzero_dalphadj = inonzero
    
    ! Derivatives of shift:
    inonzero = 0
    do j=1,^NC
       do i=1,^NC
          call get_beta(i,xprobmin^D,dummy,dbetaidj_iszero=dbetaidj_is_zero(i,j),jdir=j)
          if (.not. dbetaidj_is_zero(i,j)) inonzero = inonzero + 1
       end do
    end do
    nnonzero_dbetaidj = inonzero
    
    ! Derivatives of three-metric:
    inonzero = 0
    do k=1,^NC
       do j=1,^NC
          do i=1,^NC
             call get_g_component(i,j,xprobmin^D,dummy,dgdk_iszero=dgdk_is_zero(i,j,k),kdir=k)
             if (.not. dgdk_is_zero(i,j,k)) inonzero = inonzero + 1
          end do
       end do
    end do
    nnonzero_dgdk = inonzero
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Done collecting information on spacetime !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (mype.eq.0) print*,'-----------------------------------------------------------------------------'
    if (mype.eq.0) write(*,*) 'Initialised Meta-info on metric:'
    if (mype.eq.0) print*,'-----------------------------------------'
    if (mype.eq.0) write(*,*) 'Using ', trim(coord), ' coordinates'
    if (mype.eq.0 .and. space_is_diagonal) write(*,*) 'Space is diagonal'
    if (mype.eq.0 .and. .not.space_is_diagonal) write(*,*) 'Space is non-diagonal'
    if (mype.eq.0) write(*,*) 'Sqrt(gamma) is analytic:', sqrtgamma_is_analytic
    ! The metric will only be filled as a ^NC x ^NC matrix
    ! For the induced metrics (surfaces), it is assumed that the missing entries
    ! are diagonal 1.  If this is not what you have in mind, you should run with
    ! 3 components.  
    if (mype.eq.0.and.^NC.ne.3) write(*,*) 'Assuming left-out components are identity'

    ! print the shape of the metric to screen:
    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These metric elements are zero:'
       do i=0,^NC
          write(*,*) g_is_zero(i,:)
       end do
       write(*,*) 'Number of non-zero elements:', nnonzero_metric
    end if
    
    ! print the shape of the contravariant shift to screen:
    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These shift-components are zero:'
       write(*,*) beta_is_zero(:)
       write(*,*) 'Number of non-zero elements:', nnonzero_beta
    end if

    ! print the shape of the metric derivatives to screen:
    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These components of the metric derivative DgijDk are zero:'
       do k=1,^NC
       print*,'-----------'
       write(*,*) 'k=',k
          do i=1,^NC
             write(*,*) dgdk_is_zero(i,:,k)
          end do
       end do
       print*,'-----------'
       write(*,*) 'Number of non-zero elements:', nnonzero_dgdk
    end if

    ! print the shape of the shift derivatives:
    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These elements of shift derivatives dbetaidj are zero:'
       do i=1,^NC
          write(*,*) dbetaidj_is_zero(i,:)
       end do
       write(*,*) 'Number of non-zero elements:', nnonzero_dbetaidj
    end if

    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These elements of the lapse derivative dalphadj are zero:'
       write(*,*) dalphadj_is_zero(:)
       write(*,*) 'Number of non-zero elements:', nnonzero_dalphadj
    end if

    if (mype.eq.0) print*,'-----------------------------------------------------------------------------'
    
  end subroutine init_metric
  !=============================================================================
  subroutine LU(m)
    ! LU-decomposes the spatial metric and calculates the inverse
    ! much room for optimization. 
    use mod_lu
    ! if this breaks for you, try commenting out and use the old gfortran hack below.
    use, intrinsic :: IEEE_ARITHMETIC, only: ieee_value, ieee_quiet_nan
    !
    include 'amrvacdef.f'

    type(metric), pointer                    :: m
    ! .. local ..
    integer, dimension(^NC)                  :: indx
    integer                                  :: inonzero, i, j, code, d, ix^D
    double precision, dimension(1:^NC,1:^NC) :: a
    double precision, dimension(1:^NC)       :: b
! old gfortran (< v5) needed this hack but it breaks on gfortran 10:    
!    double precision, PARAMETER :: D_QNAN = TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
    double precision                         :: D_QNAN
    D_QNAN = ieee_value(1.,ieee_quiet_nan)
    !-----------------------------------------------------------------------------
    
    if (space_is_diagonal) then 

       do i=1,^NC
          m%gammainv(i,i)%elem(m%ixG^S) = 1.0d0/m%g(i,i)%elem(m%ixG^S)
       end do

    else ! space is diagonal

       
       ! --------------------------------------------------
       ! Initialize as quiet NaN
       ! --------------------------------------------------
       do j=1,^NC
          do i=1,j
             m%gammainv(i,j)%elem(m%ixG^S) = D_QNAN
          end do
       end do
       

       ! --------------------------------------------------
       ! re-arrange elements and do the LU by using LUDCMP
       ! from Numerical recipes
       ! --------------------------------------------------
       {^D& do ix^DB=m%ixGmin^DB, m%ixGmax^DB\}
       a(:,:) = 0.0d0
       do inonzero=1,m%nnonzero
          i = m%nonzero(inonzero)%i
          j = m%nonzero(inonzero)%j
          if (i==0 .or. j==0) cycle

          a(i,j) = m%nonzero(inonzero)%elem(ix^D)

       end do

       call LUDCMP(a,^NC,indx,d,code)

       ! --------------------------------------------------
       ! The metric becomes singular on axis (return code=1).
       ! Don't stop in these cases, but keep the NaN values.
       ! Fluxes through the axis are identically set to zero to catch this.
       ! (e.g. in subroutine tvdlf)
       ! These Nan's should thus not filter through to the solution.
       ! --------------------------------------------------
       if (code /= 1) then

          do j=1,^NC
             ! Directly calculate inverse metric:
             {b(^C)=kr(j,^C)|;}
             call LUBKSB(a,^NC,indx,b)
             do i=1,j
                m%gammainv(i,j)%elem(ix^D) = b(i)
             end do
          end do

       end if ! metric non-singular
       
       {^D& end do\}

    end if ! space is diagonal

  end subroutine LU
  !=======================================================================
  subroutine raise3_dysp(ixI^L, ixO^L, w, x, d, u)
    ! In dynamical spacetime framework
    ! Optimised square b^2=b^i b_i for contravariant input
    ! three-vectors u
     include 'amrvacdef.f'
     integer, intent(in)              :: ixI^L, ixO^L
     double precision, intent(in)     :: w(ixI^S, 1:nw)
     double precision, intent(in)     :: x(ixI^S, 1:ndim)
     double precision, intent(in)     :: d(ixI^S, 1:ndir)
     double precision, intent(out)    :: u(ixI^S, 1:ndir)
     double precision                 :: gamma_tmp(ixI^S, 1:3, 1:3)
     integer                          :: idir

     call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S, 1:3, 1:3))
     do idir = 1, ^NC
       gamma_tmp(ixO^S, idir, idir) = gamma_tmp(ixO^S, idir, idir)*w(ixO^S, psi_metric_)**4
     end do

     do idir = 1, ^NC
       u(ixO^S, idir) = d(ixO^S, idir) / gamma_tmp(ixO^S, idir, idir)
     end do
  end subroutine raise3_dysp
  !=============================================================================
  subroutine raise3(ixI^L,ixO^L,m,d,u)
    ! Takes a covariant three-vector with ndir components d (stands for down)
    ! and returns the contravariant one u with ndir components (u stands for up)
    include 'amrvacdef.f'

    integer, intent(in)                      :: ixI^L, ixO^L
    double precision, intent(in)             :: d(ixI^S,1:ndir)
    double precision, intent(out)            :: u(ixI^S,1:ndir)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: i, j
    !-----------------------------------------------------------------------------
  
    if (space_is_diagonal) then
       
       do i=1,^NC
          u(ixO^S,i) = m%gammainv(i,i)%elem(ixO^S) * d(ixO^S,i)
       end do
       
    else

       u(ixO^S,:) = 0.0d0
       do j=1,^NC
          do i=1,^NC
             u(ixO^S,i) = u(ixO^S,i) + m%gammainv(i,j)%elem(ixO^S)*d(ixO^S,j)
          end do
       end do

    end if ! space is diagonal

  end subroutine raise3

  subroutine raise3_ixD(ix^D,m,d,u)
    ! Takes a covariant three-vector with ndir components d (stands for down)
    ! and returns the contravariant one u with ndir components (u stands for up)
    include 'amrvacdef.f'

    integer, intent(in) :: ix^D
    double precision, intent(in)             :: d(1:ndir)
    double precision, intent(out)            :: u(1:ndir)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: i, j
    !-----------------------------------------------------------------------------

    if (space_is_diagonal) then

       do i=1,^NC
          u(i) = m%gammainv(i,i)%elem(ix^D) * d(i)
       end do

    else

       u(:) = 0.0d0
       do j=1,^NC
          do i=1,^NC
             u(i) = u(i) + m%gammainv(i,j)%elem(ix^D)*d(j)
          end do
       end do

    end if ! space is diagonal

  end subroutine raise3_ixD


  !=============================================================================
  subroutine raise3_point(x^D,d,u)
    ! Takes a covariant three-vector at coordinate x^D with ndir components d (stands for down)
    ! and returns the contravariant one u with ndir components (u stands for up)
    ! This is for pointwise-values!

    double precision, dimension(1:^NC), intent(in)    :: d
    double precision, dimension(1:^NC), intent(out)   :: u
    double precision, intent(in)                      :: x^D
    ! .. local ..
    double precision, dimension(1:^NC,1:^NC) :: gammainv
    integer                                  :: i, j
    !-----------------------------------------------------------------------------

    !stop 'never use raise3 pt'

    ! Component-wise. Could be more efficient.
    do j=1,^NC
       do i=1,j
          call get_gammainv_component(i,j,x^D,gammainv(i,j))
          if (i .ne. j) gammainv(j,i) = gammainv(i,j)
       end do
    end do

    ! raise:

    if (space_is_diagonal) then
       
       do i=1,^NC
          u(i) = gammainv(i,i) * d(i)
       end do

    else

       u = 0.0d0
       do i=1,^NC
          do j=1,^NC
             u(i) = u(i) + gammainv(i,j) * d(j)
          end do
       end do

    end if

  end subroutine raise3_point
  !=======================================================================
  subroutine lower3_dysp(ixI^L, ixO^L, w, x, u, d)
    ! In dynamical spacetime framework
    ! Optimised square b^2=b^i b_i for contravariant input
    ! three-vectors u
     include 'amrvacdef.f'
     integer, intent(in)              :: ixI^L, ixO^L
     double precision, intent(in)     :: w(ixI^S, 1:nw)
     double precision, intent(in)     :: x(ixI^S, 1:ndim)
     double precision, intent(in)     :: u(ixI^S, 1:ndir)
     double precision, intent(out)    :: d(ixI^S, 1:ndir)
     double precision                 :: gamma_tmp(ixI^S, 1:3, 1:3)
     integer                          :: idir

     call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S, 1:3, 1:3))
     do idir = 1, ^NC
       gamma_tmp(ixO^S, idir, idir) = gamma_tmp(ixO^S, idir, idir)*w(ixO^S, psi_metric_)**4
     end do

     do idir = 1, ^NC
       d(ixO^S, idir) = u(ixO^S, idir) * gamma_tmp(ixO^S, idir, idir)
     end do
  end subroutine lower3_dysp
  !=============================================================================
  subroutine lower3(ixI^L,ixO^L,m,u,d)
    ! Takes a contravariant three-vector with ndir components u (stands for up)
    ! and returns the covariant one d with ndir components (d stands for down)
    include 'amrvacdef.f'

    integer, intent(in)                      :: ixI^L, ixO^L
    double precision, intent(in)             :: u(ixI^S,1:ndir)
    double precision, intent(out)            :: d(ixI^S,1:ndir)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: inonzero, i, j
    !-----------------------------------------------------------------------------

    d(ixO^S,:) = 0.0d0
    
    do inonzero=1,m%nnonzero
       i = m%nonzero(inonzero)%i
       j = m%nonzero(inonzero)%j
       if (i==0 .or. j==0) cycle

       d(ixO^S,i) = d(ixO^S,i) + m%nonzero(inonzero)%elem(ixO^S) * u(ixO^S,j)
       
    end do
    
  end subroutine lower3

  subroutine lower3_ixD(ix^D,m,u,d)
    ! Takes a contravariant three-vector with ndir components u (stands for up)
    ! and returns the covariant one d with ndir components (d stands for down)
    include 'amrvacdef.f'

    integer, intent(in)                      :: ix^D
    double precision, intent(in)             :: u(1:ndir)
    double precision, intent(out)            :: d(1:ndir)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: inonzero, i, j
    !-----------------------------------------------------------------------------

    d(:) = 0.0d0
    do inonzero=1,m%nnonzero
       i = m%nonzero(inonzero)%i
       j = m%nonzero(inonzero)%j
       if (i==0 .or. j==0) cycle

       d(i) = d(i) + m%nonzero(inonzero)%elem(ix^D) * u(j)

    end do

  end subroutine lower3_ixD



    !=============================================================================
  subroutine lower3_point(x^D,u,d)
    ! Takes a contravariant three-vector at coordinate x^D with ndir components u (stands for up)
    ! and returns the covariant one d with ndir components (d stands for down)
    ! This is for pointwise-values!

    double precision, dimension(1:^NC), intent(in)    :: u
    double precision, dimension(1:^NC), intent(out)   :: d
    double precision, intent(in)                      :: x^D
    ! .. local ..
    double precision, dimension(1:^NC,1:^NC) :: g
    integer                                  :: i, j
    !-----------------------------------------------------------------------------

    !stop 'never use lower3 pt'

    ! Component-wise. Could be more efficient.
    do j=1,^NC
       do i=1,j
          call get_g_component(i,j,x^D,g(i,j))
          if (i .ne. j) g(j,i) = g(i,j)
       end do
    end do

    ! lower:

    if (space_is_diagonal) then
       
       do i=1,^NC
          d(i) = g(i,i) * u(i)
       end do

    else

       d = 0.0d0
       do i=1,^NC
          do j=1,^NC
             d(i) = d(i) + g(i,j) * u(j)
          end do
       end do

    end if

  end subroutine lower3_point
  !=======================================================================
  subroutine square3u_dysp(ixI^L, ixO^L, w, x, vec3d, vec_product)
    ! In dynamical spacetime framework
    ! Optimised square b^2=b^i b_i for contravariant input
    ! three-vectors u
     include 'amrvacdef.f'
     integer, intent(in)              :: ixI^L, ixO^L
     double precision, intent(in)     :: w(ixI^S, 1:nw)
     double precision, intent(in)     :: x(ixI^S, 1:ndim)
     double precision, intent(in)     :: vec3d(ixI^S, 1:ndir)
     double precision, intent(out)    :: vec_product(ixI^S)
     double precision                 :: covar_vec(ixI^S, 1:ndir)
     double precision                 :: gamma_tmp(ixI^S, 1:3, 1:3)
     integer                          :: idir

     call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S, 1:3, 1:3))
     do idir = 1, ndir
         gamma_tmp(ixO^S, idir, idir) = gamma_tmp(ixO^S, idir, idir)*w(ixO^S, psi_metric_)**4
     end do
     vec_product(ixO^S) = {^C& vec3d(ixO^S, ^C)**2*gamma_tmp(ixO^S, ^C, ^C)+}
  end subroutine square3u_dysp
  !=============================================================================
  subroutine square3u(ixI^L,ixO^L,m,u,u2)
    ! Optimised square b^2=b^i b_i for contravariant input
    ! three-vectors u

    include 'amrvacdef.f'

    integer, intent(in)                      :: ixI^L, ixO^L
    double precision, intent(in)             :: u(ixI^S,1:ndir)
    double precision, intent(out)            :: u2(ixI^S)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: i
    !-----------------------------------------------------------------------------

    u2(ixO^S) = zero
    
    ! Diagonal terms:
    do i = 1, ^NC
       u2(ixO^S) = u2(ixO^S) + u(ixO^S,i)**2 * m%g(i,i)%elem(ixO^S)
    end do
    
    ! Offdiagonal entries:
    {^NOONEC
    if (.not. g_is_zero(1,2)) & 
         u2(ixO^S) = u2(ixO^S) + 2.0d0 * m%g(1,2)%elem(ixO^S) * u(ixO^S,1)*u(ixO^S,2)
    }{^IFTHREEC
    if (.not. g_is_zero(1,3)) & 
         u2(ixO^S) = u2(ixO^S) + 2.0d0 * m%g(1,3)%elem(ixO^S) * u(ixO^S,1)*u(ixO^S,3)
    if (.not. g_is_zero(2,3)) & 
         u2(ixO^S) = u2(ixO^S) + 2.0d0 * m%g(2,3)%elem(ixO^S) * u(ixO^S,2)*u(ixO^S,3)
    }
    
  end subroutine square3u

  subroutine square3u_ixD(ix^D,m,u,u2)
    ! Optimised square b^2=b^i b_i for contravariant input
    ! three-vectors u
    ! for one pt of ix^D
    include 'amrvacdef.f'
    integer, intent(in) ::  ix^D
    double precision, intent(in)             :: u(1:ndir)
    double precision, intent(out)            :: u2
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: i
    !-----------------------------------------------------------------------------

    u2 = zero

    ! Diagonal terms:
    do i = 1, ^NC
       u2 = u2 + u(i)**2 * m%g(i,i)%elem(ix^D)
    end do

    ! Offdiagonal entries:
    {^NOONEC
    if (.not. g_is_zero(1,2)) &
         u2 = u2 + 2.0d0 * m%g(1,2)%elem(ix^D) * u(1)*u(2)
    }{^IFTHREEC
    if (.not. g_is_zero(1,3)) &
         u2 = u2 + 2.0d0 * m%g(1,3)%elem(ix^D) * u(1)*u(3)
    if (.not. g_is_zero(2,3)) &
         u2 = u2 + 2.0d0 * m%g(2,3)%elem(ix^D) * u(2)*u(3)
    }

  end subroutine square3u_ixD



  !=============================================================================
    subroutine acrossbU(ixI^L,ixO^L,m,idir,a,b,res)
      ! Do the crossproduct of two covariant vectors
      ! Result is contravariant
      ! res^idir = eta^{idir,j,k} 1/sqrtgamma a_j b_k
    include 'amrvacdef.f'

    integer, intent(in)                                        :: ixI^L,ixO^L
    double precision, dimension(ixI^S,1:^NC), intent(in)       :: a,b
    integer         , intent(in)                               :: idir
    type(metric), pointer                                      :: m
    double precision, intent(inout)                            :: res(ixI^S)
    ! .. local ..
    integer                                                    :: j,k
    !-----------------------------------------------------------------------------

    res(ixO^S) = zero
    do j=1,^NC
       do k=1,^NC
          if (j .eq. k .or. j .eq. idir .or. k .eq. idir) cycle
             res(ixO^S) = res(ixO^S) + &
                  1.0d0/m%sqrtgamma(ixO^S)*lvc(idir,j,k)*a(ixO^S,j)*b(ixO^S,k)
       end do
    end do

  end subroutine acrossbU
  !=============================================================================
  subroutine lower4(ixI^L,ixO^L,m,u,d)
    ! Takes a contravariant four-vector u (stands for up) 
    ! and returns the covariant one d (stands for down)
    ! Note: index of vector starts at zero!
    include 'amrvacdef.f'

    integer, intent(in)                      :: ixI^L, ixO^L
    double precision, intent(in)             :: u(ixI^S,0:ndir)
    double precision, intent(out)            :: d(ixI^S,0:ndir)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: inonzero, i, j
    !-----------------------------------------------------------------------------

    d(ixO^S,:) = 0.0d0

    do inonzero=1,nnonzero_metric
       i = m%nonzero(inonzero)%i
       j = m%nonzero(inonzero)%j

       d(ixO^S,i) = d(ixO^S,i) + m%nonzero(inonzero)%elem(ixO^S) * u(ixO^S,j)
       
    end do
    
  end subroutine lower4
  !=============================================================================

  subroutine raise4_dysp(ixI^L, ixO^L, w, x, d, u)
      ! Raise 4-dim 1-form d to 4-dim 1-vector u
      include 'amrvacdef.f'
      integer, intent(in)              :: ixI^L, ixO^L
      double precision, intent(in)     :: w(ixI^S, 1:nw)
      double precision, intent(in)     :: x(ixI^S, 1:ndim)
      double precision, intent(in)     :: d(ixI^S, 0:ndir)
      double precision, intent(out)    :: u(ixI^S, 0:ndir)
      ! local variables
      double precision                 :: gmunu_inv(ixI^S, 0:3, 0:3) 
      integer                          :: i, j 

      gmunu_inv = 0.0d0
      u = 0.0d0
      call get_gammainvij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gmunu_inv(ixI^S, 1:3, 1:3))
      gmunu_inv(ixO^S, 0, 0) = -1./w(ixO^S, alp_metric_)**2
      do i = 1, 3
          do j=1, 3
              gmunu_inv(ixO^S, i, j) = gmunu_inv(ixO^S, i, j)/w(ixO^S, psi_metric_)**4-&
                  w(ixO^S, beta_metric1_+i-1)*w(ixO^S, beta_metric1_+j-1)/w(ixO^S, alp_metric_)**2
          enddo
          gmunu_inv(ixO^S, 0, i) = w(ixO^S, beta_metric1_+i-1)/w(ixO^S, alp_metric_)**2
          gmunu_inv(ixO^S, i, 0) = w(ixO^S, beta_metric1_+i-1)/w(ixO^S, alp_metric_)**2
      end do

      do i = 0, ndir
          do j=0, ndir
              u(ixO^S, i) = u(ixO^S, i)+d(ixO^S, j)*gmunu_inv(ixO^S, i, j)
          enddo
      end do
  end subroutine raise4_dysp

  subroutine lower4_dysp(ixI^L, ixO^L, w, x, u, d)
      ! Lower 4-dim 1-vector u to 4-dim 1-form d
      include 'amrvacdef.f'
      integer, intent(in)              :: ixI^L, ixO^L
      double precision, intent(in)     :: w(ixI^S, 1:nw)
      double precision, intent(in)     :: x(ixI^S, 1:ndim)
      double precision, intent(in)     :: u(ixI^S, 0:ndir)
      double precision, intent(out)    :: d(ixI^S, 0:ndir)
      ! local variables
      double precision                 :: gmunu(ixI^S, 0:3, 0:3)
      double precision                 :: beta_3d_1form(ixI^S, 1:ndir), beta_contract(ixI^S)
      integer                          :: i, j

      gmunu = 0.0d0
      d = 0.0d0
      call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gmunu(ixI^S, 1:3, 1:3))
      call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, beta_metric1_:beta_metric3_), beta_3d_1form(ixI^S, 1:3))
      beta_contract(ixO^S) = beta_3d_1form(ixO^S, 1)*w(ixO^S, beta_metric1_)+beta_3d_1form(ixO^S, 2)* &
          w(ixO^S, beta_metric2_)+beta_3d_1form(ixO^S, 3)*w(ixO^S, beta_metric3_)
      gmunu(ixO^S, 0, 0) = -w(ixO^S, alp_metric_)**2+beta_contract(ixO^S)
      do i = 1, ndir
          do j=1, ndir
              gmunu(ixO^S, i, j) = gmunu(ixO^S, i, j)*w(ixO^S, psi_metric_)**4
          enddo
          gmunu(ixO^S, 0, i) = beta_3d_1form(ixO^S, i)
          gmunu(ixO^S, i, 0) = beta_3d_1form(ixO^S, i)
      end do

      do i = 0, ndir
          do j=0, ndir
              d(ixO^S, i) = d(ixO^S, i)+u(ixO^S, j)*gmunu(ixO^S, i, j)
          enddo
      end do
  end subroutine lower4_dysp

  subroutine get_K_ij_dysp(ixI^L, ixO^L, w, x, K_ij)
      ! Calculate 3-dim 2-form extrinsic curvature tensor K_ij
      use mod_cfc_parameters
      include 'amrvacdef.f'

      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: x(ixI^S, 1:ndim)
      double precision, intent(in)    :: w(ixI^S, 1:nw)
      double precision, intent(inout) :: K_ij(ixI^S,1:3,1:3)
      ! local variables
      {^NOONED
          double precision            :: cot_theta(ixI^S)
      }
      double precision                :: dbeta(ixI^S,1:3,1:3)
      double precision                :: gamma(ixI^S,1:3,1:3)
      integer                         :: idir, jdir

      dbeta = 0.0d0
      K_ij = 0.0d0
      ! calculate auxiliary

      call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
      do idir = 1, ndir
         gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
      end do
      do idir = 1, ndim
         do jdir = 1, ndir
            call partial_d(w(ixI^S, beta_metric1_+jdir-1), ixI^L, ixO^L, x, idir, dbeta(ixI^S,jdir,idir))
         end do
      end do
      ! calculate K_ij
      select case (coordinate)
      case (cartesian)
         K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( 2.0d0*dbeta(ixO^S,1,1) &
                        {^NOONED - dbeta(ixO^S,2,2) } {^IFTHREED - dbeta(ixO^S,3,3) } )

         K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( - dbeta(ixO^S,1,1) &
                        {^NOONED + 2.0d0*dbeta(ixO^S,2,2) } {^IFTHREED - dbeta(ixO^S,3,3) } )

         K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( - dbeta(ixO^S,1,1) &
                        {^NOONED - dbeta(ixO^S,2,2) } {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
      case (cylindrical)
         K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( 2.0d0*dbeta(ixO^S,1,1) - w(ixO^S, beta_metric1_)/x(ixO^S,r_) &
                        {^NOONED - dbeta(ixO^S,2,2) } {^IFTHREED - dbeta(ixO^S,3,3) } )

         K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( - dbeta(ixO^S,1,1) - w(ixO^S, beta_metric1_)/x(ixO^S,r_) &
                        {^NOONED + 2.0d0 * dbeta(ixO^S,2,2) } {^IFTHREED - dbeta(ixO^S,3,3) } )

         K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( - dbeta(ixO^S,1,1) + 2.0d0 * w(ixO^S, beta_metric1_)/x(ixO^S,r_) &
                        {^NOONED - dbeta(ixO^S,2,2) } {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
      case (spherical)
         {^NOONED
         cot_theta(ixO^S) = dcos(x(ixO^S,^Z))/dsin(x(ixO^S,^Z))
         !cot_theta(ixO^S) = dcos(x(ixO^S,theta_))/dsin(x(ixO^S,theta_))
         }

         K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( 2.0d0*dbeta(ixO^S,1,1) - 2.0d0*w(ixO^S, beta_metric1_)/x(ixO^S,1) &
                        {^NOONED - dbeta(ixO^S,2,2) - w(ixO^S, beta_metric2_)*cot_theta(ixO^S) } {^IFTHREED - dbeta(ixO^S,3,3) } )

         K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( -dbeta(ixO^S,1,1) + w(ixO^S, beta_metric1_)/x(ixO^S,1) &
                        {^NOONED + 2.0d0*dbeta(ixO^S,2,2) - w(ixO^S, beta_metric2_)*cot_theta(ixO^S) } {^IFTHREED - dbeta(ixO^S,3,3) } )

         K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*w(ixO^S, alp_metric_)) &
                      * ( -dbeta(ixO^S,1,1) + w(ixO^S, beta_metric1_)/x(ixO^S,1) &
                        {^NOONED - dbeta(ixO^S,2,2) + 2.0d0*w(ixO^S, beta_metric2_)*cot_theta(ixO^S) } {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
      end select
      {^NOONED
      K_ij(ixO^S,1,2) = ( dbeta(ixO^S,1,2)*gamma(ixO^S,1,1) + dbeta(ixO^S,2,1)*gamma(ixO^S,2,2) ) / (2.0d0*w(ixO^S, alp_metric_))
      K_ij(ixO^S,1,3) = ( dbeta(ixO^S,3,1)*gamma(ixO^S,3,3) + dbeta(ixO^S,1,3)*gamma(ixO^S,1,1) ) / (2.0d0*w(ixO^S, alp_metric_))
      K_ij(ixO^S,2,3) = ( dbeta(ixO^S,3,2)*gamma(ixO^S,3,3) + dbeta(ixO^S,2,3)*gamma(ixO^S,2,2) ) / (2.0d0*w(ixO^S, alp_metric_))
      ! K_ji=K_ij
      do idir=1,2
         do jdir=idir+1,3
            K_ij(ixO^S,jdir,idir) = K_ij(ixO^S,idir,jdir)
         end do
      end do
      }
  end subroutine get_K_ij_dysp

  subroutine get_Kij_dysp(ixI^L, ixO^L, w, x, Kij)
      ! Calculate 3-dim 2-vector extrinsic curvature K^ij from its 3-dim 2-form K_ij
      include 'amrvacdef.f'

      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: x(ixI^S, 1:ndim)
      double precision, intent(in)    :: w(ixI^S, 1:nw)
      double precision, intent(inout) :: Kij(ixI^S,1:3,1:3)
      ! local variables
      double precision                :: gammainv(ixI^S,1:3,1:3), K_ij(ixI^S,1:3,1:3)
      integer                         :: idir, jdir

      Kij = 0.0d0

      call get_gammainvij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gammainv(ixI^S,1:3,1:3))
      do idir = 1, ndir
         gammainv(ixO^S,idir,idir) = gammainv(ixO^S,idir,idir) / w(ixO^S, psi_metric_)**4
      end do
      call get_K_ij_dysp(ixI^L, ixO^L, w, x, K_ij)
      do idir=1, ndir
          do jdir=1, ndir
              Kij(ixO^S, idir, jdir) = K_ij(ixO^S, idir, jdir)*gammainv(ixO^S, idir, idir)*gammainv(ixO^S, jdir, jdir)
          end do
      end do
  end subroutine get_Kij_dysp

  subroutine allocate_metric(m,ixG^L,need_derivs)
    include 'amrvacdef.f'

    type(metric), pointer, intent(inout)     :: m
    integer, intent(in)                      :: ixG^L
    logical, optional                        :: need_derivs
    ! .. local ..
    integer                                  :: inonzero, i,j,k
    logical                                  :: alloc_derivs
    !-----------------------------------------------------------------------------
    if (.not.present(need_derivs)) then
       alloc_derivs = .true.
    else
       alloc_derivs = need_derivs
    end if

    
    allocate(m)
    allocate(m%nonzero(1:nnonzero_metric))
    allocate(m%nonzeroBeta(1:nnonzero_beta))
    allocate(m%g(0:^NC,0:^NC))
    allocate(m%beta(1:^NC))
    allocate(m%betaD(1:^NC))
    allocate(m%sqrtgamma(ixG^S))
    allocate(m%alpha(ixG^S))
    ! Derivatives:
    allocate(m%nonzeroDgDk(1:nnonzero_dgdk),m%DgDk(1:^NC,1:^NC,1:^NC))
    allocate(m%nonzeroDalphaDj(1:nnonzero_dalphadj),m%DalphaDj(1:^NC))
    allocate(m%nonzeroDbetaiDj(1:nnonzero_dbetaidj),m%DbetaiDj(1:^NC,1:^NC))
    ! Inverse metric:
    allocate(m%gammainv(1:^NC,1:^NC))
    
    ! Store shape of arrays:
    m%ixG^L=ixG^L;

    ! we also always allocate one element of zeros
    allocate(m%zero(ixG^S))
    m%zero = 0.0d0

    ! Inverse metric
       do i=1,^NC
          do j=1,^NC
             m%gammainv(i,j)%i = i
             m%gammainv(i,j)%j = j
             m%gammainv(i,j)%k = -1
             if (space_is_diagonal) then
                if (i .eq. j) then
                   allocate(m%gammainv(i,i)%elem(ixG^S))
                else
                   m%gammainv(i,j)%elem => m%zero
                end if
             else
                if (j .ge. i) then 
                   allocate(m%gammainv(i,j)%elem(ixG^S))
                else
                   m%gammainv(i,j)%elem => m%gammainv(j,i)%elem
                end if
             end if
          end do
       end do

       
    ! Allocate contravariant shifts m%beta
    inonzero = 0
    do j=1,^NC
       m%beta(j)%i = j
       m%beta(j)%j = j
       m%beta(j)%k = j
       if (beta_is_zero(j)) then 
          m%beta(j)%elem => m%zero
       else
          inonzero = inonzero+1
          allocate(m%beta(j)%elem(ixG^S))
          call associateIndexlist(m%nonzeroBeta(inonzero),m%beta(j))
       end if
    end do
    m%nnonzeroBeta = inonzero

    ! Covariant shifts point to metric:
    do j = 1,^NC
          nullify(m%betaD(j)%elem)
          m%betaD(j)%i = j
          m%betaD(j)%j = j
          m%betaD(j)%k = j
    end do

    ! Point them all to zero if there is at least one
    ! zero-element:
    if (any(g_is_zero(:,:))) then
       do i = 0,^NC
          do j = 0,^NC
             m%g(i,j)%i = i
             m%g(i,j)%j = j
             m%g(i,j)%k = -1
             m%g(i,j)%elem => m%zero
          end do
       end do
    else
       do i = 0,^NC
          do j = 0,^NC
             m%g(i,j)%i = i
             m%g(i,j)%j = j
             m%g(i,j)%k = -1
             nullify(m%g(i,j)%elem)
          end do
       end do
    end if

    !------------------------------------------------
    ! Initialize indices and allocate the metric derivatives:
    !------------------------------------------------
    inonzero = 0
    do k=1,^NC
       do i=1,^NC
          do j=1,^NC
             m%DgDk(i,j,k)%i = i
             m%DgDk(i,j,k)%j = j
             m%DgDk(i,j,k)%k = k
             nullify(m%DgDk(i,j,k)%elem)
             if (.not.dgdk_is_zero(i,j,k)) then
                inonzero = inonzero + 1
                m%nonzeroDgDk(inonzero)%i = i
                m%nonzeroDgDk(inonzero)%j = j
                m%nonzeroDgDk(inonzero)%k = k
                ! We need only to allocate the upper triangle with j>=i
                ! Link the other components to the transposed DgDk(j,i,k) structure
                if (j .ge. i) then
                   if (alloc_derivs) allocate(m%nonzeroDgDk(inonzero)%elem(ixG^S))
                   m%DgDk(i,j,k)%elem => m%nonzeroDgDk(inonzero)%elem
                else
                   m%nonzeroDgDk(inonzero)%elem => m%DgDk(j,i,k)%elem 
                   m%DgDk(i,j,k)%elem => m%DgDk(j,i,k)%elem
                end if
             else
                m%DgDk(i,j,k)%elem => m%zero
             end if
          end do
       end do
    end do
    m%nnonzeroDgDk = inonzero

    !------------------------------------------------
    ! Initialize indices and allocate the shift derivatives:
    !------------------------------------------------
    inonzero = 0
    do j=1,^NC
       do i=1,^NC
          m%DbetaiDj(i,j)%i = i
          m%DbetaiDj(i,j)%j = j
          m%DbetaiDj(i,j)%k = -1
          nullify(m%DbetaiDj(i,j)%elem)
          if (.not. dbetaidj_is_zero(i,j)) then
             inonzero = inonzero + 1
             m%nonzeroDbetaiDj(inonzero)%i = i
             m%nonzeroDbetaiDj(inonzero)%j = j
             m%nonzeroDbetaiDj(inonzero)%k = -1
             if (alloc_derivs) allocate(m%nonzeroDbetaiDj(inonzero)%elem(ixG^S))
             m%DbetaiDj(i,j)%elem => m%nonzeroDbetaiDj(inonzero)%elem
          else
             m%DbetaiDj(i,j)%elem => m%zero
          end if
       end do
    end do
    m%nnonzeroDbetaiDj = inonzero

    !------------------------------------------------
    ! Initialize indices and allocate the lapse derivatives:
    !------------------------------------------------
    inonzero = 0
    do j=1,^NC
       m%DalphaDj(j)%i = j
       m%DalphaDj(j)%j = j
       m%DalphaDj(j)%k = j
       nullify(m%DalphaDj(j)%elem)
       if (.not. dalphadj_is_zero(j)) then
             inonzero = inonzero + 1
             m%nonzeroDalphaDj(inonzero)%i = j
             m%nonzeroDalphaDj(inonzero)%j = j
             m%nonzeroDalphaDj(inonzero)%k = j
             if (alloc_derivs) allocate(m%nonzeroDalphaDj(inonzero)%elem(ixG^S))
             m%DalphaDj(j)%elem => m%nonzeroDalphaDj(inonzero)%elem
          else
             m%DalphaDj(j)%elem => m%zero
          end if
    end do
    m%nnonzeroDalphaDj = inonzero

    
    ! ==============================
    ! Allocate and link the non-zero metric components
    ! ==============================
    
    inonzero=0
    if (.not.g_is_zero(0,0)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixG^S))
       m%nonzero(inonzero)%i = 0
       m%nonzero(inonzero)%j = 0
       m%nonzero(inonzero)%k = -1
       m%g(0,0)%elem => m%nonzero(inonzero)%elem
    end if
    
    {^C&
    if (.not.g_is_zero(0,^C)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixG^S))
       m%nonzero(inonzero)%i = 0
       m%nonzero(inonzero)%j = ^C
       m%nonzero(inonzero)%k = -1
       m%g(0,^C)%elem => m%nonzero(inonzero)%elem
       ! Lower triangle:
       inonzero = inonzero + 1       
       m%nonzero(inonzero)%i = ^C
       m%nonzero(inonzero)%j = 0
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(0,^C)%elem
       m%g(^C,0)%elem => m%g(0,^C)%elem
    end if
    m%betaD(^C)%elem => m%g(0,^C)%elem
    |\}

    
    {^C&
    if (.not.g_is_zero(1,^C)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixG^S))
       m%nonzero(inonzero)%i = 1
       m%nonzero(inonzero)%j = ^C
       m%nonzero(inonzero)%k = -1
       m%g(1,^C)%elem => m%nonzero(inonzero)%elem
    end if
    |\}
       
    {#IFDEF C2
    if (.not.g_is_zero(2,2)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixG^S))
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 2
       m%nonzero(inonzero)%k = -1
       m%g(2,2)%elem => m%nonzero(inonzero)%elem
    end if
    if (.not.g_is_zero(2,1)) then
       inonzero = inonzero + 1
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 1
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(1,2)%elem
       m%g(2,1)%elem => m%g(1,2)%elem
    end if
    }
    
    {#IFDEF C3
    if (.not.g_is_zero(2,2)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixG^S))
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 2
       m%nonzero(inonzero)%k = -1
       m%g(2,2)%elem => m%nonzero(inonzero)%elem
    end if
    if (.not.g_is_zero(2,1)) then
       inonzero = inonzero + 1
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 1
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(1,2)%elem
       m%g(2,1)%elem => m%g(1,2)%elem
    end if
    if (.not.g_is_zero(2,3)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixG^S))
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 3
       m%nonzero(inonzero)%k = -1
       m%g(2,3)%elem => m%nonzero(inonzero)%elem
       ! Lower triangle
       inonzero = inonzero + 1
       m%nonzero(inonzero)%i = 3
       m%nonzero(inonzero)%j = 2
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(2,3)%elem
       m%g(3,2)%elem => m%g(2,3)%elem
    end if
    
    if (.not.g_is_zero(3,3)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixG^S))
       m%nonzero(inonzero)%i = 3
       m%nonzero(inonzero)%j = 3
       m%nonzero(inonzero)%k = -1
       m%g(3,3)%elem => m%nonzero(inonzero)%elem
    end if

    if (.not.g_is_zero(3,1)) then
       inonzero = inonzero + 1
       m%nonzero(inonzero)%i = 3
       m%nonzero(inonzero)%j = 1
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(1,3)%elem
       m%g(3,1)%elem => m%g(1,3)%elem
    end if
    }

    m%nnonzero  = inonzero

    if (inonzero .ne. nnonzero_metric) &
         call mpistop('allocate_metric: wrong number of non-zero elements')
    
  end subroutine allocate_metric
  !=============================================================================
  subroutine associateIndexlist(src,tgt)
    include 'amrvacdef.f'
    
    type(indexlist), intent(inout)                    :: src,tgt
    !-----------------------------------------------------------------------------

    src%elem => tgt%elem
    src%i = tgt%i
    src%j = tgt%j
    src%k = tgt%k

  end subroutine associateIndexlist
  !=============================================================================
  subroutine deallocate_metric(m)
    include 'amrvacdef.f'

    type(metric), pointer, intent(inout)   :: m
    integer                                :: i,j, inonzero
    !----------------------------------------
    ! Free metric (and covariant shift) memory:
    !----------------------------------------
    do inonzero=1,m%nnonzero
       ! Only deallocate upper triangle, the others are just links
       if (m%nonzero(inonzero)%j .ge. m%nonzero(inonzero)%i) then
          deallocate(m%nonzero(inonzero)%elem)
       end if
       nullify(m%nonzero(inonzero)%elem)
    end do
    deallocate(m%nonzero); nullify(m%nonzero)

    ! Nullify covariant shifts (elements of four-metric), lived in m%nonzero
    deallocate(m%betaD); nullify(m%betaD)
    
    ! Nullify four-metric, lived in m%nonzero
    deallocate(m%g); nullify(m%g)
    !----------------------------------------

    
    !----------------------------------------
    ! Metric derivatives:
    !----------------------------------------
    do inonzero=1,m%nnonzeroDgDk
       ! Deallocate only upper triangle in i,j, the others are just links
       if (m%nonzeroDgDk(inonzero)%j .ge. m%nonzeroDgDk(inonzero)%i &
            .and. associated(m%nonzeroDgDk(inonzero)%elem)) then
          deallocate(m%nonzeroDgDk(inonzero)%elem)
       end if
       nullify(m%nonzeroDgDk(inonzero)%elem)
    end do
    deallocate(m%nonzeroDgDk); nullify(m%nonzeroDgDk)
    deallocate(m%DgDk); nullify(m%DgDk)
    !----------------------------------------

    
    !----------------------------------------
    ! Lapse derivatives:
    !----------------------------------------
    do inonzero=1,m%nnonzeroDalphaDj
       if (associated(m%nonzeroDalphaDj(inonzero)%elem)) deallocate(m%nonzeroDalphaDj(inonzero)%elem)
       nullify(m%nonzeroDalphaDj(inonzero)%elem)
    end do
    deallocate(m%nonzeroDalphaDj); nullify(m%nonzeroDalphaDj)
    deallocate(m%DalphaDj); nullify(m%DalphaDj)
    !----------------------------------------

    
    !----------------------------------------
    ! Shift derivatives:
    !----------------------------------------
    do inonzero=1,m%nnonzeroDbetaiDj
       if (associated(m%nonzeroDbetaiDj(inonzero)%elem)) deallocate(m%nonzeroDbetaiDj(inonzero)%elem)
       nullify(m%nonzeroDbetaiDj(inonzero)%elem)
    end do
    deallocate(m%nonzeroDbetaiDj); nullify(m%nonzeroDbetaiDj)
    deallocate(m%DbetaiDj); nullify(m%DbetaiDj)
    !----------------------------------------
    
    
    !----------------------------------------
    ! deallocate or just nullify contravariant shift:
    !----------------------------------------
    do i=1,^NC
       if (.not.beta_is_zero(i)) then 
          deallocate(m%beta(i)%elem); nullify(m%beta(i)%elem)
       else
          nullify(m%beta(i)%elem)
       end if
    end do
    deallocate(m%nonzeroBeta); nullify(m%nonzeroBeta)
    deallocate(m%beta); nullify(m%beta)
    !----------------------------------------

    
    !----------------------------------------
    ! Lapse and shift:
    !----------------------------------------
    deallocate(m%alpha,m%sqrtgamma); nullify(m%alpha,m%sqrtgamma)
    !----------------------------------------

    
    !----------------------------------------
    ! inverse Metric:
    !----------------------------------------
    if (space_is_diagonal) then
       do i=1,^NC
          deallocate(m%gammainv(i,i)%elem); nullify(m%gammainv(i,i)%elem)
       end do
    else
       do i=1,^NC
          do j=1,^NC
             ! Deallocate only upper triangle in i,j, the others are just links
             if (j .ge. i) then
                deallocate(m%gammainv(i,j)%elem)
             end if
             nullify(m%gammainv(i,j)%elem)
          end do
       end do
    end if
    deallocate(m%gammainv); nullify(m%gammainv)
    !----------------------------------------


    deallocate(m%zero); nullify(m%zero)

    deallocate(m); nullify(m)

  end subroutine deallocate_metric
  !=============================================================================
  subroutine fill_metric(m,ixGext^L,ixO^L,x,w,need_derivs)
    include 'amrvacdef.f'

    type(metric), pointer                                   :: m
    integer, intent(in)                                     :: ixGext^L, ixO^L
    double precision, dimension(ixGext^S,1:ndim), intent(in):: x
    logical, optional                                       :: need_derivs
    double precision, dimension(ixGext^S,1:nw) , intent(in) :: w

    ! .. local ..
    integer                                                 :: ix^D, i, j, k, inonzero
    double precision, dimension(ixGext^S,1:ndir)            :: betaD ! covariant shift
    double precision, dimension(ixGext^S,1:ndir)            :: betaU ! contravariant sift
    double precision, dimension(ixGext^S)                   :: beta2
    double precision                                        :: dummy
    double precision, dimension(ixGext^S,1:^NC,1:^NC,1:^NC) :: dgammainvdk
    double precision, dimension(ixGext^S,0:^NC,1:^NC)       :: dbetaDownjDk
    logical                                                 :: fill_derivs
    !-----------------------------------------------------------------------------
    if (.not.present(need_derivs)) then
       fill_derivs = .true.
    else
       fill_derivs = need_derivs
    end if
   
    {#IFDEF DY_SP
        fill_derivs = .false.
    }
 
    !--------------------------------------------------
    if (.not. init_from_g4) then
    !--------------------------------------------------


       ! set the lapse:
       {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}

       {#IFDEF DY_SP
       call get_alpha({x(ix^DD,^D)},m%alpha(ix^D),w_pt=w(ix^D,1:nw))
       }
       {#IFNDEF DY_SP
       call get_alpha({x(ix^DD,^D)},m%alpha(ix^D))
       }

       {^D& end do \}

       ! set the shift (contravariant):
       {^C&
            if (.not. beta_is_zero(^C)) then
       {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}

       {#IFDEF DY_SP
       call get_beta(^C,{x(ix^DD,^D)},m%beta(^C)%elem(ix^D),w_pt=w(ix^D,1:nw))
       }
       {#IFNDEF DY_SP
       call get_beta(^C,{x(ix^DD,^D)},m%beta(^C)%elem(ix^D))
       }

       {^D& end do \}
       end if
       \}
       
       ! set the spatial part:
       do inonzero=1,nnonzero_metric
          i = m%nonzero(inonzero)%i
          j = m%nonzero(inonzero)%j

          ! Just spatial and only upper triangle (the others are already linking here)
          if (i == 0 .or. j == 0 .or. i .gt. j) cycle

          {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}

          {#IFDEF DY_SP
          call get_g_component(i,j,{x(ix^DD,^D)},m%nonzero(inonzero)%elem(ix^D),w_pt=w(ix^D,1:nw))
          }
          {#IFNDEF DY_SP
          call get_g_component(i,j,{x(ix^DD,^D)},m%nonzero(inonzero)%elem(ix^D))
          }

          {^D& end do \}
       end do

       ! set the space-time part (covariant shifts):
       {^C& betaU(ixO^S,^C) = m%beta(^C)%elem(ixO^S)\}
       call lower3(ixGext^L,ixO^L,m,betaU,betaD)
       do i=1,^NC
          if (.not. g_is_zero(0,i)) &
               m%g(0,i)%elem(ixO^S) = betaD(ixO^S,i)
       end do
       
       ! set the 00-component from the shift and lapse
       beta2(ixO^S) = {^C& betaD(ixO^S,^C)*betaU(ixO^S,^C) |+}
       m%g(0,0)%elem(ixO^S) = - m%alpha(ixO^S)**2 + beta2(ixO^S)

       ! set the square root of the spatial determinant
       {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}

       {#IFDEF DY_SP
       call get_sqrtgamma({x(ix^DD,^D)},m%sqrtgamma(ix^D),w_pt=w(ix^D,1:nw))
       }
       {#IFNDEF DY_SP
       call get_sqrtgamma({x(ix^DD,^D)},m%sqrtgamma(ix^D))
       }

       {^D& end do \}

       ! set the inverse of the three-metric
       call LU(m)

       !------------------------------------------------------------
       ! Fill the derivative data:
       !------------------------------------------------------------
       if (fill_derivs) then
          {#IFDEF DY_SP
            call mpistop('cfc does not need fill_derivs')
          }

          ! set the lapse:
          do inonzero=1, m%nnonzeroDalphaDj
             j = m%nonzeroDalphaDj(inonzero)%j
             {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
             call get_alpha({x(ix^DD,^D)},dummy,dalphadj=m%nonzeroDalphaDj(inonzero)%elem(ix^D),jdir=j)
             {^D& end do \}
          end do

          ! set the shift (contravariant):
          do inonzero=1,m%nnonzeroDbetaiDj
             i = m%nonzeroDbetaiDj(inonzero)%i
             j = m%nonzeroDbetaiDj(inonzero)%j
             {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
             call get_beta(i,{x(ix^DD,^D)},dummy,dbetaidj=m%nonzeroDbetaiDj(inonzero)%elem(ix^D),jdir=j)
             {^D& end do \}
          end do

          ! set the spatial part:
          do inonzero=1,m%nnonzeroDgDk
             i = m%nonzeroDgDk(inonzero)%i
             j = m%nonzeroDgDk(inonzero)%j
             k = m%nonzeroDgDk(inonzero)%k

             ! Just spatial and only upper triangle (the others are already linking here)
             if (i == 0 .or. j == 0 .or. i .gt. j) cycle

             {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
             call get_g_component(i,j,{x(ix^DD,^D)},dummy,dgdk=m%nonzeroDgDk(inonzero)%elem(ix^D),kdir=k)
             {^D& end do \}
          end do

       end if

    !--------------------------------------------------
    else ! init_from_g4
    !--------------------------------------------------
          {#IFDEF DY_SP
            call mpistop('cfc does not need init_from_g4')
          }

       ! set the full four metric:
       do inonzero=1,nnonzero_metric
          i = m%nonzero(inonzero)%i
          j = m%nonzero(inonzero)%j

          ! Only upper triangle (the others are already linking here)
          if (i .gt. j) cycle

          {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
          call get_g_component(i,j,{x(ix^DD,^D)},m%nonzero(inonzero)%elem(ix^D))

          {^D& end do \}
       end do
       
       ! set the inverse of the three-metric
       call LU(m)

       ! raise the shift vector and fill beta:
       do i=1,^NC
          betaD(ixO^S,i) = m%g(0,i)%elem(ixO^S)
       end do
       call raise3(ixGext^L,ixO^L,m,betaD,betaU)

       {^C& m%beta(^C)%elem(ixO^S) = betaU(ixO^S,^C) \}

       ! Obtain the lapse from shift and g00:
       beta2(ixO^S) = {^C& betaD(ixO^S,^C)*betaU(ixO^S,^C) |+}

       ! Freeze regions where the lapse becomes negative
       ! (such region should be protected by a horizon)
       m%alpha(ixO^S) = sqrt(max(beta2(ixO^S)-m%g(0,0)%elem(ixO^S),0.0d0))

       ! set the square root of the spatial determinant
       {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
       call get_sqrtgamma({x(ix^DD,^D)},m%sqrtgamma(ix^D))
       {^D& end do \}
       
       !------------------------------------------------------------
       ! Fill the derivative data:
       !------------------------------------------------------------

       if (fill_derivs) then
          {#IFDEF DY_SP
            call mpistop('cfc does not need fill_derivs')
          }

          ! set the spatial part:
          do inonzero=1,m%nnonzeroDgDk
             i = m%nonzeroDgDk(inonzero)%i
             j = m%nonzeroDgDk(inonzero)%j
             k = m%nonzeroDgDk(inonzero)%k

             ! Just spatial and only upper triangle (the others are already linking here)
             if (i == 0 .or. j == 0 .or. i .gt. j) cycle

             {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
             call get_g_component(i,j,{x(ix^DD,^D)},dummy,dgdk=m%nonzeroDgDk(inonzero)%elem(ix^D),kdir=k)
             {^D& end do \}
          end do


          !------------------------------
          ! Calculate derivatives of lapse and contravariant shift:
          !------------------------------

          ! Numerically compute the derivatives of the inverse three-metric (later not needed by scheme):
          do k=1,^NC
             {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
             call get_dgammainvdk({x(ix^DD,^D)},dgammainvdk(ix^D,1:^NC,1:^NC,k),k)
             {^D& end do \}
          end do

          ! Store the derivatives of the covariant shift and gtt (later not needed by scheme):
          do j=0,^NC
             do k=1,^NC
                {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
                call get_g_component(0,j,{x(ix^DD,^D)},dummy,dgdk=dbetaDownjDk(ix^D,j,k),kdir=k)
                {^D& end do \}
             end do
          end do

          ! Get derivative of contravariant shift:
          ! d_j beta^i = d_j gamma^{ik}*beta_k
          do inonzero=1,m%nnonzeroDbetaiDj
             i = m%nonzeroDbetaiDj(inonzero)%i
             j = m%nonzeroDbetaiDj(inonzero)%j
             m%nonzeroDbetaiDj(inonzero)%elem(ixO^S) = &
                  + {^C& betaD(ixO^S,^C)*dgammainvdk(ixO^S,i,^C,j)|+} &
                  + {^C& m%gammainv(i,^C)%elem(ixO^S)*dbetaDownjDk(ixO^S,^C,j)|+}
          end do

          ! Derivative of lapse:
          ! d_j alpha=1/2*(beta^k*d_j beta_k + beta_k*d_j beta^k - d_j g_{00})
          do inonzero=1, m%nnonzeroDalphaDj
             j = m%nonzeroDalphaDj(inonzero)%j

             m%nonzeroDalphaDj(inonzero)%elem(ixO^S) = &
                  0.5d0*({^C& betaU(ixO^S,^C)*dbetaDownjDk(ixO^S,^C,j) |+} &
                  + {^C& betaD(ixO^S,^C)*m%DbetaiDj(^C,j)%elem(ixO^S)|+} &
                  - dbetaDownjDk(ixO^S,0,j) &
                  )
          end do

       end if

    !--------------------------------------------------
    end if
    !--------------------------------------------------

  end subroutine fill_metric
  !=============================================================================

! xloc is different for sqrtgamma, it should be different for w(ix^D,1:nw) for different xloc
  subroutine int_surface(ixI^L,ixO^L,idims,x,dx^D,func,integral,w)
    ! integrates the function func over the cell-surface using 
    ! Simpsons rule. The cell extends over (x-dx/2, x+dx/2)
    include 'amrvacdef.f'
    integer, intent(in)                                               :: ixI^L, ixO^L, idims
    double precision, dimension(ixI^S,1:^ND), intent(in)              :: x
    double precision, intent(in)                                      :: dx^D
    double precision, dimension(ixI^S), intent(out)                   :: integral
    double precision                                                  :: func
    double precision, dimension(ixI^S,1:nw) , intent(in)              :: w

    ! .. local ..
    integer                :: is^D, ix^D
    logical, save          :: initialize=.true.
    double precision       :: tmp, xloc^D
{^NOONED     double precision, save :: s({^DE&-1:1}) }
{#IFDEF D3 double precision :: stmp(-1:1)}
    !-----------------------------------------------------------------------------
{^NOONED
    if (initialize) then 
       ! fill matrices for Simpsons rule:
       {#IFDEF D2 s = (/1.0d0,4.0d0,1.0d0/)}
       {#IFDEF D3 stmp = (/1.0d0,4.0d0,1.0d0/)
       do is2=-1,1
          do is1=-1,1
             s(is1,is2)=stmp(is1)*stmp(is2)
          end do
       end do
       }
       initialize = .false.
    end if
}

{#IFDEF D1
do ix1 = ixOmin1,ixOmax1
   {#IFDEF DY_SP
   call get_sqrtgamma(x(ix1,1),tmp,w_pt=w(ix^D,1:nw))
   }
   {#IFNDEF DY_SP
   call get_sqrtgamma(x(ix1,1),tmp)
   }
   integral(ix^D) = tmp * func(x(ix1,1))
end do
}

{#IFDEF D2
    select case(idims)
       case(1)
       {^D& do ix^DB = ixO^LIM^DB\}
       integral(ix^D) = 0.0d0
       xloc1 = x(ix^D,1)
       do is2 = -1,1
          xloc2 = x(ix^D,2)+dx2/2.0d0*dble(is2)
          {#IFDEF DY_SP
          call get_sqrtgamma(xloc^D,tmp,w_pt=w(ix^D,1:nw))
          }
          {#IFNDEF DY_SP
          call get_sqrtgamma(xloc^D,tmp)
          }
          integral(ix^D) = integral(ix^D) + tmp*func(xloc^D)*s(is2)
       end do
       integral(ix^D) = integral(ix^D) * dx2/6.0d0
       {^D& end do\}

       case(2)
       {^D& do ix^DB = ixO^LIM^DB\}
       integral(ix^D) = 0.0d0
       xloc2 = x(ix^D,2)
       do is1 = -1,1
          xloc1 = x(ix^D,1)+dx1/2.0d0*dble(is1)
          {#IFDEF DY_SP
          call get_sqrtgamma(xloc^D,tmp,w_pt=w(ix^D,1:nw))
          }
          {#IFNDEF DY_SP
          call get_sqrtgamma(xloc^D,tmp)
          }
          integral(ix^D) = integral(ix^D) + tmp*func(xloc^D)*s(is1)
       end do
       integral(ix^D) = integral(ix^D) * dx1/6.0d0
       {^D& end do\}
    end select
}
    
{#IFDEF D3
    select case(idims)
       case(1)
       {^D& do ix^DB = ixO^LIM^DB\}
       integral(ix^D) = 0.0d0
       xloc1 = x(ix^D,1)
       do is2 = -1,1
          xloc2 = x(ix^D,2)+dx2/2.0d0*dble(is2)
          do is3 = -1,1
             xloc3 = x(ix^D,3)+dx3/2.0d0*dble(is3)
             {#IFDEF DY_SP
             call get_sqrtgamma(xloc^D,tmp,w_pt=w(ix^D,1:nw))
             }
             {#IFNDEF DY_SP
             call get_sqrtgamma(xloc^D,tmp)
             }
             integral(ix^D) = integral(ix^D) + tmp*func(xloc^D)*s(is2,is3)
          end do
       end do
       integral(ix^D) = integral(ix^D) * dx2*dx3/36.0d0
       {^D& end do\}

       case(2)
       {^D& do ix^DB = ixO^LIM^DB\}
       integral(ix^D) = 0.0d0
       xloc2 = x(ix^D,2)
       do is1 = -1,1
          xloc1 = x(ix^D,1)+dx1/2.0d0*dble(is1)
          do is3 = -1,1
             xloc3 = x(ix^D,3)+dx3/2.0d0*dble(is3)
             {#IFDEF DY_SP
             call get_sqrtgamma(xloc^D,tmp,w_pt=w(ix^D,1:nw))
             }
             {#IFNDEF DY_SP
             call get_sqrtgamma(xloc^D,tmp)
             }
             integral(ix^D) = integral(ix^D) + tmp*func(xloc^D)*s(is1,is3)
          end do
       end do
       integral(ix^D) = integral(ix^D) * dx1*dx3/36.0d0
       {^D& end do\}

       case(3)
       {^D& do ix^DB = ixO^LIM^DB\}
       integral(ix^D) = 0.0d0
       xloc3 = x(ix^D,3)
       do is1 = -1,1
          xloc1 = x(ix^D,1)+dx1/2.0d0*dble(is1)
          do is2 = -1,1
             xloc2 = x(ix^D,2)+dx2/2.0d0*dble(is2)
             {#IFDEF DY_SP
             call get_sqrtgamma(xloc^D,tmp,w_pt=w(ix^D,1:nw))
             }
             {#IFNDEF DY_SP
             call get_sqrtgamma(xloc^D,tmp)
             }
             integral(ix^D) = integral(ix^D) + tmp*func(xloc^D)*s(is1,is2)
          end do
       end do
       integral(ix^D) = integral(ix^D) * dx1*dx2/36.0d0
       {^D& end do\}
    end select
}

  end subroutine int_surface
  !=============================================================================
  subroutine int_volume(ixI^L,ixO^L,x,dx^D,func,integral,w)
    include 'amrvacdef.f'
    ! integrates the function func over the cell-volume using
    ! Simpsons rule.  The cell extends over (x-dx/2, x+dx/2)
    integer, intent(in)                                               :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)              :: x
    double precision, intent(in)                                      :: dx^D
    double precision, dimension(ixI^S), intent(out)                   :: integral
    double precision                                                  :: func
    double precision, dimension(ixI^S,1:nw) , intent(in)              :: w

    ! .. local ..
    integer                    :: is^D, ix^D
    double precision           :: tmp
    double precision,save      :: s1D(-1:1) {#IFDEF D1 , s(-1:1)} {#IFDEF D2 , s(-1:1,-1:1)}
    {#IFDEF D3 double precision,save      ::  s2D(-1:1,-1:1), s(-1:1,-1:1,-1:1)}
    logical, save                         :: initialize=.true.
    !-----------------------------------------------------------------------------
    if (initialize) then
       ! fill matrices for Simpsons rule:
       s1D = (/1.0d0,4.0d0,1.0d0/)
       {#IFDEF D1 s(:) = s1D(:)}
       {^NOONED
       do is2=-1,1
          do is1=-1,1
             {#IFDEF D2
             s(is1,is2)=s1D(is1)*s1D(is2)
             }{#IFDEF D3
             s2D(is1,is2)=s1D(is1)*s1D(is2)
             }
          end do
       end do
       }{^IFTHREED
       {do is^DB=-1,1\}
                s(is^D)=s2D(is1,is2)*s1D(is3)
       {end do\}
       }
       initialize = .false.
    end if

    !--------------------------------
    ! Now integrate the function over the volume
    ! Thus obtain an approximation for
    ! I=\int\int\int sqrt(gamma)*func dx dy dz
    !--------------------------------
    
    {^D& do ix^DB = ixO^LIM^DB\}
    integral(ix^D) = 0.0d0

    {^D& do is^DB=-1,1\}
    {#IFDEF DY_SP
    call get_sqrtgamma({x(ix^DD,^D)+dx^D/2.0d0*dble(is^D)},tmp,w_pt=w(ix^D,1:nw))
    }
    {#IFNDEF DY_SP
    call get_sqrtgamma({x(ix^DD,^D)+dx^D/2.0d0*dble(is^D)},tmp)
    }
    integral(ix^D) = integral(ix^D) &
         + tmp*func({x(ix^DD,^D)+dx^D/2.0d0*dble(is^D)})*s(is^D)
    {^D& end do \}

    integral(ix^D) = integral(ix^D) * {^D& dx^D/6.0d0|*}

    {^D& end do\}
    
  end subroutine int_volume
  !=============================================================================
  subroutine complex_derivative(x^D,func,jdir,derivative)
    ! Performs complex-step derivative for accurate numerical derivatives 
    ! of an analytic real-valued function func.  
    ! This simply takes advantage of the Cauchy-Riemann relations.  
    ! 
    ! Interface for func must look like this:
    !
    ! double complex function func(x^D)
    ! double complex, intent(in)        :: x^D

    include 'amrvacdef.f'
    
    double precision, intent(in)        :: x^D
    integer, intent(in)                 :: jdir
    double precision, intent(out)       :: derivative

    interface
       double complex function func(x^D)
         double complex, intent(in)     :: x^D
       end function func
    end interface
    
    !  ..local..
    double complex                      :: xloc^D
    !-----------------------------------------------------------------------------

    xloc^D=cmplx(x^D , kr(jdir,^D)*smalldouble , kind(1.d0));
    
    derivative = aimag(func(xloc^D))/smalldouble
    
  end subroutine complex_derivative
  !=============================================================================
  double precision function ones(x^D)
    double precision        :: x^D
    !-----------------------------------------------------------------------------
    
    ones = 1.0d0

  end function ones
  !=============================================================================
  {^D&
  double precision function coordinate_x^D(x^DD)
    double precision        :: x^DD
    !-----------------------------------------------------------------------------
    
    coordinate_x^D = x^D

  end function coordinate_x^D
  !=============================================================================}
  subroutine fillgeo_covariant(pgeogrid,pwgrid,ixG^L,ixGext^L,xmin^D,dx^D,need_only_volume)
    use mod_interpolate
    include 'amrvacdef.f'

    type(geoalloc) :: pgeogrid
    type(walloc) :: pwgrid
    integer, intent(in) :: ixG^L, ixGext^L
    double precision, intent(in) :: xmin^D, dx^D
    logical, intent(in) :: need_only_volume
    ! .. local ..
    integer :: idims, idir, ixM^L, ix^L, ixC^L, ixF^L, ix^D, ix,&
               ixGsurf^L, ixMsurf^L
    double precision :: x(ixGextmin^D-1:ixGextmax^D,1:ndim),xi(ixGextmin^D-1:ixGextmax^D,1:ndim)
    double precision :: wi(ixGextmin^D-1:ixGextmax^D,1:nw)  !cell interface metric
    double precision :: wC(ixGextmin^D-1:ixGextmax^D,1:nw)

    !-----------------------------------------------------------------------------
    ixM^L=ixG^L^LSUBdixB;
    ix^L=ixG^L^LSUB1;
    ixGsurfmin^D=ixGextmin^D-1;
    ixGsurfmax^D=ixGextmax^D;
    ixMsurf^L=ixGsurf^L^LSUBdixB;

    {^D& 
    if(ixGsurfmin^D .eq. 0) then
      ixGsurfmin^D  = 1
    end if 
    \}
    
    !--------------------------------------------------
    ! Cell center positions:
    !--------------------------------------------------
    do idims=1,ndim
       select case(idims)
          {case(^D)
          do ix = ixGsurf^LIM^D
             x(ix^D%ixGsurf^S,^D)=xmin^D+(dble(ix-dixB)-half)*dx^D
          end do\}
       end select
    end do
   
    {#IFDEF DY_SP
    ! As surf is Gextmin - 1
    ! FIXME: which one we actually need 
    ! ixGsurf is correct
     wC(ixGext^S,1:nw) = pwgrid%w(ixGext^S,1:nw)
    !wC(ixGsurf^S,1:nw) = pwgrid%w(ixGsurf^S,1:nw)
    }

    !--------------------------------------------------
    ! Calculate the cell-volume with Simpsons rule:
    !--------------------------------------------------
    call int_volume(ixGext^L,ixGext^L,x(ixGext^S,1:ndim),dx^D,ones,pgeogrid%dvolume,wC(ixGext^S,1:nw))     

    !--------------------------------------------------
    ! Calculate the barycenter, also with Simpsons-rule:
    !--------------------------------------------------
    {^D&
    call int_volume(ixGext^L,ixGext^L,x(ixGext^S,1:ndim),dx^DD,coordinate_x^D,pgeogrid%xbar(ixGext^S,^D),wC(ixGext^S,1:nw))
    pgeogrid%xbar(ixGext^S,^D) = pgeogrid%xbar(ixGext^S,^D) / pgeogrid%dvolume(ixGext^S)
    \}

    !--------------------------------------------------
    ! Fill metric at the barycenter position:
    !--------------------------------------------------
    call fill_metric(pgeogrid%m,ixGext^L,ixGext^L,pgeogrid%xbar,wC(ixGext^S,1:nw))

    !--------------------------------------------------
    ! Get out when only volume is required
    !--------------------------------------------------
    if (need_only_volume) return


    !--------------------------------------------------
    ! Fill metric and surfaces at the interface/barycenter position:
    !--------------------------------------------------
    
    do idims=1,ndim
       select case(idims)
          {case(^D)
             !  not use CFC, you can turn interpolation off
             if (metric_vars_interpolation) then
              {^DD& xi(ixGsurf^S,^DD) = x(ixGsurf^S,^DD) + 0.5d0*dx^DD*kr(^DD,idims)\}
            
              ! Assume the first 1-2 ghostcell surface/ vol will not be used 
              ! call metric_interpolation(ixGsurf^L,ixGsurf^L,idims,wC(ixGsurf^S,1:nw),&  !not working output ixGsurf
              !           x(ixGsurf^S,1:ndim),wi(ixGsurf^S,1:nw),xi(ixGsurf^S,1:ndim))
              !  Fill number 3, 4 ghostcell of wi for surface
               call metric_interpolation(ixGsurf^L,ixMsurf^L^LADD2,idims,wC(ixGsurf^S,1:nw),&
                         x(ixGsurf^S,1:ndim),wi(ixGsurf^S,1:nw),xi(ixGsurf^S,1:ndim))

              ixCmin^DD=ixmin^DD-kr(^D,^DD); ixCmax^DD=ixmax^DD;
              ixF^L=ixC^L^LADD1;
              ! Metric at interface xi:
              call fill_metric(pgeogrid%mSurface^D,ixGsurf^L,ixGsurf^L,xi(ixGsurf^S,1:ndim),wi(ixGsurf^S,1:nw),need_derivs=.false.)
              ! Surface at interface xi:
              call int_surface(ixF^L,ixF^L,idims,xi(ixF^S,1:ndim),dx^DD,ones,pgeogrid%surfaceC^D,wi(ixF^S,1:nw))

              ! Use wC, as pgeogrid%surface is cell center surface, not interface
              ! Surface at barycenter:
              call int_surface(ixC^L,ixC^L,idims,pgeogrid%xbar(ixC^S,1:ndim),dx^DD,ones,pgeogrid%surface^D,wC(ixC^S,1:nw))
             else
              {^DD& xi(ixGsurf^S,^DD) = x(ixGsurf^S,^DD) + 0.5d0*dx^DD*kr(^DD,idims)\}
              
              ixCmin^DD=ixmin^DD-kr(^D,^DD); ixCmax^DD=ixmax^DD;
              ixF^L=ixC^L^LADD1;
              ! Metric at interface xi:
              call fill_metric(pgeogrid%mSurface^D,ixGsurf^L,ixGsurf^L,xi(ixGsurf^S,1:ndim),wC(ixGsurf^S,1:nw),need_derivs=.false.)
              ! Surface at interface xi:
              call int_surface(ixF^L,ixF^L,idims,xi(ixF^S,1:ndim),dx^DD,ones,pgeogrid%surfaceC^D,wC(ixF^S,1:nw))
              ! Surface at barycenter:
              call int_surface(ixC^L,ixC^L,idims,pgeogrid%xbar(ixC^S,1:ndim),dx^DD,ones,pgeogrid%surface^D,wC(ixC^S,1:nw))

             endif
          \}
       end select
    end do

    !--------------------------------------------------
    ! Fill the dx:
    !--------------------------------------------------
    ^D&pgeogrid%dx(ixGext^S,^D)=dx^D;

  end subroutine fillgeo_covariant
  !=============================================================================
end module mod_metric
!=============================================================================
