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

module mod_physicaldata
  use mod_indices, only: ngridshi
  implicit none
  save

  !=============================================================================

  type walloc
     double precision, dimension(:^D&,:), pointer :: w=>Null()
     integer                                      :: {^D& ixGmin^D=-1,ixGmax^D=-1},nwmin=-1,nwmax=-1 ! allocated index-range
     logical                                      :: allocated=.false.
  end type walloc

  ! General state datastructure:
  type state
     integer                     :: igrid=-1       ! comes in handy, with igrid, we know everything about the AMR
     integer                     :: iwpos=0        ! location of the w-array, in case its reconstructed, etc.
                                                   ! -1: unknown,
                                                   !  0: w is center (default),
                                                   !  1: w is interface 1
                                                   !  2: w is interface 2
                                                   !  3: w is interface 3
     logical                     :: w_is_primitive=.false.  ! are we in primitive state?
     logical                     :: is_coarse=.false.  ! are we a coarse buffer (psCoarse)?
     type(xalloc), pointer       :: x=>Null()      ! Center positions
     type(walloc), pointer       :: w=>Null()      ! Variables, normally center
     type(walloc), pointer       :: ws=>Null()     ! Staggered variables, always on interface
     type(walloc), pointer       :: we=>Null()     ! Edge variables, always on edge
     type(walloc), pointer       :: wc=>Null()     ! Corner variables, always on corner
     type(geoalloc), pointer     :: geo=>Null()    ! Points to the geometry for the block, geo%m is center-metric
  end type state
  
  !> allocate m1 microphysics array type
  !> store the rates: e.g. wradimpl(kappa_a_ ) ....
  !{#IFDEF M1
    type m1_walloc
      double precision, dimension(:^D&,:), pointer :: wradimpl=>Null()
      integer                                      :: {^D& ixGmin^D=-1,ixGmax^D=-1},nwm1min=-1,nwm1max=-1 !nm1rad_easmin=-1,nm1rad_easmax=-1 ! allocated index-range
      logical                                      :: allocated=.false.
    end type m1_walloc
    !-----------------
    type m1_impl_state
      logical                     :: calc_rates =.true.  !> whether to calculate rates in advect or not
      logical                     :: calc_tau_path =.true.  !> whether to calculate optical depth tau_path in advect or not
      logical                     :: calc_eqrates =.true.  !> whether to calculate equilibrium rates in advect or not
      integer                     :: igrid=-1        !> comes in handy, with igrid, we know everything about the AMR
      type(m1_walloc), pointer    :: pwrad=>Null()   !> Variables, center
    end type m1_impl_state  
  !}


  {^NOONED
  type walloc_sub
     double precision, dimension(:^DE&,:), pointer :: w=>Null()
  end type walloc_sub
  }
  {^IFONED
  type walloc_sub
     double precision, dimension(:), pointer :: w=>Null()
  end type walloc_sub
  }

  ! The centered fluid variables live here:
  type(walloc), target, dimension(ngridshi) :: pw, pwold, pw1, pw2, pw3, pw4, pwres
  type(walloc), target, dimension(ngridshi) :: pwCoarse, pwCoCo

  {#IFDEF STAGGERED
  ! The staggered fluid variables live here:
  type(walloc), target, dimension(ngridshi) :: pws, pwsold, pws1, pws2, pws3, pws4, pwsres
  type(walloc), target, dimension(ngridshi) :: pwsCoarse, pwsCoCo
  }

  ! The meta-structure state:
  type(state), dimension(ngridshi) :: ps, psold, ps1, ps2, ps3, ps4, psres
  type(state), dimension(ngridshi) :: psCoarse, psCoCo

  ! For GW backreaction
  type(walloc), target, dimension(ngridshi) :: pw_t_dot_old, pw_t_dot_old2
  type(state), dimension(ngridshi) :: ps_t_dot_old, ps_t_dot_old2

  type(walloc), dimension(ngridshi) :: pwio
  {#IFDEF STAGGERED
  type(walloc), dimension(ngridshi) :: pwsio
  }
  type(walloc), dimension(ngridshi), target :: pB0_cell,  pB0_face^D
  type(walloc), pointer :: myB0_cell=>Null(), {myB0_face^D=>Null()}, myB0=>Null()
  type(walloc_sub), dimension(ngridshi) :: pw_sub  ! For the center values on the slice
  type(walloc_sub), dimension(ngridshi) :: pwC_sub ! For the corner values on the slice
  
  ! For M1 rates, new argument of advect1()
  type(m1_walloc), target, dimension(ngridshi) :: pwm1,pwrad
  type(m1_impl_state), dimension(ngridshi) :: psm1

  {^IFONED
  double precision, dimension(:), allocatable :: collapsedData
  }
  {^NOONED
  double precision, dimension(:^DE&,:), allocatable :: collapsedData
  }

  type xalloc
     double precision, dimension(:^D&,:), pointer :: x=>Null()
     integer                                      :: {^D& ixGmin^D=-1,ixGmax^D=-1},ndimmin=-1,ndimmax=-1 ! allocated index-range
     logical                                      :: allocated=.false.
  end type xalloc

  {^NOONED
  type xalloc_sub
     double precision, dimension(:^DE&,:), pointer :: x=>Null()
  end type xalloc_sub
  }
  {^IFONED
  type xalloc_sub
     double precision, dimension(:), pointer :: x=>Null()
  end type xalloc_sub
  }
  type(xalloc), dimension(ngridshi), target :: px, pxCoarse
  type(xalloc_sub), dimension(ngridshi) :: px_sub  ! For the centerpositions on the slice
  type(xalloc_sub), dimension(ngridshi) :: pxC_sub ! For the cornerpositions on the slice

  type indexlist
     integer                                       :: i=-1
     integer                                       :: j=-1
     integer                                       :: k=-1
     double precision, dimension(:^D&), pointer    :: elem => Null() ! Data is stored here, everything else are pointers.
  end type indexlist

  type metric
     integer                                       :: nnonzero=-1, nnonzeroBeta=-1, nnonzeroDgDk=-1, nnonzeroDalphaDj=-1
     integer                                       :: nnonzeroDbetaiDj=-1
     integer                                       :: {^D& ixGmin^D=-1,ixGmax^D=-1}
     !Derivatives:
     type(indexlist), dimension(:), pointer        :: nonzeroDgDk => Null()! Non-zero elements: Metric derivatives
     type(indexlist), dimension(:), pointer        :: nonzeroDalphaDj => Null()! Non-zero elements: Lapse derivatives
     type(indexlist), dimension(:), pointer        :: nonzeroDbetaiDj => Null()! Non-zero elements: Shift derivatives
     type(indexlist), dimension(:,:,:), pointer    :: DgDk => Null()       ! Metric derivative
     type(indexlist), dimension(:,:), pointer      :: DbetaiDj => Null()   ! Shift derivative
     type(indexlist), dimension(:), pointer        :: DalphaDj => Null()   ! Lapse derivative
     !Primary elements:
     type(indexlist), dimension(:,:), pointer      :: g => Null()          ! Array of 4-metric elements
     type(indexlist), dimension(:), pointer        :: nonzero => Null()    ! Coordinatelist of non-zero elements
     type(indexlist), dimension(:), pointer        :: nonzeroBeta => Null()! Non-zero elements: Shift vector
     !Inverse of the three-metric:
     type(indexlist), dimension(:,:), pointer      :: gammainv => Null() ! Array of inverse 3-metric elements
     double precision, dimension(:^D&), pointer    :: alpha => Null()      ! lapse
     double precision, dimension(:^D&), pointer    :: zero => Null()       ! auxillary zeroes to point to
     double precision, dimension(:^D&), pointer    :: sqrtgamma => Null()  ! square root of the spatial determinant
     type(indexlist), dimension(:), pointer        :: beta => Null()       ! shift vector (contravariant)
     type(indexlist), dimension(:), pointer        :: betaD => Null()      ! shift vector (covariant, as in metric)
     type(indexlist), dimension(:,:), pointer      :: kij => Null()        ! Extrinsic curvature tensor
  end type metric

  type geoalloc
     double precision, dimension(:^D&), pointer :: dvolume=>Null()
     double precision, dimension(:^D&), pointer :: {surfaceC^D=>Null()},{surface^D=>Null()}
     double precision, dimension(:^D&,:), pointer :: dx=>Null(), xbar=>Null() ! xbar is the barycenter-position
     type(metric), pointer                        :: m=>Null(), {mSurface^D=>Null()}
  end type geoalloc

  type(geoalloc), dimension(ngridshi), target :: pgeo, pgeoCoarse, pgeoCoCo
  type(geoalloc), pointer                     :: mygeo => Null()
  type(metric), pointer                       :: myM => Null() ! points to metric on current face or center.

  ! ------------------------------ OMP directives ------------------------------
  !$OMP THREADPRIVATE(myB0_cell,myB0_face^D,myB0,mygeo,myM)
  ! ----------------------------------------------------------------------------
  
  !=============================================================================
contains
  !=============================================================================
  subroutine alloc_pw(pw,ixG^L,nwmin,nwmax)

    type(walloc)              :: pw
    integer, intent(in)       :: ixG^L, nwmin, nwmax
    ! .. local ..
    logical                   :: match
    !-----------------------------------------------------------------------------

    ! Did we previously use this structure?
    if (pw%allocated .eqv. .true.) then
       ! Check if the index-ranges match
       match = .true.
       if ({^D& pw%ixGmin^D .ne. ixGmin^D |.or.}) match = .false.
       if ({^D& pw%ixGmax^D .ne. ixGmax^D |.or.}) match = .false.
       if (pw%nwmin .ne. nwmin .or. pw%nwmax .ne. nwmax) match = .false.

       if (match) then
          ! previously allocated and indices match, return
          return
       else
          ! Indices changed, deallocate and newly allocate below
          call dealloc_pw(pw)
       end if
    end if

    {pw%ixGmin^D=ixGmin^D;}
    {pw%ixGmax^D=ixGmax^D;}
    pw%nwmin=nwmin; pw%nwmax=nwmax
    allocate(pw%w(ixG^S,nwmin:nwmax))
    pw%allocated = .true.

  end subroutine alloc_pw
  !=============================================================================
  subroutine dealloc_pw(pw)

    type(walloc)              :: pw
    !-----------------------------------------------------------------------------

    deallocate(pw%w)
    pw%allocated = .false.
    {^D& pw%ixGmin^D=-1;pw%ixGmax^D=-1|;}
    pw%nwmin=-1
    pw%nwmax=-1 

  end subroutine dealloc_pw
  !=============================================================================
  subroutine alloc_px(px,ixG^L,ndimmin,ndimmax)

    type(xalloc)              :: px
    integer, intent(in)       :: ixG^L, ndimmin, ndimmax
    ! .. local ..
    logical                   :: match
    !-----------------------------------------------------------------------------

    ! Did we previously use this structure?
    if (px%allocated .eqv. .true.) then
       ! Check if the index-ranges match
       match = .true.
       if ({^D& px%ixGmin^D .ne. ixGmin^D |.or.}) match = .false.
       if ({^D& px%ixGmax^D .ne. ixGmax^D |.or.}) match = .false.
       if (px%ndimmin .ne. ndimmin .or. px%ndimmax .ne. ndimmax) match = .false.

       if (match) then
          ! previously allocated and indices match, return
          return
       else
          ! Indices changed, deallocate and newly allocate below
          call dealloc_px(px)
       end if
    end if

    {px%ixGmin^D=ixGmin^D;}
    {px%ixGmax^D=ixGmax^D;}
    px%ndimmin=ndimmin; px%ndimmax=ndimmax
    allocate(px%x(ixG^S,ndimmin:ndimmax))
    px%allocated = .true.

  end subroutine alloc_px
  !=============================================================================
  subroutine dealloc_px(px)

    type(xalloc)              :: px
    !-----------------------------------------------------------------------------

    deallocate(px%x)
    px%allocated = .false.
    {px%ixGmin^D=-1;}
    {px%ixGmax^D=-1;}
    px%ndimmin=-1; px%ndimmax=-1

  end subroutine dealloc_px
  !=============================================================================
  subroutine alloc_state(ps,ixG^L)
    ! Allocate the solution arrays within the state structure
    ! Using info from the physics module

    !-------------------------------------------------
    ! DEFINITIONS OF GLOBAL PARAMETERS AND VARIABLES
    !-------------------------------------------------
    use amrvacpar
    INTEGER,PARAMETER:: r_=1, phi_=^PHI, z_=^Z
    INTEGER,PARAMETER:: pphi_=^PPHI, zz_=^ZZ
    
!    include 'amrvacpar.f'

    !-------------------------------------------------

    type(state)                                   :: ps
    integer, intent(in)                           :: ixG^L 
    ! .. local ..
    {#IFDEF STAGGERED
    integer                                       :: ixGs^L 
    }
    !-----------------------------------------------------------------------------

    call alloc_pw(ps%w,ixG^L,1,nw)

    {#IFDEF STAGGERED
    {^D& ixGsmin^D = ixGmin^D-1; ixGsmax^D = ixGmax^D|;}
    call alloc_pw(ps%ws,ixGs^L,1,nws)
    }

  end subroutine alloc_state
  !=============================================================================
  subroutine dealloc_state(ps)
    ! Dellocate the solution arrays within the state structure
    type(state)                                   :: ps
    !-----------------------------------------------------------------------------

    call dealloc_pw(ps%w)
    {#IFDEF STAGGERED
    call dealloc_pw(ps%ws)
    }

  end subroutine dealloc_state
  !=============================================================================
  subroutine copy_state(a,b)
    ! Copys state a to state b as in b=a

    type(state), intent(in)                        :: a
    type(state), intent(inout)                     :: b
    !-----------------------------------------------------------------------------

    b%igrid = a%igrid
    b%iwpos = a%iwpos
    b%w_is_primitive = a%w_is_primitive

    call copy_pw(a%w,b%w)
    {#IFDEF STAGGERED
    call copy_pw(a%ws,b%ws)
    }

  end subroutine copy_state
  !=============================================================================
  subroutine copy_pw(a,b)
    ! Copys walloc structure a to new walloc structure b as in b=a
    type(walloc), intent(in)                       :: a
    type(walloc), intent(inout)                    :: b
    ! .. local ..
    !-----------------------------------------------------------------------------

    call alloc_pw(b,a%ixG^L,a%nwmin,a%nwmax)

    b%w = a%w

  end subroutine copy_pw
 !=============================================================================
 {#IFDEF M1  ! radiation variables of M1 for implicit part
  !=============================================================================
   subroutine alloc_m1_pw(pwm1,ixG^L,nwm1min,nwm1max)
     type(m1_walloc)              :: pwm1
     ! nwmin = 1, nwmax = nm1rad_eas
     integer, intent(in)       :: ixG^L, nwm1min, nwm1max 
     ! .. local ..
     logical                   :: match
     !-----------------------------------------------------------------------------
     ! Did we previously use this structure?
     if (pwm1%allocated .eqv. .true.) then
        ! Check if the index-ranges match
        match = .true.
        if ({^D& pwm1%ixGmin^D .ne. ixGmin^D |.or.}) match = .false.
        if ({^D& pwm1%ixGmax^D .ne. ixGmax^D |.or.}) match = .false.
        if (pwm1%nwm1min .ne. nwm1min .or. pwm1%nwm1max .ne. nwm1max) match = .false.
        if (match) then
           ! previously allocated and indices match, return
           return
        else
           ! Indices changed, deallocate and newly allocate below
           call dealloc_m1_pw(pwm1)
        end if
     end if

     {pwm1%ixGmin^D=ixGmin^D;}
     {pwm1%ixGmax^D=ixGmax^D;}
     pwm1%nwm1min=nwm1min; pwm1%nwm1max=nwm1max

     ! allocate wradimpl !
     allocate(pwm1%wradimpl(ixG^S,nwm1min:nwm1max))
     pwm1%allocated = .true.

   end subroutine alloc_m1_pw
 !=============================================================================
   subroutine dealloc_m1_pw(pwm1)
     type(m1_walloc)              :: pwm1
     !-----------------------------------------------------------------------------
     deallocate(pwm1%wradimpl)
     pwm1%allocated = .false.
     {^D& pwm1%ixGmin^D=-1;pwm1%ixGmax^D=-1|;}
     pwm1%nwm1min=-1
     pwm1%nwm1max=-1 
   end subroutine dealloc_m1_pw
 
 !=============================================================================
   subroutine alloc_m1_state(psm1,ixG^L) 
     ! called in amr_solution_node
     ! Allocate the radiation arrays within the m1_impl_state structure
     ! Using info from the physics module
     use amrvacpar
     INTEGER,PARAMETER:: r_=1, phi_=^PHI, z_=^Z
     INTEGER,PARAMETER:: pphi_=^PPHI, zz_=^ZZ  
     type(m1_impl_state) ,intent(inout)        :: psm1
     integer, intent(in)                       :: ixG^L 
     !-----------------------------------------------------------------------------
     call alloc_m1_pw(psm1%pwrad,ixG^L,1,nm1rad_eas) !from 1 to nm1rad_eas!

   end subroutine alloc_m1_state
   !=============================================================================
   subroutine dealloc_m1_state(psm1)
     ! Dellocate the solution arrays within the state structure
     type(m1_impl_state)                 :: psm1
     call dealloc_m1_pw(psm1%pwrad)
   end subroutine dealloc_m1_state
  !=============================================================================  
 } !end IFDEF M1
 !=============================================================================  
end module mod_physicaldata
!=============================================================================
