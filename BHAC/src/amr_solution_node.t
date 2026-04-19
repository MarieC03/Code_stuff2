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
integer function getnode(ipe)
use mod_forest, only: igrid_inuse
include 'amrvacdef.f'

! getnode = get first available igrid on processor ipe

integer, intent(in) :: ipe

integer :: igrid, igrid_available
!----------------------------------------------------------------------------
igrid_available=0

do igrid=1,ngridshi
   if (igrid_inuse(igrid,ipe)) cycle

   igrid_available=igrid
   exit
end do

if (igrid_available==0) then
   write(unitterm,*) " out of nodal space - allowed ",ngridshi," grids"
   call mpistop("")
else
   getnode=igrid_available
   igrid_inuse(igrid,ipe)=.true.
end if

if (ipe==mype) then
   ! initialize nodal block
   node(1:nodehi,getnode) = 0
   rnode(1:rnodehi,getnode) = zero
end if

end function getnode
!=============================================================================
subroutine putnode(igrid,ipe)
use mod_forest
implicit none

! putnode = return igrid node on processor ipe
 
integer, intent(in) :: igrid, ipe
!----------------------------------------------------------------------------
igrid_inuse(igrid,ipe)=.false.

end subroutine putnode
!=============================================================================
subroutine alloc_node(igrid)
use mod_forest
{#IFDEF DY_SP
use mod_cfc_parameters
}
include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: level, ig^D, ixCoG^L, ixCoCoG^L, ix
double precision :: rXmin^D, dx^D
!-----------------------------------------------------------------------------

ixCoGmin^D=1;
ixCoGmax^D=ixGhi^D/2+dixB;

! set level information
level=igrid_to_node(igrid,mype)%node%level
ig^D=igrid_to_node(igrid,mype)%node%ig^D;

node(plevel_,igrid)=level
^D&node(pig^D_,igrid)=ig^D\

! set dx information
^D&rnode(rpdx^D_,igrid)=dx(^D,level)\

! determine the minimal and maximal corners
^D&rnode(rpxmin^D_,igrid)=xprobmin^D+dble(ig^D-1)*dg^D(level)\
^D&rnode(rpxmax^D_,igrid)=xprobmax^D-dble(ng^D(level)-ig^D)*dg^D(level)\

! Allocate and fill the position arrays
call alloc_px(px(igrid),ixG^LL,1,ndim)
call alloc_px(pxCoarse(igrid),ixCoG^L,1,ndim)

^D&dx^D=rnode(rpdx^D_,igrid)\
^D&rXmin^D=rnode(rpxmin^D_,igrid)-dixB*dx^D\
{do ix=ixGlo^D,ixGhi^D
   px(igrid)%x(ix^D%ixG^T,^D)=rXmin^D+(dble(ix)-half)*dx^D
end do\}
^D&dx^D=2.0d0*rnode(rpdx^D_,igrid)\
^D&rXmin^D=rnode(rpxmin^D_,igrid)-dixB*dx^D\
{do ix=ixCoGmin^D,ixCoGmax^D
   pxCoarse(igrid)%x(ix^D%ixCoG^S,^D)=rXmin^D+(dble(ix)-half)*dx^D
end do\}

{#IFDEF DY_SP
! initialize solution space
if (use_t_dependent_output .or. (initialize_metric .and. snapshotini==-1) &
    .or. restart_init_metric) then
    !.or. restart_init_metric .or. gw_br_use_I3_ij) then
  call alloc_state(psold(igrid),ixG^LL) 
endif 
}
{#IFNDEF DY_SP
if (use_t_dependent_output) then
  call alloc_state(psold(igrid),ixG^LL)
endif 
}
call alloc_state(ps(igrid),ixG^LL)
call alloc_state(psCoarse(igrid),ixCoG^L)
{#IFDEF DY_SP
if ( (initialize_metric .and. snapshotini==-1) .or. restart_init_metric) then
  call alloc_state(ps1(igrid),ixG^LL)
endif
}

{#IFDEF DY_SP
if (use_gw_br .and. gw_br_use_I3_ij) then
  ! method I, does not need
  ! method II
  ! alloc pw_t_dot_old and pw_t_dot_old2 with size of w(ixG^T, 1:6)
  ! dot of I11, I12, I13, I22, I23, I33 (representing the order of index from 1 to 6)
  ! 7: sqrt{gamma} D,  8: eps ,  9-11: vel(1:3),  12: U_br   
  !  Warning, these are not dotted quantities
  ! old = t-1 timeslice ; old2 = t-2 timeslice
  call alloc_pw(ps_t_dot_old(igrid)%w,ixG^LL,1,12)
  call alloc_pw(ps_t_dot_old2(igrid)%w,ixG^LL,1,12)
endif
}

if(residmin>smalldouble) then
   call alloc_state(psres(igrid),ixG^LL)
end if

if (errorestimate==1) then
   call mpistop('Never use errorestimate==1')
   ixCoCoGmin^D=1;
   ixCoCoGmax^D=ixCoGmax^D/2+dixB;
   call alloc_state(psCoCo(igrid),ixCoCoG^L)
end if

if (.not.slab) call getgridgeo(igrid)

if (B0field) then   
   call alloc_B0_grid(igrid)
   call set_B0_grid(igrid)
end if

!{#IFDEF M1
!  call alloc_m1_state(psm1(igrid),ixG^LL)
!}
end subroutine alloc_node
!=============================================================================
subroutine dealloc_node(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (igrid==0) then
   call mpistop("trying to delete a non-existing grid in dealloc_node")
end if

! We don't deallocate the solution arrays any more 
! but do some smart checking when trying to re-allocate
! This helps to avoid memory fragmentation.  

! We should extend this strategy to the geometry datastructures.

if (.not.slab) call putgridgeo(igrid)

if (B0field) call dealloc_B0_grid(igrid)

end subroutine dealloc_node
!=============================================================================
