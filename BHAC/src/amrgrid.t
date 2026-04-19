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
subroutine settree

include 'amrvacdef.f'

! create and initialize grids on all levels > 1. On entry, all
! level=1 grids have been formed and initialized. This subroutine
! creates and initializes new level grids

integer :: igrid, iigrid, levnew
type(walloc) :: pwtmp
!----------------------------------------------------------------------------
! when only one level allowed, there is nothing to do anymore

if (mxnest == 1) return

call getbc(t,ps,psCoarse)
do levnew=2,mxnest
   if (errorestimate==1.or.errorestimate==2) then
      call setdt
      call advance(0)
   end if

   call errest

   if (errorestimate==1.or.errorestimate==2) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
!{#IFDEF DY_SP
!         !pw1(igrid)%w   => pwold(igrid)%w
!         pwtmp%w        => pwold(igrid)%w
!         pwold(igrid)%w => pw(igrid)%w
!         !pw(igrid)%w    => pw1(igrid)%w
!         pw(igrid)%w    => pwtmp%w
!}
!{#IFNDEF DY_SP
         pwtmp%w => pwold(igrid)%w
         pwold(igrid)%w => pw(igrid)%w
         pw(igrid)%w => pwtmp%w
!}
      end do
   end if

   call amr_coarsen_refine
   
   if (.not.resetgrid) then
     ! if no finer level grids created: exit
     if (levmax/=levnew) exit
   end if
end do

end subroutine settree
!=============================================================================
subroutine resettree
{#IFDEF STAGGERED
use mod_amr_fct}
{#IFDEF DY_SP
use mod_coord_cfc
}
use mod_eos, only: total_cell_prev, first_time
include 'amrvacdef.f'
integer :: iigrid, igrid, total_cell
integer :: ixGext^L
double precision :: xmin^D, dx^D
!-----------------------------------------------------------------------------

if (levmax>levmin) call deallocateBflux
{#IFDEF STAGGERED
call deallocateBfaces}

      call errest

      call amr_coarsen_refine

! set up boundary flux conservation arrays
if (levmax>levmin) call allocateBflux

end subroutine resettree
!=============================================================================
subroutine resettree_convert

include 'amrvacdef.f'
integer  :: igrid,iigrid, my_levmin, my_levmax
!-----------------------------------------------------------------------------
if (level_io > 0) then
   my_levmin = level_io
   my_levmax = level_io
else
   my_levmin = max(1,level_io_min)
   my_levmax = min(mxnest,level_io_max)
end if


do while(levmin<my_levmin.or.levmax>my_levmax)
 call getbc(t,ps,psCoarse)
 do iigrid=1,igridstail; igrid=igrids(iigrid);
    call forcedrefine_grid_io(igrid,pw(igrid)%w)
 end do

 call amr_coarsen_refine
end do
! set up boundary flux conservation arrays

end subroutine resettree_convert
!=============================================================================
