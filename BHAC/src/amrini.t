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
subroutine initlevelone
use mod_metric
{#IFDEF DY_SP
use mod_cfc_parameters, only: use_multigrid
}
use mod_multigrid_coupling

include 'amrvacdef.f'

integer :: iigrid, igrid
double precision :: max_global(10), max_local(10)
double precision :: min_global(1), min_local(1), lfac_max_local, lfacmax
!-----------------------------------------------------------------------------

{#IFDEF DEBUGMEE
write(*,*) "inside initlevelone"
}

levmin=1
levmax=1

call init_forest_root

call getigrids
call build_connectivity

! fill solution space of all root grids
{#IFNDEF DY_SP
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call alloc_node(igrid)
   call initial_condition(igrid)  !read the metric vars from initial data first
end do
!$OMP END PARALLEL DO
}

{#IFDEF DY_SP

!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call alloc_node(igrid)
   call initial_condition(igrid)  !read the metric vars from initial data first
enddo
!$OMP END PARALLEL DO

! When calculating the div.B, we need to multiply sqrtgamma to mygeo
call getbc(t,ps,psCoarse)
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call set_tmpGlobals(igrid)
   call get_B_field(ixG^LL,ixM^LL,ps(igrid))
enddo
!$OMP END PARALLEL DO
}

call getbc(t,ps,psCoarse)

call selectgrids 

{#IFDEF DEBUGMEE
write(*,*) "after selectgridss"
}
end subroutine initlevelone
!=============================================================================
subroutine initial_condition(igrid)

! Need only to set the mesh values (can leave ghost cells untouched)

include 'amrvacdef.f'

integer, intent(in) :: igrid
! .. local ..
integer  :: ixCoG^L, ixCoCoG^L
{#IFDEF STAGGERED
integer                                       :: ixGs^L 
}
external initonegrid_usr
!----------------------------------------------------------------------------
ixCoGmin^D=1;
ixCoGmax^D=ixGhi^D/2+dixB;


call set_tmpGlobals(igrid)

pw(igrid)%w(ixG^T,1:nw)=zero
{#IFDEF STAGGERED
{ixGsmin^D=pws(igrid)%ixGmin^D;ixGsmax^D=pws(igrid)%ixGmax^D;}
pws(igrid)%w(ixGs^S,1:nws)=zero
}

call initonegrid_usr(ixG^LL,ixM^LL,ps(igrid))

{#IFDEF DY_SP
!fill psCoarsen for fillgeo(pgeoCoarse)
!code test, do not do this?

! ixO^L = ixCoG^L will read ghostzone
!call initonegrid_usr(ixCoG^L,ixCoG^L,psCoarse(igrid))
}
end subroutine initial_condition
!=============================================================================
subroutine modify_IC

include 'amrvacdef.f'

integer :: iigrid, igrid

external initonegrid_usr
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call set_tmpGlobals(igrid)
   call initonegrid_usr(ixG^LL,ixM^LL,ps(igrid))
end do
!$OMP END PARALLEL DO

end subroutine modify_IC
!=============================================================================
