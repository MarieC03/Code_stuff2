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
subroutine comm_start

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
call MPI_INIT(ierrmpi)
call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierrmpi)
call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierrmpi)

icomm=MPI_COMM_WORLD

end subroutine comm_start
!=============================================================================
subroutine comm_finalize

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
call MPI_FINALIZE(ierrmpi)

end subroutine comm_finalize
!=============================================================================
subroutine init_comm_types

include 'amrvacdef.f'

integer, dimension(ndim+1) :: sizes, subsizes, start
!integer :: i^D, ic^D, nx^D, nxCo^D, size_double
integer :: i^D, ic^D, nx^D, nxCo^D
integer(kind=MPI_ADDRESS_KIND):: size_double, lb
{#IFDEF STAGGERED integer :: idir}
!-----------------------------------------------------------------------------
nx^D=ixMhi^D-ixMlo^D+1;
nxCo^D=nx^D/2;

^D&sizes(^D)=ixGhi^D;
sizes(ndim+1)=nw
^D&subsizes(^D)=nx^D;
subsizes(ndim+1)=nwflux
^D&start(^D)=ixMlo^D-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block,ierrmpi)
call MPI_TYPE_COMMIT(type_block,ierrmpi)

^D&sizes(^D)=ixGhi^D/2+dixB;
sizes(ndim+1)=nw
^D&subsizes(^D)=nxCo^D;
subsizes(ndim+1)=nw
^D&start(^D)=ixMlo^D-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_coarse_block,ierrmpi)
call MPI_TYPE_COMMIT(type_coarse_block,ierrmpi)

^D&sizes(^D)=ixGhi^D;
sizes(ndim+1)=nw
{do ic^DB=1,2\}
   ^D&subsizes(^D)=nxCo^D;
   subsizes(ndim+1)=nw
   ^D&start(^D)=ixMlo^D-1+(ic^D-1)*nxCo^D;
   start(ndim+1)=0
   call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                                 MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                 type_sub_block(ic^D),ierrmpi)
   call MPI_TYPE_COMMIT(type_sub_block(ic^D),ierrmpi)
{end do\}

^D&sizes(^D)=ixGhi^D;
sizes(ndim+1)=nw
^D&subsizes(^D)=nx^D;
subsizes(ndim+1)=nw
^D&start(^D)=ixMlo^D-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_io,ierrmpi)

{#IFDEF STAGGERED
!----------------------------------------
! Types for staggered grid:
!----------------------------------------

^D&sizes(^D)=ixGhi^D+1;
sizes(ndim+1)=nws
^D&subsizes(^D)=nx^D+1;
subsizes(ndim+1)=nws
^D&start(^D)=ixMlo^D-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_staggered_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_staggered_io,ierrmpi)

!--------------------------------------------------------
! Types for sending coarse blocks for coarsening
^D&sizes(^D)=ixGhi^D/2+dixB+1;
sizes(ndim+1)=nws

{do ic^DB=1,2\}
do idir=1,^ND
  ^D&subsizes(^D)=nxCo^D+kr(ic^D,1)*kr(idir,^D);
  subsizes(ndim+1)=1
  ^D&start(^D)=ixMlo^D-kr(ic^D,1)*kr(idir,^D);
  start(ndim+1)=idir-1

  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                                MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                type_coarse_block_stg(idir,ic^D),ierrmpi)
  call MPI_TYPE_COMMIT(type_coarse_block_stg(idir,ic^D),ierrmpi)

end do
{end do\}

! Types for receiving for coarsening
^D&sizes(^D)=ixGhi^D+1;
sizes(ndim+1)=nws
{do ic^DB=1,2\}
do idir=1,^ND
   ^D&subsizes(^D)=nxCo^D+kr(ic^D,1)*kr(idir,^D);
   subsizes(ndim+1)=1
   ^D&start(^D)=ixMlo^D-kr(ic^D,1)*kr(idir,^D)+(ic^D-1)*nxCo^D;
   start(ndim+1)=idir-1
   call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                                 MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                 type_sub_block_stg(idir,ic^D),ierrmpi)
   call MPI_TYPE_COMMIT(type_sub_block_stg(idir,ic^D),ierrmpi)
end do
{end do\}

call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
size_block_stg_io={(nx^D+1)*}*nws*size_double

!----------------------------------------
}

!call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
size_block_io={nx^D*}*nw*size_double

{#IFDEF TRANSFORMW
if(nwtf>0) then
  ^D&sizes(^D)=ixGhi^D;
  sizes(ndim+1)=nwtf
  ^D&subsizes(^D)=nx^D;
  subsizes(ndim+1)=nwtf
  ^D&start(^D)=ixMlo^D-1;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                                MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                type_block_io_tf,ierrmpi)
  call MPI_TYPE_COMMIT(type_block_io_tf,ierrmpi)
  
  !call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
  call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
  size_block_io_tf={nx^D*}*nwtf*size_double
endif
}

^D&sizes(^D)=ixMhi^D-ixMlo^D+1;
sizes(ndim+1)=^ND
^D&subsizes(^D)=sizes(^D);
subsizes(ndim+1)=^ND
^D&start(^D)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_xcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xcc_io,ierrmpi)

^D&sizes(^D)=ixMhi^D-ixMlo^D+2;
sizes(ndim+1)=^ND
^D&subsizes(^D)=sizes(^D);
subsizes(ndim+1)=^ND
^D&start(^D)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_xc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xc_io,ierrmpi)

^D&sizes(^D)=ixMhi^D-ixMlo^D+1;
sizes(ndim+1)=nw+nwauxio
^D&subsizes(^D)=sizes(^D);
subsizes(ndim+1)=nw+nwauxio
^D&start(^D)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_wcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wcc_io,ierrmpi)

^D&sizes(^D)=ixMhi^D-ixMlo^D+2;
sizes(ndim+1)=nw+nwauxio
^D&subsizes(^D)=sizes(^D);
subsizes(ndim+1)=nw+nwauxio
^D&start(^D)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_wc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wc_io,ierrmpi)

end subroutine init_comm_types
!=============================================================================
subroutine mpistop(message)

include 'amrvacdef.f'

character(len=*), intent(in) :: message

integer :: ierrcode
!------------------------------------------------------------------------------
write(*,*) "ERROR for processor",mype,":"
write(*,*) message
call MPI_ABORT(icomm,ierrcode,ierrmpi)

end subroutine mpistop
!==============================================================================
