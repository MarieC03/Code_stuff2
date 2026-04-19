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

subroutine init_comm_gc(timein,psuse,psuseCo)
!! This routine initializes the ghost cell communications:
!! calculates buffer sizes, allocates buffers
!! and calls the routines 'make_task_' that place
!! all the receive requests.

use mod_comm_gc
include 'amrvacdef.f'

double precision :: timein
type(state), dimension(ngridshi), target    :: psuse, psuseCo

! ... local ...
! Counters and auxiliaries
integer :: i^D, ic^D, inc^D, idir, igrid, iigrid, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------

!! Set pointers to the state structures
time=timein
pstate => psuse
pstateCo => psuseCo

!! Define block extents. ixM and ixG from amrvacdef.f ----
ixCoGmin^D=1;
ixCoGmax^D=ixGhi^D/2+dixB;
ixCoM^L=ixCoG^L^LSUBdixB;

nx^D=ixMhi^D-ixMlo^D+1;
nxCo^D=nx^D/2;

select case (typeghostfill)
case ("copy")
   interpolation_order=1
case ("linear","unlimit")
   interpolation_order=2
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select
dixBCo=int((dixB+1)/2)

if (dixBCo+interpolation_order-1>dixB) then
   call mpistop("interpolation order for prolongation in getbc too high")
end if

!! -------------- Limits for communications -------------------
!! ... same resolution level ...
!! Cell-centred variables
{

ixS_srl_min^D(-1)=ixMlo^D
ixS_srl_max^D(-1)=ixMlo^D-1+dixB
ixS_srl_min^D(0) =ixMlo^D
ixS_srl_max^D(0) =ixMhi^D
ixS_srl_min^D(1) =ixMhi^D+1-dixB
ixS_srl_max^D(1) =ixMhi^D

ixR_srl_min^D(-1)=1
ixR_srl_max^D(-1)=dixB
ixR_srl_min^D(0) =ixMlo^D
ixR_srl_max^D(0) =ixMhi^D
ixR_srl_min^D(1) =ixMhi^D+1
ixR_srl_max^D(1) =ixGhi^D

\}
{#IFDEF STAGGERED
!! Staggered (face-allocated) variables
do idir=1,^ND
{

   ixS_srl_stg_min^D(idir,-1)=ixMlo^D-kr(idir,^D)
   ixS_srl_stg_max^D(idir,-1)=ixMlo^D-1+dixB
   ixS_srl_stg_min^D(idir,0) =ixMlo^D-kr(idir,^D)
   ixS_srl_stg_max^D(idir,0) =ixMhi^D
   ixS_srl_stg_min^D(idir,1) =ixMhi^D-dixB+1-kr(idir,^D)
   ixS_srl_stg_max^D(idir,1) =ixMhi^D
   
   ixR_srl_stg_min^D(idir,-1)=1-kr(idir,^D)
   ixR_srl_stg_max^D(idir,-1)=dixB
   ixR_srl_stg_min^D(idir,0) =ixMlo^D-kr(idir,^D)
   ixR_srl_stg_max^D(idir,0) =ixMhi^D
   ixR_srl_stg_min^D(idir,1) =ixMhi^D+1-kr(idir,^D)
   ixR_srl_stg_max^D(idir,1) =ixGhi^D

\}
end do
}

!! Sizes for srl communications
{do i^DB=-1,1\}

   ! Cell-centred variables
   sizes_srl_send(i^D)=(nwflux+nwaux)*{(ixS_srl_max^D(i^D)-ixS_srl_min^D(i^D)+1)|*}
   sizes_srl_recv(i^D)=(nwflux+nwaux)*{(ixR_srl_max^D(i^D)-ixR_srl_min^D(i^D)+1)|*}
   sizes_srl_send_total(i^D)=sizes_srl_send(i^D)
   sizes_srl_recv_total(i^D)=sizes_srl_recv(i^D)


{#IFDEF STAGGERED
   ! Staggered (face-allocated) variables   
   do idir=1,^ND
     sizes_srl_send_stg(idir,i^D)={(ixS_srl_stg_max^D(idir,i^D)-ixS_srl_stg_min^D(idir,i^D)+1)|*}
     sizes_srl_recv_stg(idir,i^D)={(ixR_srl_stg_max^D(idir,i^D)-ixR_srl_stg_min^D(idir,i^D)+1)|*}

     sizes_srl_send_total(i^D)=sizes_srl_send_total(i^D)+sizes_srl_send_stg(idir,i^D)
     sizes_srl_recv_total(i^D)=sizes_srl_recv_total(i^D)+sizes_srl_recv_stg(idir,i^D)

   end do
}

{end do\}

if (levmin/=levmax) then
!! ... restriction ...
{
   ixS_r_min^D(-1)=ixCoMmin^D
   ixS_r_min^D(0) =ixCoMmin^D
   ixS_r_min^D(1) =ixCoMmax^D+1-dixB
   ixS_r_max^D(-1)=ixCoMmin^D-1+dixB
   ixS_r_max^D(0) =ixCoMmax^D
   ixS_r_max^D(1) =ixCoMmax^D

   ixR_r_min^D(0)=1
   ixR_r_min^D(1)=ixMlo^D
   ixR_r_min^D(2)=ixMlo^D+nxCo^D
   ixR_r_min^D(3)=ixMhi^D+1
   ixR_r_max^D(0)=dixB
   ixR_r_max^D(1)=ixMlo^D-1+nxCo^D
   ixR_r_max^D(2)=ixMhi^D
   ixR_r_max^D(3)=ixGhi^D
\}

{#IFDEF STAGGERED
!! Staggered (face-allocated) variables
do idir=1,^ND
{
   ixS_r_stg_min^D(idir,-1)=ixCoMmin^D-kr(idir,^D)
   ixS_r_stg_max^D(idir,-1)=ixCoMmin^D-1+dixB
   ixS_r_stg_min^D(idir,0) =ixCoMmin^D-kr(idir,^D)
   ixS_r_stg_max^D(idir,0) =ixCoMmax^D
   ixS_r_stg_min^D(idir,1) =ixCoMmax^D+1-dixB-kr(idir,^D)
   ixS_r_stg_max^D(idir,1) =ixCoMmax^D
 
   ixR_r_stg_min^D(idir,0)=1-kr(idir,^D)
   ixR_r_stg_max^D(idir,0)=dixB
   ixR_r_stg_min^D(idir,1)=ixMlo^D-kr(idir,^D)
   ixR_r_stg_max^D(idir,1)=ixMlo^D-1+nxCo^D
   ixR_r_stg_min^D(idir,2)=ixMlo^D+nxCo^D-kr(idir,^D)
   ixR_r_stg_max^D(idir,2)=ixMhi^D
   ixR_r_stg_min^D(idir,3)=ixMhi^D+1-kr(idir,^D)
   ixR_r_stg_max^D(idir,3)=ixGhi^D

\}
end do
}

!! ... prolongation
{
   ixS_p_min^D(0)=ixMlo^D
   ixS_p_max^D(0)=ixMlo^D-1+dixBCo+(interpolation_order-1)
   ixS_p_min^D(1)=ixMlo^D
   ixS_p_max^D(1)=ixMlo^D-1+nxCo^D+dixBCo+(interpolation_order-1)
   ixS_p_min^D(2)=ixMlo^D+nxCo^D-dixBCo-(interpolation_order-1)
   ixS_p_max^D(2)=ixMhi^D
   ixS_p_min^D(3)=ixMhi^D+1-dixBCo-(interpolation_order-1)
   ixS_p_max^D(3)=ixMhi^D

   ixR_p_min^D(0)=ixCoMmin^D-dixBCo-(interpolation_order-1)
   ixR_p_max^D(0)=dixB
   ixR_p_min^D(1)=ixCoMmin^D
   ixR_p_max^D(1)=ixCoMmax^D+dixBCo+(interpolation_order-1)
   ixR_p_min^D(2)=ixCoMmin^D-dixBCo-(interpolation_order-1)
   ixR_p_max^D(2)=ixCoMmax^D
   ixR_p_min^D(3)=ixCoMmax^D+1
   ixR_p_max^D(3)=ixCoMmax^D+dixBCo+(interpolation_order-1)
\}

{#IFDEF STAGGERED
   do idir = 1, ^ND
   {if (idir.eq.^D) then
     ! Parallel components
{
     ixS_p_stg_min^D(idir,0)=ixMlo^D-1 ! -1 to make redundant 
     ixS_p_stg_max^D(idir,0)=ixMlo^D-1+dixBCo
     ixS_p_stg_min^D(idir,1)=ixMlo^D-1 ! -1 to make redundant 
     ixS_p_stg_max^D(idir,1)=ixMlo^D-1+nxCo^D+dixBCo
     ixS_p_stg_min^D(idir,2)=ixMhi^D-nxCo^D-dixBCo
     ixS_p_stg_max^D(idir,2)=ixMhi^D
     ixS_p_stg_min^D(idir,3)=ixMhi^D-dixBCo
     ixS_p_stg_max^D(idir,3)=ixMhi^D

     ixR_p_stg_min^D(idir,0)=ixCoMmin^D-1-dixBCo
     ixR_p_stg_max^D(idir,0)=ixCoMmin^D-1
     ixR_p_stg_min^D(idir,1)=ixCoMmin^D-1 ! -1 to make redundant 
     ixR_p_stg_max^D(idir,1)=ixCoMmax^D+dixBCo
     ixR_p_stg_min^D(idir,2)=ixCoMmin^D-1-dixBCo
     ixR_p_stg_max^D(idir,2)=ixCoMmax^D
     ixR_p_stg_min^D(idir,3)=ixCoMmax^D+1-1 ! -1 to make redundant 
     ixR_p_stg_max^D(idir,3)=ixCoMmax^D+dixBCo
\}
    else
{
     ! Perpendicular component
     ixS_p_stg_min^D(idir,0)=ixMlo^D
     ixS_p_stg_max^D(idir,0)=ixMlo^D-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min^D(idir,1)=ixMlo^D
     ixS_p_stg_max^D(idir,1)=ixMlo^D-1+nxCo^D+dixBCo+(interpolation_order-1)
     ixS_p_stg_min^D(idir,2)=ixMhi^D+1-nxCo^D-dixBCo-(interpolation_order-1)
     ixS_p_stg_max^D(idir,2)=ixMhi^D
     ixS_p_stg_min^D(idir,3)=ixMhi^D+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max^D(idir,3)=ixMhi^D
 
     ixR_p_stg_min^D(idir,0)=ixCoMmin^D-dixBCo-(interpolation_order-1)
     ixR_p_stg_max^D(idir,0)=ixCoMmin^D-1
     ixR_p_stg_min^D(idir,1)=ixCoMmin^D
     ixR_p_stg_max^D(idir,1)=ixCoMmax^D+dixBCo+(interpolation_order-1)
     ixR_p_stg_min^D(idir,2)=ixCoMmin^D-dixBCo-(interpolation_order-1)
     ixR_p_stg_max^D(idir,2)=ixCoMmax^D
     ixR_p_stg_min^D(idir,3)=ixCoMmax^D+1
     ixR_p_stg_max^D(idir,3)=ixCoMmax^D+dixBCo+(interpolation_order-1)
\}
    end if
   }
   end do
}

!! Sizes for multi-resolution communications
{do i^DB=-1,1\}

   ! Cell-centred variables
   sizes_r_send(i^D)=(nwflux+nwaux)*{(ixS_r_max^D(i^D)-ixS_r_min^D(i^D)+1)|*}
   sizes_r_send_total(i^D)=sizes_r_send(i^D)

{#IFDEF STAGGERED
   ! Staggered (face-allocated) variables
   do idir=1,^ND

     sizes_r_send_stg(idir,i^D)={(ixS_r_stg_max^D(idir,i^D)-ixS_r_stg_min^D(idir,i^D)+1)|*}

     sizes_r_send_total(i^D)=sizes_r_send_total(i^D)+sizes_r_send_stg(idir,i^D)

   end do
}

{end do\}

{do i^DB=0,3\}

   ! Cell-centred variables
   sizes_r_recv(i^D)=(nwflux+nwaux)*{(ixR_r_max^D(i^D)-ixR_r_min^D(i^D)+1)|*}
   sizes_p_send(i^D)=(nwflux+nwaux)*{(ixS_p_max^D(i^D)-ixS_p_min^D(i^D)+1)|*}
   sizes_p_recv(i^D)=(nwflux+nwaux)*{(ixR_p_max^D(i^D)-ixR_p_min^D(i^D)+1)|*}

   sizes_r_recv_total(i^D)=sizes_r_recv(i^D)
   sizes_p_send_total(i^D)=sizes_p_send(i^D)
   sizes_p_recv_total(i^D)=sizes_p_recv(i^D)

{#IFDEF STAGGERED
   ! Staggered (face-allocated) variables
   do idir=1,^ND

     sizes_r_recv_stg(idir,i^D)={(ixR_r_stg_max^D(idir,i^D)-ixR_r_stg_min^D(idir,i^D)+1)|*}
     sizes_p_send_stg(idir,i^D)={(ixS_p_stg_max^D(idir,i^D)-ixS_p_stg_min^D(idir,i^D)+1)|*}
     sizes_p_recv_stg(idir,i^D)={(ixR_p_stg_max^D(idir,i^D)-ixR_p_stg_min^D(idir,i^D)+1)|*}

     sizes_r_recv_total(i^D)=sizes_r_recv_total(i^D)+sizes_r_recv_stg(idir,i^D)
     sizes_p_send_total(i^D)=sizes_p_send_total(i^D)+sizes_p_send_stg(idir,i^D)
     sizes_p_recv_total(i^D)=sizes_p_recv_total(i^D)+sizes_p_recv_stg(idir,i^D)

   end do
}

{end do\}

end if !! levmin/=levmax

!! Calculate size of communication buffers
!! Loop over local grids and see what is needed
!! This could be done in build_connectivity, in connectivity.t

nrecv_gc_srl=0
nsend_gc_srl=0
nbuff_gc_send_srl=0
nbuff_gc_recv_srl=0

nrecv_gc_r=0
nsend_gc_r=0
nbuff_gc_send_r=0
nbuff_gc_recv_r=0

nrecv_gc_p=0
nsend_gc_p=0
nbuff_gc_send_p=0
nbuff_gc_recv_p=0

do iigrid=1,igridstail; igrid=igrids(iigrid);
   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle
      my_neighbor_type=neighbor_type(i^D,igrid)
      ineighbor=neighbor(1,i^D,igrid)
      ipe_neighbor=neighbor(2,i^D,igrid)

      select case(my_neighbor_type)
      case(3) !! same resolution level
        if (ipe_neighbor.ne.mype) then
        nsend_gc_srl=nsend_gc_srl+1
        nrecv_gc_srl=nrecv_gc_srl+1
        nbuff_gc_send_srl=nbuff_gc_send_srl+sizes_srl_send_total(i^D)
        nbuff_gc_recv_srl=nbuff_gc_recv_srl+sizes_srl_recv_total(i^D)
        end if
      case(2) !! Coarser neighbor
        if (ipe_neighbor.ne.mype) then
          ic^D=1+modulo(node(pig^D_,igrid)-1,2);
          ! This is the local index of the grid.
          ! Depending on it, sometimes communication does not happen.
          ! Consider for example the configuration
          !
          !       F2|
          !       --| C
          !       F1|
          !
          ! The upper corner of F1 does not need to be communicated
          ! to C because there is already F2.

          if ({(i^D==0.or.i^D==2*ic^D-3)|.and.}) then
            nsend_gc_r=nsend_gc_r+1
            nrecv_gc_p=nrecv_gc_p+1
            nbuff_gc_send_r=nbuff_gc_send_r+sizes_r_send_total(i^D)
          ! This is the local index of the prolonged ghost zone
            inc^D=ic^D+i^D;
            nbuff_gc_recv_p=nbuff_gc_recv_p+sizes_p_recv_total(inc^D)
          end if
        end if
      case(4) !! Finer neighbor
        ! Loop over the local indices of children ic^D
        ! and calculate local indices of ghost zone inc^D.
        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
           if (ipe_neighbor.ne.mype) then
             nsend_gc_p=nsend_gc_p+1
             nrecv_gc_r=nrecv_gc_r+1
             ! Although the indices change, the size of the send buffer
             ! does not depend on whether there is a pole.
             nbuff_gc_send_p=nbuff_gc_send_p+sizes_p_send_total(inc^D)
             nbuff_gc_recv_r=nbuff_gc_recv_r+sizes_r_recv_total(inc^D)
           end if
        {end do\}
      end select
   {end do\}
end do

!! Allocate communication buffers, initialize MPI requests
!! ... same resolution level ...
if (nrecv_gc_srl>0) then
   allocate(recvbuffer_srl(nbuff_gc_recv_srl))
   allocate(recvstatus_srl(MPI_STATUS_SIZE,nrecv_gc_srl),recvrequest_srl(nrecv_gc_srl))
   recvrequest_srl=MPI_REQUEST_NULL

   !! Make 'task' fill ghostcells srl

   ibuf_recv_srl=1
   irecv_srl=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_fill_gc_srl(igrid)
   end do

end if

if (nsend_gc_srl>0) then
   allocate(sendbuffer_srl(nbuff_gc_send_srl))
   allocate(sendstatus_srl(MPI_STATUS_SIZE,nsend_gc_srl),sendrequest_srl(nsend_gc_srl))
   sendrequest_srl=MPI_REQUEST_NULL

   ibuf_send_srl=1
   isend_srl=0

end if

! ... restrict ...
if (nrecv_gc_r>0) then
   allocate(recvbuffer_r(nbuff_gc_recv_r))
   allocate(recvstatus_r(MPI_STATUS_SIZE,nrecv_gc_r),recvrequest_r(nrecv_gc_r))
   recvrequest_r=MPI_REQUEST_NULL

   !! Make 'task' fill ghostcells srl

   ibuf_recv_r=1
   irecv_r=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_fill_gc_r(igrid)
   end do

end if

if (nsend_gc_r>0) then
   allocate(sendbuffer_r(nbuff_gc_send_r))
   allocate(sendstatus_r(MPI_STATUS_SIZE,nsend_gc_r),sendrequest_r(nsend_gc_r))
   sendrequest_r=MPI_REQUEST_NULL

   ibuf_send_r=1
   isend_r=0

end if

! ... prolong ...
if (nrecv_gc_p>0) then
   allocate(recvbuffer_p(nbuff_gc_recv_p))
   allocate(recvstatus_p(MPI_STATUS_SIZE,nrecv_gc_p),recvrequest_p(nrecv_gc_p))
   recvrequest_p=MPI_REQUEST_NULL

   !! Make 'task' fill ghostcells srl

   ibuf_recv_p=1
   irecv_p=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_fill_gc_p(igrid)
   end do

end if

if (nsend_gc_p>0) then
   allocate(sendbuffer_p(nbuff_gc_send_p))
   allocate(sendstatus_p(MPI_STATUS_SIZE,nsend_gc_p),sendrequest_p(nsend_gc_p))
   sendrequest_p=MPI_REQUEST_NULL

   ibuf_send_p=1
   isend_p=0

end if


{^IFPHI
! Allocate buffers for reversing phi at poles
allocate(pole_buf%w(ixG^T,1:nwflux+nwaux))
{#IFDEF STAGGERED
allocate(pole_buf_stg%w(ixGlo^D-1:ixGhi^D,1:^ND))
}
}


end subroutine init_comm_gc
!===============================================================================
subroutine make_task_fill_gc_srl(igrid)
use mod_comm_gc
include 'amrvacdef.f'

! In a loop, place MPI receive requests, updating the part of the receive
! buffer where ghost cells will be written
! Note: ibuf_recv and irecv_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i^D, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

{do i^DB=-1,1\}
   if (i^D==0|.and.) cycle
   my_neighbor_type=neighbor_type(i^D,igrid)
   select case (my_neighbor_type) !! Could the restrict receive be also here?
   case (3) !!! Case 3 = srl
     ipe_neighbor=neighbor(2,i^D,igrid)
       if (ipe_neighbor/=mype) then
         irecv_srl=irecv_srl+1
         itag=(3**^ND+4**^ND)*(igrid-1)+{(i^D+1)*3**(^D-1)+}
         call MPI_IRECV(recvbuffer_srl(ibuf_recv_srl),sizes_srl_recv_total(i^D), &
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                        icomm,recvrequest_srl(irecv_srl),ierrmpi)
         !call add_task_to_list() 
         ibuf_recv_srl=ibuf_recv_srl+sizes_srl_recv_total(i^D)
       end if
   end select
{end do\}

end subroutine make_task_fill_gc_srl
!===============================================================================
subroutine fill_gc_srl
use mod_comm_gc
include 'amrvacdef.f'

integer :: i^D, igrid, iigrid
integer :: ipole, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
! Wait for the receive buffer to be complete

if (nrecv_gc_srl>0) then
   call MPI_WAITALL(nrecv_gc_srl,recvrequest_srl,recvstatus_srl,ierrmpi)
   ibuf_recv_srl=1
end if
if (nsend_gc_srl>0) then
   call MPI_WAITALL(nsend_gc_srl,sendrequest_srl,sendstatus_srl,ierrmpi)
end if

! In a loop with the same order as that in make_task_fill_gc_srl,
! directly fill the ghost cells that are in the same cpu and 
! unpack the receive buffer to fill those that are not.

do iigrid=1,igridstail; igrid=igrids(iigrid);
   call set_tmpGlobals(igrid)
      
  {do i^DB=-1,1\}
     if (i^D==0|.and.) cycle
     my_neighbor_type=neighbor_type(i^D,igrid)
     select case (my_neighbor_type) !! Could the restrict receive be also here?
     case (3) !!! Case 3 = srl
        call bc_fill_srl
     end select
  {end do\}
end do

! Deallocate communication buffers
if (nrecv_gc_srl>0) deallocate(recvbuffer_srl,recvstatus_srl,recvrequest_srl)
if (nsend_gc_srl>0) deallocate(sendbuffer_srl,sendstatus_srl,sendrequest_srl)

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_fill_srl
integer :: ixS^L,ixR^L,n_i^D,ixSsync^L,ixRsync^L
integer :: ibufnext,idir
{#IFDEF STAGGERED
double precision, dimension(pstate(igrid)%ws%ixG^S) :: tmp
}
integer :: idirect 
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i^D,igrid)
ipe_neighbor=neighbor(2,i^D,igrid)
ipole=neighbor_pole(i^D,igrid)
idirect={abs(i^D)|+}


!! Now the special treatment of the pole is done here, at the receive step
if (ipole.eq.0) then    
  n_i^D=-i^D;
  ixR^L=ixR_srl_^L(i^D);
  ixS^L=ixS_srl_^L(n_i^D);

  if (ipe_neighbor.eq.mype) then
    !! Just copy from the other block
      pstate(igrid)%w%w(ixR^S,1:nwflux+nwaux)=pstate(ineighbor)%w%w(ixS^S,1:nwflux+nwaux)
   {#IFDEF STAGGERED
   do idir=1,^ND
      ixS^L=ixS_srl_stg_^L(idir,n_i^D);
      ixR^L=ixR_srl_stg_^L(idir,i^D);
      
      ! opedit: in principle, averaging together would be desirable.
      ! when setting pstate(ineighbor) = pstate(igrid) here this leads to a bug though, which is
      ! not fully understood.  
      ! So for now we just exchange twice, meaning we keep either the left or the right representation
      ! of the staggerd field between blocks. In any case, the states are sychronized this way.
!      if (idirect .eq. 1) then
!         call indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L) ! Overwrites ixR, ixS
!         pstate(igrid)%ws%w(ixRsync^S,idir) = &
!              half*(pstate(igrid)%ws%w(ixRsync^S,idir) + pstate(ineighbor)%ws%w(ixSsync^S,idir))
!      end if
      
      pstate(igrid)%ws%w(ixR^S,idir) = pstate(ineighbor)%ws%w(ixS^S,idir)
   end do
   }
  else
    !! Unpack the buffer and fill the ghost cells
    ibufnext=ibuf_recv_srl+sizes_srl_recv(i^D)
    pstate(igrid)%w%w(ixR^S,1:nwflux+nwaux)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibufnext-1),shape=shape(pstate(igrid)%w%w(ixR^S,1:nwflux+nwaux)))
  
    ibuf_recv_srl=ibufnext
  
    {#IFDEF STAGGERED
    do idir=1,^ND
       ixS^L=ixS_srl_stg_^L(idir,n_i^D);
       ixR^L=ixR_srl_stg_^L(idir,i^D);

       ibufnext=ibuf_recv_srl+sizes_srl_recv_stg(idir,i^D)
       
       tmp(ixS^S) = reshape(source=recvbuffer_srl(ibuf_recv_srl:ibufnext-1),shape=shape(pstate(igrid)%ws%w(ixS^S,idir)))       

       if (idirect .eq. 1) then
          call indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L) ! Overwrites ixR, ixS
          pstate(igrid)%ws%w(ixRsync^S,idir) = &
               half*(tmp(ixSsync^S) + pstate(igrid)%ws%w(ixRsync^S,idir))
       end if

       pstate(igrid)%ws%w(ixR^S,idir) = tmp(ixS^S)
    
      ibuf_recv_srl=ibufnext
   end do
    }
  end if

else ! There is a pole
  ixR^L=ixR_srl_^L(i^D);
  select case (ipole)
  {case (^D)
     n_i^D=i^D^D%n_i^DD=-i^DD;\}
  end select
  ixS^L=ixS_srl_^L(n_i^D);

  if (ipe_neighbor==mype) then
    !! Fill ghost cells
    call pole_copy(pstate(igrid)%w,ixR^L,pstate(ineighbor)%w,ixS^L,ipole,i^D)

    {#IFDEF STAGGERED
    do idir=1,^ND
       ixR^L=ixR_srl_stg_^L(idir,i^D);
       ixS^L=ixS_srl_stg_^L(idir,n_i^D);
       !! Fill ghost cells
       call pole_copy_stg(pstate(igrid)%ws,ixR^L,pstate(ineighbor)%ws,ixS^L,ipole,i^D,idir)
    end do
    }

  else
    !! Unpack the buffer and fill an auxiliary array
    ibufnext=ibuf_recv_srl+sizes_srl_recv(i^D)
    pole_buf%w=zero
    pole_buf%w(ixS^S,1:nwflux+nwaux)=&
      reshape(source=recvbuffer_srl(ibuf_recv_srl:ibufnext-1),&
             shape=shape(pstate(igrid)%w%w(ixS^S,1:nwflux+nwaux)))
    ibuf_recv_srl=ibufnext
 
    !! Fill ghost cells
    call pole_copy(pstate(igrid)%w,ixR^L,pole_buf,ixS^L,ipole,i^D)

    {#IFDEF STAGGERED
    pole_buf_stg%w=zero
    do idir=1,^ND
       ixR^L=ixR_srl_stg_^L(idir,i^D);
       ixS^L=ixS_srl_stg_^L(idir,n_i^D);

       ibufnext=ibuf_recv_srl+sizes_srl_recv_stg(idir,i^D)
       pole_buf_stg%w(ixS^S,idir)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibufnext-1),&
         shape=shape(pstate(igrid)%ws%w(ixS^S,idir)))
       ibuf_recv_srl=ibufnext

       call pole_copy_stg(pstate(igrid)%ws,ixR^L,pole_buf_stg,ixS^L,ipole,i^D,idir)

    end do
    }


  end if
end if

end subroutine bc_fill_srl
!=============================================================================
subroutine indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L)
  integer, intent(in)       :: i^D,idir
  integer, intent(inout)    :: ixR^L,ixS^L
  integer, intent(out)      :: ixRsync^L,ixSsync^L
  ! .. local ..
  !-----------------------------------------------------------------------------

  ixRsync^L=ixR^L;
  ixSsync^L=ixS^L;
  
  {
  if (i^D .eq. -1 .and. idir .eq. ^D) then
     ixRsyncmin^D = ixRmax^D
     ixRsyncmax^D = ixRmax^D
     ixSsyncmin^D = ixSmax^D
     ixSsyncmax^D = ixSmax^D
     ixRmax^D = ixRmax^D - 1
     ixSmax^D = ixSmax^D - 1
  else if (i^D .eq. 1 .and. idir .eq. ^D) then
     ixRsyncmin^D = ixRmin^D
     ixRsyncmax^D = ixRmin^D
     ixSsyncmin^D = ixSmin^D
     ixSsyncmax^D = ixSmin^D
     ixRmin^D = ixRmin^D + 1
     ixSmin^D = ixSmin^D + 1
  end if
  \}

end subroutine indices_for_syncing
!=============================================================================
end subroutine fill_gc_srl
!===============================================================================
subroutine send_gc_srl(igrid)
use mod_comm_gc
include 'amrvacdef.f'
! ibuf_send and isend_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer :: igrid, i^D, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------

call set_tmpGlobals(igrid)

{do i^DB=-1,1\}
   if (i^D==0|.and.) cycle
   my_neighbor_type=neighbor_type(i^D,igrid)
   if (my_neighbor_type==3) call bc_send_srl
{end do\}



contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_send_srl
integer :: n_i^D,idir,ixS^L
integer :: ibufaux,ibufnext,ipole
 ! Auxialiary array to avoid problems with the preprocessor...
 ! (It keeps splitting the line in the wrong place)
integer, dimension(1) :: sizes
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i^D,igrid)
ipe_neighbor=neighbor(2,i^D,igrid)
ipole=neighbor_pole(i^D,igrid)

!If the neighbor is in another cpu, ...
if (ipe_neighbor.ne.mype) then

  ! The ghost region of the neighbor changes
  ! if there is a pole
  select case (ipole)
  case(0) ! No pole
     n_i^D=-i^D;
 {case (^D) ! Pole in this direction
     n_i^D=i^D^D%n_i^DD=-i^DD;\}
  end select

  ! fill corresponding part of the send buffer...

  ixS^L=ixS_srl_^L(i^D);
  ibufaux=ibuf_send_srl
  ibufnext=ibufaux+sizes_srl_send(i^D)

  sizes=(/sizes_srl_send(i^D)/)

  sendbuffer_srl(ibufaux:ibufnext-1)=reshape(pstate(igrid)%w%w(ixS^S,1:nwflux+nwaux),sizes)

  ibufaux=ibufnext
{#IFDEF STAGGERED
  do idir=1,^ND
    ixS^L=ixS_srl_stg_^L(idir,i^D);
    ibufnext=ibufaux+sizes_srl_send_stg(idir,i^D)
    sizes=(/sizes_srl_send_stg(idir,i^D)/)
    sendbuffer_srl(ibufaux:ibufnext-1)=&
      reshape(pstate(igrid)%ws%w(ixS^S,idir),sizes)   
    ibufaux=ibufnext
  end do
}

  ! ...and send
  itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
  isend_srl=isend_srl+1
  call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),sizes_srl_send_total(i^D), &
                 MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                 icomm,sendrequest_srl(isend_srl),ierrmpi)

  ibuf_send_srl=ibufnext

end if
end subroutine bc_send_srl
end subroutine send_gc_srl
!===============================================================================
subroutine make_task_fill_gc_r(igrid)
use mod_comm_gc
include 'amrvacdef.f'

! In a loop, place MPI receive requests, updating the part of the receive
! buffer where ghost cells will be written
! Note: ibuf_recv_r and irecv_r were already initialized in init_comm_gc,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i^D, ic^D, inc^D, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

{do i^DB=-1,1\}
   if (i^D==0|.and.) cycle
   my_neighbor_type=neighbor_type(i^D,igrid)
   if (my_neighbor_type.ne.4) cycle
   ! Loop over the local indices of children ic^D
   ! and calculate local indices of ghost zone inc^D.
   {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
      inc^DB=2*i^DB+ic^DB\}
      ipe_neighbor=neighbor_child(2,inc^D,igrid)
      if (ipe_neighbor.ne.mype) then
        irecv_r=irecv_r+1
        itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
        call MPI_IRECV(recvbuffer_r(ibuf_recv_r),sizes_r_recv_total(inc^D), &
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                       icomm,recvrequest_r(irecv_r),ierrmpi)
        ibuf_recv_r=ibuf_recv_r+sizes_r_recv_total(inc^D)
      end if
   {end do\}
{end do\}

end subroutine make_task_fill_gc_r
!===============================================================================
subroutine fill_gc_r
use mod_comm_gc
include 'amrvacdef.f'

integer :: i^D, igrid, iigrid
integer :: ipole, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
! Wait for the receive buffer to be complete

if (nrecv_gc_r>0) then
   call MPI_WAITALL(nrecv_gc_r,recvrequest_r,recvstatus_r,ierrmpi)
   ibuf_recv_r=1
end if

! In a loop with the same order as that in make_task_fill_gc_r,
! directly fill the ghost cells that are in the same cpu and 
! unpack the receive buffer to fill those that are not.

do iigrid=1,igridstail; igrid=igrids(iigrid);
  call set_tmpGlobals(igrid)
     
 {do i^DB=-1,1\}
    if (i^D==0|.and.) cycle
    my_neighbor_type=neighbor_type(i^D,igrid)
    if (my_neighbor_type.eq.4) then
      call bc_fill_r
    end if
 {end do\}
end do


! Wait for the sends to complete and deallocate communication buffers
if (nrecv_gc_r>0) deallocate(recvbuffer_r,recvstatus_r,recvrequest_r)

if (nsend_gc_r>0) then
   call MPI_WAITALL(nsend_gc_r,sendrequest_r,sendstatus_r,ierrmpi)
   deallocate(sendbuffer_r,sendstatus_r,sendrequest_r)
end if

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_fill_r
integer :: ixS^L,ixR^L,ic^D,inc^D,n_i^D
integer :: ibufnext,idir
!-----------------------------------------------------------------------------
ipole=neighbor_pole(i^D,igrid)

!write(*,*) 'passed bc_fill r'

! Check first if there is pole
if (ipole.eq.0) then
! Loop over the children ic^D to and their neighbors inc^D
{do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
   inc^DB=2*i^DB+ic^DB\}
   n_i^D=-i^D;
   ineighbor=neighbor_child(1,inc^D,igrid)
   ipe_neighbor=neighbor_child(2,inc^D,igrid)
   if (ipe_neighbor.eq.mype) then ! Same processor
     ixS^L=ixS_r_^L(n_i^D);
     ixR^L=ixR_r_^L(inc^D);
     pstate(igrid)%w%w(ixR^S,1:nwflux+nwaux)=pstateCo(ineighbor)%w%w(ixS^S,1:nwflux+nwaux)


   {#IFDEF STAGGERED
     do idir=1,^ND
        ixS^L=ixS_r_stg_^L(idir,n_i^D);
        ixR^L=ixR_r_stg_^L(idir,inc^D);
        pstate(igrid)%ws%w(ixR^S,idir)=pstateCo(ineighbor)%ws%w(ixS^S,idir)
    !    print *, 'idir',idir,'n_i',n_i^D
    !    print *, 'ixS',ixS^L
    !    print *, 'idir',idir,'inc',inc^D
    !    print *, 'ixR',ixR^L
    !    call printarray(ixR^L,ixR^L,pstate(igrid)%ws%w(ixR^S,idir))
     end do
   }

   else ! Different processor
     ixR^L=ixR_r_^L(inc^D);
     !! Unpack the buffer and fill the ghost cells
     ibufnext=ibuf_recv_r+sizes_r_recv(inc^D)
     pstate(igrid)%w%w(ixR^S,1:nwflux+nwaux)=reshape(source=recvbuffer_r(ibuf_recv_r:ibufnext-1),shape=shape(pstate(igrid)%w%w(ixR^S,1:nwflux+nwaux)))
   
     ibuf_recv_r=ibufnext

    {#IFDEF STAGGERED
     do idir=1,^ND
       ixR^L=ixR_r_stg_^L(idir,inc^D);
       ibufnext=ibuf_recv_r+sizes_r_recv_stg(idir,inc^D)
       pstate(igrid)%ws%w(ixR^S,idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibufnext-1),shape=shape(pstate(igrid)%ws%w(ixR^S,idir)))
       ibuf_recv_r=ibufnext
     end do
    }
   end if
{end do\}

else !! There is a pole

{do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
   inc^DB=2*i^DB+ic^DB\}
   select case(ipole)
  {case (^D)
     n_i^D=i^D^D%n_i^DD=-i^DD;\}
   end select
   ineighbor=neighbor_child(1,inc^D,igrid)
   ipe_neighbor=neighbor_child(2,inc^D,igrid)
   if (ipe_neighbor.eq.mype) then ! Same processor
     ixS^L=ixS_r_^L(n_i^D);
     ixR^L=ixR_r_^L(inc^D);
     !! Fill ghost cells
     call pole_copy(pstate(igrid)%w,ixR^L,pstateCo(ineighbor)%w,ixS^L,ipole,i^D)

   {#IFDEF STAGGERED
     do idir=1,^ND
        ixS^L=ixS_r_stg_^L(idir,n_i^D);
        ixR^L=ixR_r_stg_^L(idir,inc^D);
        !! Fill ghost cells
        call pole_copy_stg(pstate(igrid)%ws,ixR^L,pstateCo(ineighbor)%ws,ixS^L,ipole,i^D,idir)

     end do
   }
   else ! Different processor
     !! Unpack the buffer and fill an auxiliary array
     ixS^L=ixS_r_^L(n_i^D);
     ixR^L=ixR_r_^L(inc^D);

     ibufnext=ibuf_recv_r+sizes_r_recv(inc^D)
     pole_buf%w=zero
     pole_buf%w(ixR^S,1:nwflux+nwaux)=&
       reshape(source=recvbuffer_r(ibuf_recv_r:ibufnext-1),&
              shape=shape(pstate(igrid)%w%w(ixR^S,1:nwflux+nwaux)))
     ibuf_recv_r=ibufnext
 
     !! Fill ghost cells
     call pole_copy(pstate(igrid)%w,ixR^L,pole_buf,ixR^L,ipole,i^D)

     {#IFDEF STAGGERED
     pole_buf_stg%w=zero
     do idir=1,^ND
        ixS^L=ixS_r_stg_^L(idir,n_i^D);
        ixR^L=ixR_r_stg_^L(idir,inc^D);

        ibufnext=ibuf_recv_r+sizes_r_recv_stg(idir,inc^D)
        pole_buf_stg%w(ixR^S,idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibufnext-1),&
          shape=shape(pstate(igrid)%ws%w(ixR^S,idir)))
        ibuf_recv_r=ibufnext

        call pole_copy_stg(pstate(igrid)%ws,ixR^L,pole_buf_stg,ixR^L,ipole,i^D,idir)

     end do
     }


   end if
{end do\}

end if !! ipole == 0


end subroutine bc_fill_r
end subroutine fill_gc_r
!===============================================================================
subroutine gc_restrict(igrid)
! This subroutine fills the coarse representation of the given block.
! 
use mod_comm_gc
include 'amrvacdef.f'
! ibuf_send and isend_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i^D, ineighbor, ipe_neighbor, my_neighbor_type, ixCoR^L, ixR^L
!-------------------------------------------------------------------------------

associate(x=>pstate(igrid)%x%x,xCo=>pstateCo(igrid)%x%x,&
          pgeoFi=>pstate(igrid)%geo,pgeoCo=>pstateCo(igrid)%geo)

call set_tmpGlobals(igrid)

if (any(neighbor_type(:^D&,igrid).eq. 2)) then

   {do i^DB=-1,1\}
      my_neighbor_type=neighbor_type(i^D,igrid)      
      ! Restriction is necessary in the physical extent
      ! for sending and in srl ghost regions for later
      ! prolongation
      if (my_neighbor_type .eq. 0 .or. &
           my_neighbor_type .eq. 3) then
      
         ixR^L=ixR_srl_^L(i^D);
         {ixCoRmin^D=int((ixRmin^D+dixB+1)/2);}
         {ixCoRmax^D=int((ixRmax^D+dixB+1)/2);}
         call coarsen_grid(pstate(igrid),x,ixG^LL,ixR^L,pstateCo(igrid),xCo, &
              ixCoG^L,ixCoR^L,pgeoFi,pgeoCo, &
              coarsenprimitive,.true.)

      end if
   {end do\}
end if

end associate

end subroutine gc_restrict
!===============================================================================
subroutine gc_prolong(igrid)
use mod_comm_gc
include 'amrvacdef.f'

integer, intent(in) :: igrid
! ... local ...
integer :: i^D,iside,idims,my_neighbor_type
{#IFDEF STAGGERED
logical,dimension(-1:1^D&)   :: NeedProlong
integer                      :: idim
}
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

{#IFNDEF STAGGERED
if (any(neighbor_type(:^D&,igrid)==2)) then
   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle
      my_neighbor_type=neighbor_type(i^D,igrid)
      if (my_neighbor_type==2) then
         call bc_prolong
        ! NeedProlong(i^D)=.true.
      end if
   {end do\}
end if
}
{#IFDEF STAGGERED
if (any(neighbor_type(:^D&,igrid)==2)) then
   ! Fill the array NeedProlong
   ! and prolong cell-centered variables
   NeedProlong=.false.

   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle
      my_neighbor_type=neighbor_type(i^D,igrid)
      if (my_neighbor_type==2) then
         call bc_prolong
         NeedProlong(i^D)=.true.
      end if
   {end do\}


   ! Ghost cell prolongation for staggered variables
   ! must be done in a specific order.
   ! First the first neighbours, which have 2 indices=0 in 3D
   ! or one index=0 in 2D
   do idim=1,^ND
     i^D=0;
     select case(idim)
    {case(^D)
       do i^D=-1,1,2
         if (NeedProlong(i^DD)) call bc_prolong_stg(NeedProlong)
       end do
     \}
     end select
   end do
   ! Then the second neighbours which have 1 index=0 in 3D
   ! (Only in 3D)
   {^IFTHREED
   i1=0;
   do i2=-1,1,2
   do i3=-1,1,2
   if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
   end do
   end do
   i2=0;
   do i3=-1,1,2
   do i1=-1,1,2
   if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
   end do
   end do
   i3=0;
   do i1=-1,1,2
   do i2=-1,1,2
   if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
   end do
   end do
   }
   ! Finally, the corners, that have no index=0
  {do i^D=-1,1,2\}
   if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
 {end do\}

end if
}


contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_prolong
integer :: ixFi^L,ixCo^L,ii^D
double precision :: dxFi^D, dxCo^D, xFimin^D, xComin^D, invdxCo^D
integer :: ixB^L
logical :: skip_primitive(ixG^T)
!--- for debugging ----
integer :: idx^D
!-----------------------------------------------------------------------------
associate(pgeoCo=>pstateCo(igrid)%geo)



ixFi^L=ixR_srl_^L(i^D);

dxFi^D=rnode(rpdx^D_,igrid);
dxCo^D=two*dxFi^D;
invdxCo^D=1.d0/dxCo^D;

xFimin^D=rnode(rpxmin^D_,igrid)-dble(dixB)*dxFi^D;
xComin^D=rnode(rpxmin^D_,igrid)-dble(dixB)*dxCo^D;

! moved the physical boundary filling here, to only fill the
! part needed

ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;


do idims=1,ndim
   do iside=1,2
      ii^D=kr(^D,idims)*(2*iside-3);

      if (neighbor_type(ii^D,igrid)/=1) cycle

      if  (( {(iside==1.and.idims==^D.and.ixComin^D<ixCoGmin^D+dixB)|.or.} ) &
       .or.( {(iside==2.and.idims==^D.and.ixComax^D>ixCoGmax^D-dixB)|.or. }))then
        {ixBmin^D=merge(ixCoGmin^D,ixComin^D,idims==^D);}
        {ixBmax^D=merge(ixCoGmax^D,ixComax^D,idims==^D);}
        if(.not.slab)mygeo=>pgeoCo
        {#IFNDEF DY_SP
        if(covariant)myM => mygeo%m
        }

        call bc_phys(iside,idims,time,pstateCo(igrid),ixB^L)
      end if
   end do
end do

if(.not.slab)mygeo=>pgeoCo
{#IFNDEF DY_SP
if(covariant)myM => mygeo%m
}
if (amrentropy) then
   call e_to_rhos(ixCoG^L,ixCo^L,pstateCo(igrid)%w%w,pxCoarse(igrid)%x)
else if (prolongprimitive) then
   ! Convert to primitives only in the region with meaningful values.
   skip_primitive=.true.
{  do idx^DB=ixComin^DB,ixComax^DB\}
{^IFTWOD
      if (((ixComin1+2-interpolation_order.lt.idx1).and.(idx1.lt.ixComax1-2+interpolation_order)).or.& 
          ((ixComin2+2-interpolation_order.lt.idx2).and.(idx2.lt.ixComax2-2+interpolation_order))) then 
}
{^IFTHREED
      if ((((ixComin1+2-interpolation_order.lt.idx1).and.(idx1.lt.ixComax1-2+interpolation_order)).and.&
           ((ixComin2+2-interpolation_order.lt.idx2).and.(idx2.lt.ixComax2-2+interpolation_order))).or.&
          (((ixComin2+2-interpolation_order.lt.idx2).and.(idx2.lt.ixComax2-2+interpolation_order)).and.&
           ((ixComin3+2-interpolation_order.lt.idx3).and.(idx3.lt.ixComax3-2+interpolation_order))).or.&
          (((ixComin3+2-interpolation_order.lt.idx3).and.(idx3.lt.ixComax3-2+interpolation_order)).and.&
           ((ixComin1+2-interpolation_order.lt.idx1).and.(idx1.lt.ixComax1-2+interpolation_order)))) then
}
      skip_primitive(idx^D)=.false.
      !call primitive(ixCoG^L,idx^D,idx^D,pstateCo(igrid)%w%w,pxCoarse(igrid)%x)

{^NOONED
end if
}
{  end do\}
end if

!not sure add this: 
! amrvac3.0 only update faces2center when using prolongprimitive
! gmunu update it when not prolongprim, with a call of p2c
! bhac does not do anything
!if (.not. prolongprimitive) then
!   ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
!   ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;
!   ! fixme: if added p2c here, --> rho>rhomax
!   call conserven(ixCoG^L,ixCo^L,pstateCo(igrid)%w%w,pxCoarse(igrid)%x,patchfalse)
!   {#IFDEF STAGGERED
!   call faces2centers(ixCo^L,pstateCo(igrid))
!   }
!endif

select case (typeghostfill)
case ("linear")
   call interpolation_linear(pstate(igrid),ixFi^L,dxFi^D,xFimin^D, &
                           pstateCo(igrid),ixCo^L,dxCo^D,invdxCo^D,xComin^D)
case ("copy")
   call interpolation_copy(pstate(igrid)%w,ixFi^L,dxFi^D,xFimin^D, &
                           pstateCo(igrid)%w,dxCo^D,invdxCo^D,xComin^D)
case ("unlimit")
   call interpolation_unlimit(pstate(igrid)%w,ixFi^L,dxFi^D,xFimin^D, &
                           pstateCo(igrid)%w,dxCo^D,invdxCo^D,xComin^D)
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select

if(.not.slab)mygeo=>pgeoCo
{#IFNDEF DY_SP
if(covariant)myM => mygeo%m
}
! code test erase first
{#IFNDEF DY_SP
if (amrentropy) then
    call rhos_to_e(ixCoG^L,ixCo^L,pstateCo(igrid)%w%w,pxCoarse(igrid)%x)
else if (prolongprimitive) then
    call conserve(ixCoG^L,ixCo^L,pstateCo(igrid)%w%w,pxCoarse(igrid)%x,skip_primitive)
end if
}
end associate
end subroutine bc_prolong
{#IFDEF STAGGERED
!=============================================================================
subroutine bc_prolong_stg(NeedProlong)
use mod_amr_fct
logical,dimension(-1:1^D&) :: NeedProlong
!... local ...
logical                    :: fine_^Lin
integer                    :: ixFi^L,ixCo^L
double precision           :: dxFi^D,dxCo^D,xFimin^D,xComin^D,invdxCo^D
!-----------------------------------------------------------------------------
! Check what is already at the desired level
fine_^Lin=.false.;
{
if(i^D.gt.-1) fine_min^Din=(.not.NeedProlong(i^DD-kr(^D,^DD)).and.neighbor_type(i^DD-kr(^D,^DD),igrid)/=1)
if(i^D.lt.1)  fine_max^Din=(.not.NeedProlong(i^DD+kr(^D,^DD)).and.neighbor_type(i^DD+kr(^D,^DD),igrid)/=1)
\}

ixFi^L=ixR_srl_^L(i^D);

dxFi^D=rnode(rpdx^D_,igrid);
dxCo^D=two*dxFi^D;
invdxCo^D=1.d0/dxCo^D;

xFimin^D=rnode(rpxmin^D_,igrid)-dble(dixB)*dxFi^D;
xComin^D=rnode(rpxmin^D_,igrid)-dble(dixB)*dxCo^D;

! moved the physical boundary filling here, to only fill the
! part needed

ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;

call prolong_2nd_stg(pstateCo(igrid),pstate(igrid),ixCo^L,ixFi^L,dxCo^D,xComin^D,dxFi^D,xFimin^D,.true.,fine_^Lin)

! The current region has already been refined, so it doesn t need to be prolonged again

 NeedProlong(i^D)=.false. 

end subroutine bc_prolong_stg
}
!=============================================================================
!=== Prolongation routines (Should be made independent and moved elsewhere) ==
!=============================================================================
subroutine interpolation_linear(psFi,ixFi^L,dxFi^D,xFimin^D, &
                                psCo,ixCo^L,dxCo^D,invdxCo^D,xComin^D)

integer, intent(in) :: ixFi^L
double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D
type(state) :: psCo, psFi
type(walloc) :: pwCo, pwFi
integer :: ixCo^L,ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, iw, idims
double precision :: xCo^D, xFi^D, eta^D
double precision :: slopeL, slopeR, slopeC, signC, signR, slope(nwflux+nwaux,ndim)
!-----------------------------------------------------------------------------
associate(pgeoFi=>pstate(igrid)%geo)

pwCo=psCo%w
pwFi=psFi%w
{do ixFi^DB = ixFi^LIM^DB
   ! cell-centered coordinates of fine grid point
   xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB

   ! indices of coarse cell which contains the fine cell
   ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1

   ! cell-centered coordinate for coarse cell
   xCo^DB=xComin^DB+(dble(ixCo^DB)-half)*dxCo^DB\}

   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta^D=(xFi^D-xCo^D)*invdxCo^D;
   else
      ix^D=2*int((ixFi^D+ixMlo^D)/2)-ixMlo^D;
      {eta^D=(xFi^D-xCo^D)*invdxCo^D &
            *two*(one-pgeoFi%dvolume(ixFi^DD) &
            /sum(pgeoFi%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
   end if

   
   do idims=1,ndim
      hxCo^D=ixCo^D-kr(^D,idims)\
      jxCo^D=ixCo^D+kr(^D,idims)\

      do iw=1,nwflux+nwaux
         slopeL=pwCo%w(ixCo^D,iw)-pwCo%w(hxCo^D,iw)
         slopeR=pwCo%w(jxCo^D,iw)-pwCo%w(ixCo^D,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idims)=signR*max(zero,min(dabs(slopeR), &
                                             signR*slopeL))
         case('woodward')
           slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), &
                              signR*slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idims)=signR*max(zero,min(mcbeta*dabs(slopeR), &
                              mcbeta*signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idims)=signR*max(zero,min(two*signR*slopeL, &
            (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idims)=signC*max(zero,min(dabs(slopeC), &
                             signC*slopeL,signC*slopeR))
         end select
      end do
   end do

   ! Interpolate from coarse cell using limited slopes
   pwFi%w(ixFi^D,1:nwflux+nwaux)=pwCo%w(ixCo^D,1:nwflux+nwaux)+{(slope(1:nwflux+nwaux,^D)*eta^D)+}

{end do\}

if(.not.slab)mygeo=>pgeoFi
{#IFNDEF DY_SP
if(covariant)myM => mygeo%m
}
! code test erase first
{#IFNDEF DY_SP
if (amrentropy) then
   call rhos_to_e(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x,patchfalse)
end if
}

end associate
end subroutine interpolation_linear
!=============================================================================
subroutine interpolation_copy(pwFi,ixFi^L,dxFi^D,xFimin^D, &
                              pwCo,dxCo^D,invdxCo^D,xComin^D)

integer, intent(in) :: ixFi^L
double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D
type(walloc) :: pwCo, pwFi

integer :: ixCo^D, ixFi^D
double precision :: xFi^D
!-----------------------------------------------------------------------------
associate(pgeoFi=>pstate(igrid)%geo)

{do ixFi^DB = ixFi^LIM^DB
   ! cell-centered coordinates of fine grid point
   xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB

   ! indices of coarse cell which contains the fine cell
   ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1\}

   ! Copy from coarse cell
   pwFi%w(ixFi^D,1:nwflux+nwaux)=pwCo%w(ixCo^D,1:nwflux+nwaux)

{end do\}



if(.not.slab)mygeo=>pgeoFi
{#IFNDEF DY_SP
if(covariant)myM => mygeo%m
}
! code test erase first
{#IFNDEF DY_SP
if (amrentropy) then
   call rhos_to_e(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x,patchfalse)
end if
}
end associate
end subroutine interpolation_copy
!=============================================================================
subroutine interpolation_unlimit(pwFi,ixFi^L,dxFi^D,xFimin^D, &
                                 pwCo,dxCo^D,invdxCo^D,xComin^D)

integer, intent(in) :: ixFi^L
double precision, intent(in) :: dxFi^D, xFimin^D, dxCo^D,invdxCo^D, xComin^D
type(walloc) :: pwCo, pwFi

integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, idims
double precision :: xCo^D, xFi^D, eta^D
double precision :: slope(nwflux+nwaux,ndim)
!-----------------------------------------------------------------------------
associate(pgeoFi=>pstate(igrid)%geo)

{do ixFi^DB = ixFi^LIM^DB
   ! cell-centered coordinates of fine grid point
   xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB

   ! indices of coarse cell which contains the fine cell
   ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1

   ! cell-centered coordinate for coarse cell
   xCo^DB=xComin^DB+(dble(ixCo^DB)-half)*dxCo^DB\}

   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta^D=(xFi^D-xCo^D)*invdxCo^D;
   else
      ix^D=2*int((ixFi^D+ixMlo^D)/2)-ixMlo^D;
      {eta^D=(xFi^D-xCo^D)*invdxCo^D &
            *two*(one-pgeoFi%dvolume(ixFi^DD) &
            /sum(pgeoFi%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
   end if

   do idims=1,ndim
      hxCo^D=ixCo^D-kr(^D,idims)\
      jxCo^D=ixCo^D+kr(^D,idims)\

      ! get centered slope
      slope(1:nwflux+nwaux,idims)=half*(pwCo%w(jxCo^D,1:nwflux+nwaux)-pwCo%w(hxCo^D,1:nwflux+nwaux))
   end do

   ! Interpolate from coarse cell using centered slopes
   pwFi%w(ixFi^D,1:nwflux+nwaux)=pwCo%w(ixCo^D,1:nwflux+nwaux)+{(slope(1:nwflux+nwaux,^D)*eta^D)+}
{end do\}



if(.not.slab)mygeo=>pgeoFi
{#IFNDEF DY_SP
if(covariant)myM => mygeo%m
}
! code test erase first
{#IFNDEF DY_SP
if (amrentropy) then
   call rhos_to_e(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x,patchfalse)
end if
}

end associate
end subroutine interpolation_unlimit
!=============================================================================


end subroutine gc_prolong
!===============================================================================
subroutine send_gc_r(igrid)
! Before this routine is called,
! the coarse representation of the block must be first filled
use mod_comm_gc
include 'amrvacdef.f'
! ibuf_send and isend_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i^D, ineighbor, ipe_neighbor, my_neighbor_type
{^IFPHI integer :: ipole}
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

{do i^DB=-1,1\}
   if (i^D==0|.and.) cycle
   my_neighbor_type=neighbor_type(i^D,igrid)
   if (my_neighbor_type==2) call bc_send_restrict
{end do\}

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_send_restrict
integer :: ic^D,inc^D,n_i^D,n_inc^D,ixS^L,idir,ibufaux,ibufnext
 ! Auxialiary array to avoid problems with the preprocessor...
 ! (It keeps splitting the line in the wrong place)
integer, dimension(1) :: sizes
!-----------------------------------------------------------------------------
ic^D=1+modulo(node(pig^D_,igrid)-1,2);
if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

ineighbor=neighbor(1,i^D,igrid)
ipe_neighbor=neighbor(2,i^D,igrid)
ipole=neighbor_pole(i^D,igrid)

if (ipe_neighbor.ne.mype) then

  inc^D=i^D+ic^D;
  ! The ghost region of the neighbor changes
  ! if there is a pole
  select case (ipole)
  case(0) ! No pole
    n_inc^D=-2*i^D+ic^D;
 {case (^D) ! Pole in this direction
    n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
  end select

  ! fill corresponding part of the send buffer...
  ixS^L=ixS_r_^L(i^D);
  ibufaux=ibuf_send_r
  ibufnext=ibufaux+sizes_r_send(i^D)

  sizes=(/sizes_r_send(i^D)/)

  sendbuffer_r(ibufaux:ibufnext-1)=reshape(pstateCo(igrid)%w%w(ixS^S,1:nwflux+nwaux),sizes)

  ibufaux=ibufnext
{#IFDEF STAGGERED
  do idir=1,^ND
    ixS^L=ixS_r_stg_^L(idir,i^D);
    ibufnext=ibufaux+sizes_r_send_stg(idir,i^D)
    sizes=(/sizes_r_send_stg(idir,i^D)/)
    sendbuffer_r(ibufaux:ibufnext-1)=&
      reshape(pstateCo(igrid)%ws%w(ixS^S,idir),sizes)   
    ibufaux=ibufnext
  end do
}

  ! ...and send
  itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
  isend_r=isend_r+1
  call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i^D), &
                 MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                 icomm,sendrequest_r(isend_r),ierrmpi)

  ibuf_send_r=ibufnext

end if

end subroutine bc_send_restrict
end subroutine send_gc_r
!===============================================================================
subroutine make_task_fill_gc_p(igrid)
use mod_comm_gc
include 'amrvacdef.f'

! In a loop, place MPI receive requests, updating the part of the receive
! buffer where ghost cells will be written
! Note: ibuf_recv_p and irecv_p were already initialized in init_comm_gc,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i^D, ic^D, inc^D, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------

call set_tmpGlobals(igrid)

{do i^DB=-1,1\}
   if (i^D==0|.and.) cycle
   my_neighbor_type=neighbor_type(i^D,igrid)
   if (my_neighbor_type.ne.2) cycle
   ic^D=1+modulo(node(pig^D_,igrid)-1,2);
   if ({(i^D==0.or.i^D==2*ic^D-3)|.and.}) then
     ipe_neighbor=neighbor(2,i^D,igrid)
     if (ipe_neighbor.ne.mype) then
     inc^D=ic^D+i^D;
     irecv_p=irecv_p+1
     itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
!     print *, 'recv tag',itag,'mype',mype,'igrid',igrid
       call MPI_IRECV(recvbuffer_p(ibuf_recv_p),sizes_p_recv_total(inc^D), &
                      MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                      icomm,recvrequest_p(irecv_p),ierrmpi)
       !call add_task_to_list() 
       ibuf_recv_p=ibuf_recv_p+sizes_p_recv_total(inc^D)
     end if
   end if
{end do\}

end subroutine make_task_fill_gc_p
!===============================================================================
subroutine fill_gc_p
use mod_comm_gc
include 'amrvacdef.f'

integer :: i^D, igrid, iigrid
integer :: ipole, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
! Wait for the receive buffer to be complete

if (nrecv_gc_p>0) then
   call MPI_WAITALL(nrecv_gc_p,recvrequest_p,recvstatus_p,ierrmpi)
   ibuf_recv_p=1
end if

! In a loop with the same order as that in make_task_fill_gc_p,
! directly fill the ghost cells that are in the same cpu and 
! unpack the receive buffer to fill those that are not.

do iigrid=1,igridstail; igrid=igrids(iigrid);
  call set_tmpGlobals(igrid)
     
 {do i^DB=-1,1\}
    if (i^D==0|.and.) cycle
    my_neighbor_type=neighbor_type(i^D,igrid)
    if (my_neighbor_type.eq.2) then
      call bc_fill_p
    end if
 {end do\}
end do


! Wait for the sends to complete and deallocate communication buffers
if (nrecv_gc_p>0) deallocate(recvbuffer_p,recvstatus_p,recvrequest_p)

if (nsend_gc_p>0) then
   call MPI_WAITALL(nsend_gc_p,sendrequest_p,sendstatus_p,ierrmpi)
   deallocate(sendbuffer_p,sendstatus_p,sendrequest_p)
end if


{^IFPHI
!! Deallocate pole buffers, now that all communications have ended
deallocate(pole_buf%w)
{#IFDEF STAGGERED
deallocate(pole_buf_stg%w)
}
}

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_fill_p
integer :: ic^D,inc^D,n_inc^D,ixS^L,ixR^L,idir,ibufnext
!--- for debugging ----
integer :: idx^D
!-----------------------------------------------------------------------------
ic^D=1+modulo(node(pig^D_,igrid)-1,2);
if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

ineighbor=neighbor(1,i^D,igrid)
ipe_neighbor=neighbor(2,i^D,igrid)
ipole=neighbor_pole(i^D,igrid)

if (ipole.eq.0) then   !! There is no pole 

inc^D=ic^D+i^D;
ixR^L=ixR_p_^L(inc^D);
if (ipe_neighbor.eq.mype) then !! Same processor
  n_inc^D=-2*i^D+ic^D;
  ixS^L=ixS_p_^L(n_inc^D);
  pstateCo(igrid)%w%w(ixR^S,1:nwflux+nwaux) &
          =pstate(ineighbor)%w%w(ixS^S,1:nwflux+nwaux)

  {#IFDEF STAGGERED
  do idir=1,^ND
     ixS^L=ixS_p_stg_^L(idir,n_inc^D);
     ixR^L=ixR_p_stg_^L(idir,inc^D);
     pstateCo(igrid)%ws%w(ixR^S,idir)=pstate(ineighbor)%ws%w(ixS^S,idir)
  end do
  }
else !! Different processor
  !! Unpack the buffer and fill the ghost cells
  ibufnext=ibuf_recv_p+sizes_p_recv(inc^D)
  pstateCo(igrid)%w%w(ixR^S,1:nwflux+nwaux)=reshape(source=recvbuffer_p(ibuf_recv_p:ibufnext-1),shape=shape(pstateCo(igrid)%w%w(ixR^S,1:nwflux+nwaux)))

  ibuf_recv_p=ibufnext

  {#IFDEF STAGGERED
  do idir=1,^ND
    ixR^L=ixR_p_stg_^L(idir,inc^D);
!      ixS^L=ixS_srl_stg_^L(idir,n_i^D);

    ibufnext=ibuf_recv_p+sizes_p_recv_stg(idir,inc^D)

    pstateCo(igrid)%ws%w(ixR^S,idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibufnext-1),shape=shape(pstateCo(igrid)%ws%w(ixR^S,idir)))
  
    ibuf_recv_p=ibufnext
  end do
  }
end if

else !! There is a pole
inc^D=ic^D+i^D;
!n_inc^D=inc^D; !! (Hope this is correct...)
select case (ipole)
{case (^D)
   n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
end select

!print *, 'n_inc',n_inc^D,'<--inc',inc^D


if (ipe_neighbor.eq.mype) then
  
  ixS^L=ixS_p_^L(n_inc^D);
  ixR^L=ixR_p_^L(inc^D);
  
  call pole_copy(pstateCo(igrid)%w,ixR^L,pstate(ineighbor)%w,ixS^L,ipole,i^D)
  
  {#IFDEF STAGGERED
  do idir=1,^ND
     ixS^L=ixS_p_stg_^L(idir,n_inc^D);
     ixR^L=ixR_p_stg_^L(idir,inc^D);
  
     call pole_copy_stg(pstateCo(igrid)%ws,ixR^L,pstate(ineighbor)%ws,ixS^L,ipole,i^D,idir)
  
  end do
  }

else
  ixR^L=ixR_p_^L(inc^D);
  !! Unpack the buffer and fill an auxiliary array
  ibufnext=ibuf_recv_p+sizes_p_recv(inc^D)
  pole_buf%w=zero
  pole_buf%w(ixR^S,1:nwflux+nwaux)=&
    reshape(source=recvbuffer_p(ibuf_recv_p:ibufnext-1),&
           shape=shape(pstateCo(igrid)%w%w(ixR^S,1:nwflux+nwaux)))
  ibuf_recv_p=ibufnext
 
  !! Fill ghost cells
  call pole_copy(pstateCo(igrid)%w,ixR^L,pole_buf,ixR^L,ipole,i^D)

  {#IFDEF STAGGERED
  pole_buf_stg%w=zero
  do idir=1,^ND
     ixR^L=ixR_p_stg_^L(idir,inc^D);

     ibufnext=ibuf_recv_p+sizes_p_recv_stg(idir,inc^D)
     pole_buf_stg%w(ixR^S,idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibufnext-1),&
       shape=shape(pstateCo(igrid)%ws%w(ixR^S,idir)))
     ibuf_recv_p=ibufnext

     call pole_copy_stg(pstateCo(igrid)%ws,ixR^L,pole_buf_stg,ixR^L,ipole,i^D,idir)

  end do
  }

end if

end if !! ipole == 0

end subroutine bc_fill_p
end subroutine fill_gc_p
!===============================================================================
subroutine send_gc_p(igrid)
use mod_comm_gc
include 'amrvacdef.f'
! ibuf_send and isend_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer :: igrid, i^D, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

{do i^DB=-1,1\}
   if (i^D==0|.and.) cycle
   my_neighbor_type=neighbor_type(i^D,igrid)
   if (my_neighbor_type==4) call bc_send_prolong
{end do\}



contains
!=============================================================================
subroutine bc_send_prolong
integer :: ic^D,inc^D,n_i^D,n_inc^D,ii^D,idir,ibufaux,ibufnext,ipole
integer :: ixS^L
 ! Auxialiary array to avoid problems with the preprocessor...
 ! (It keeps splitting the line in the wrong place)
integer, dimension(1) :: sizes
!-----------------------------------------------------------------------------

ipole=neighbor_pole(i^D,igrid)


{do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
  inc^DB=2*i^DB+ic^DB\}

  ineighbor=neighbor_child(1,inc^D,igrid)
  ipe_neighbor=neighbor_child(2,inc^D,igrid)

  if (ipe_neighbor.ne.mype) then
    select case (ipole)
    case(0) ! No pole
       n_i^D=-i^D;
       n_inc^D=ic^D+n_i^D;
   {case (^D) ! Pole in this direction
       n_inc^D=inc^D^D%n_inc^DD=ic^DD-i^DD;
!       print *, 'inc',inc^DD,'-->n_inc',n_inc^DD
       \}
    end select
  
   ! fill corresponding part of the send buffer...
  
    ixS^L=ixS_p_^L(inc^D);
  
    ibufaux=ibuf_send_p
    ibufnext=ibufaux+sizes_p_send(inc^D)
  
    sizes=(/sizes_p_send(inc^D)/)
  
    sendbuffer_p(ibufaux:ibufnext-1)=reshape(pstate(igrid)%w%w(ixS^S,1:nwflux+nwaux),sizes)
  
    ibufaux=ibufnext
  {#IFDEF STAGGERED
    do idir=1,^ND
      ixS^L=ixS_p_stg_^L(idir,inc^D);
      ibufnext=ibufaux+sizes_p_send_stg(idir,inc^D)
      sizes=(/sizes_p_send_stg(idir,inc^D)/)
      sendbuffer_p(ibufaux:ibufnext-1)=&
        reshape(pstate(igrid)%ws%w(ixS^S,idir),sizes)   
      ibufaux=ibufnext
    end do
  }

    ! ...and send
    itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}

    isend_p=isend_p+1
    call MPI_ISEND(sendbuffer_p(ibuf_send_p),sizes_p_send_total(inc^D), &
                   MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                   icomm,sendrequest_p(isend_p),ierrmpi)
  
    ibuf_send_p=ibufnext

  end if

{end do\}

end subroutine bc_send_prolong
end subroutine send_gc_p
!===============================================================================
subroutine fill_boundary(igrid)
! Physical boundary conditions
use mod_comm_gc
include 'amrvacdef.f'
integer, intent(in) :: igrid
! ... local ...
integer :: idims,iside,i^D,k^L,ixB^L
logical :: isphysbound
!-----------------------------------------------------------------------------

   call set_tmpGlobals(igrid)

   do idims=1,ndim
      ! to avoid using as yet unknown corner info in more than 1D, we
      ! fill only interior mesh ranges of the ghost cell ranges at first,
      ! and progressively enlarge the ranges to include corners later
      kmin1=0; kmax1=0;
      {^IFTWOD
       kmin2=merge(1, 0,  idims .lt. 2 .and. neighbor_type(0,-1,igrid)==1)
       kmax2=merge(1, 0,  idims .lt. 2 .and. neighbor_type(0, 1,igrid)==1)}
      {^IFTHREED
       kmin2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0,-1,0,igrid)==1)
       kmax2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0, 1,0,igrid)==1)
       kmin3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0,-1,igrid)==1)
       kmax3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0, 1,igrid)==1)}
      ixBmin^D=ixGlo^D+kmin^D*dixB;
      ixBmax^D=ixGhi^D-kmax^D*dixB;
      
      do iside=1,2
         i^D=kr(^D,idims)*(2*iside-3);
         if (neighbor_type(i^D,igrid)/=1) cycle
         
         call bc_phys(iside,idims,time,pstate(igrid),ixB^L)
         
      end do
         
   end do

end subroutine fill_boundary
!=============================================================================
subroutine physbound(i^D,igrid,isphysbound)
use mod_forest
include 'amrvacdef.f'

integer, intent(in)  :: i^D, igrid
logical, intent(out) :: isphysbound
type(tree_node_ptr)  :: tree
integer              :: level, ig^D, ign^D
!-----------------------------------------------------------------------------
isphysbound = .false.

tree%node => igrid_to_node(igrid,mype)%node
level = tree%node%level
{ig^D = tree%node%ig^D; }

{ign^D = ig^D + i^D; }
if ({ign^D .gt. ng^D(level) .or. ign^D .lt. 1|.or.}) isphysbound = .true.

end subroutine physbound
!=============================================================================
subroutine pole_copy(pwrecv,ixR^L,pwsend,ixS^L,ipole,i^D)
include 'amrvacdef.f'

integer, intent(in) :: ixR^L, ixS^L, ipole, i^D
type(walloc) :: pwrecv, pwsend

integer :: iw, iB, iside
!-----------------------------------------------------------------------------
select case (ipole)
{case (^D)
   iside=int((i^D+3)/2)
   iB=2*(^D-1)+iside
   do iw=1,nwflux+nwaux
      select case (typeB(iw,iB))
      case ("symm","polefix")
         pwrecv%w(ixR^S,iw) = pwsend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
      case ("asymm")
         pwrecv%w(ixR^S,iw) =-pwsend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
      case default
         call mpistop("Boundary condition at pole should be symm,asymm or polefix")
      end select
   end do \}
end select

end subroutine pole_copy
!=============================================================================
subroutine pole_copy_stg(pwrecv,ixR^L,pwsend,ixS^L,ipole,i^D,idir)
include 'amrvacdef.f'

integer, intent(in) :: ixR^L, ixS^L, ipole, i^D, idir
type(walloc) :: pwrecv, pwsend

integer :: iw, iB, iside
!-----------------------------------------------------------------------------

iw=idir+b0_

select case (ipole)
{case (^D)
   iside=int((i^D+3)/2)
   iB=2*(^D-1)+iside
   select case (typeB(iw,iB))
      case ("symm","polefix")
         pwrecv%w(ixR^S,idir) = pwsend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,idir)
      case ("asymm")
         pwrecv%w(ixR^S,idir) =-pwsend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,idir)
      case default
         call mpistop("Boundary condition at pole should be symm,asymm or polefix")
   end select
\}
end select
end subroutine pole_copy_stg
!=============================================================================
subroutine fix_auxiliary(igrid)
use mod_comm_gc
include 'amrvacdef.f'

integer, intent(in) :: igrid
! ... local ...
! Counters and auxiliaries
integer :: ix^L,i^D
!-----------------------------------------------------------------------------

associate(x=>pstate(igrid)%x%x,pgeoFi=>pstate(igrid)%geo)

call set_tmpGlobals(igrid)

saveigrid=igrid

      
{do i^DB=-1,1\}
! same-level auxiliaries have been sent and are not altered:
   if ((i^D==0|.and.).or.neighbor_type(i^D,igrid).eq.3) cycle

   ix^L=ixR_srl_^L(i^D);
   if(.not.slab)mygeo=>pgeoFi
   {#IFNDEF DY_SP
   if(covariant)myM => mygeo%m
   }
!code test
   call getaux(.true.,pstate(igrid)%w%w,x,ixG^LL,ix^L,"bc")
{end do\}

end associate

end subroutine fix_auxiliary
!=============================================================================


