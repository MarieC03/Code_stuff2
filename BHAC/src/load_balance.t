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
subroutine load_balance
use mod_forest
include 'amrvacdef.f'

integer :: Morton_no, recv_igrid, recv_ipe, send_igrid, send_ipe, igrid, ipe

integer, external :: getnode
!-----------------------------------------------------------------------------
if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
call get_Morton_range_active

if (npe==1) then
   sfc_to_igrid(:)=sfc(1,Morton_start(mype):Morton_stop(mype))
   return
end if

irecv=0
isend=0

do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
   recv_ipe=ipe

   send_igrid=sfc(1,Morton_no)
   send_ipe=sfc(2,Morton_no)

   if (recv_ipe/=send_ipe) then
      recv_igrid=getnode(recv_ipe)
            call change_ipe_tree_leaf(recv_igrid,recv_ipe,send_igrid,send_ipe)
      if (recv_ipe==mype) call lb_recv
      if (send_ipe==mype) call lb_send
   end if
   if (recv_ipe==mype) then
      if (recv_ipe==send_ipe) then
         sfc_to_igrid(Morton_no)=send_igrid
      else
         sfc_to_igrid(Morton_no)=recv_igrid
      end if
   end if
end do; end do

if (irecv>0) call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
{#IFDEF STAGGERED
if (irecv>0) call MPI_WAITALL(irecv,recvrequest_stg,recvstatus_stg,ierrmpi)
}

if (isend>0) call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
{#IFDEF STAGGERED
if (isend>0) call MPI_WAITALL(isend,sendrequest_stg,sendstatus_stg,ierrmpi)
}

! post processing
do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
   recv_ipe=ipe

   send_igrid=sfc(1,Morton_no)
   send_ipe=sfc(2,Morton_no)

   if (recv_ipe/=send_ipe) then
      if (send_ipe==mype) call dealloc_node(send_igrid)
      call putnode(send_igrid,send_ipe)
   end if
end do; end do

call amr_Morton_order

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine lb_recv
!-----------------------------------------------------------------------------
call alloc_node(recv_igrid)

itag=recv_igrid
irecv=irecv+1
call MPI_IRECV(pw(recv_igrid)%w,1,type_block_io,send_ipe,itag, &
     icomm,recvrequest(irecv),ierrmpi)

{#IFDEF STAGGERED
itag=recv_igrid+ngridshi
call MPI_IRECV(pws(recv_igrid)%w,1,type_block_staggered_io,send_ipe,itag, &
     icomm,recvrequest_stg(irecv),ierrmpi)
}

end subroutine lb_recv
!=============================================================================
subroutine lb_send
!-----------------------------------------------------------------------------
itag=recv_igrid
isend=isend+1
call MPI_ISEND(pw(send_igrid)%w,1,type_block_io,recv_ipe,itag, &
               icomm,sendrequest(isend),ierrmpi)

{#IFDEF STAGGERED
itag=recv_igrid+ngridshi
call MPI_ISEND(pws(send_igrid)%w,1,type_block_staggered_io,recv_ipe,itag, &
               icomm,sendrequest_stg(isend),ierrmpi)
}

end subroutine lb_send
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine load_balance
!=============================================================================
subroutine amr_Morton_order
! Construct Morton-order as a global recursive lexicographic ordering.

use mod_forest
include 'amrvacdef.f'

integer :: ig^D, Morton_no, isfc
!-----------------------------------------------------------------------------
Morton_no=0

!code test
nglev1={ng^D(1)*}
do isfc=1,nglev1
   ig^D=sfc_iglevel1(^D,isfc)\
   call get_Morton_number(tree_root(ig^D))
end do



!{#IFDEF SFC_ON_BASE
!do isfc=1,nglev1
!   ig^D=sfc_iglevel1(^D,isfc)\ 
!   call get_Morton_number(tree_root(ig^D))
!end do
!}{#IFNDEF SFC_ON_BASE
!{do ig^DB=1,ng^DB(1)\}
!   call get_Morton_number(tree_root(ig^D))
!{end do\}
!
!}

if (Morton_no/=nleafs) then
   call mpistop("error in amr_Morton_order: Morton_no/=nleafs")
end if

contains
!=============================================================================
! internal procedures
!=============================================================================
recursive subroutine get_Morton_number(tree)

type(tree_node_ptr) :: tree

integer :: ic^D
!-----------------------------------------------------------------------------
if (tree%node%leaf) then
   Morton_no=Morton_no+1
   sfc(1,Morton_no)=tree%node%igrid
   sfc(2,Morton_no)=tree%node%ipe
   if (tree%node%active) then 
      sfc(3,Morton_no)=1 
   else 
      sfc(3,Morton_no)=0 
   end if
   igrid_to_sfc(tree%node%igrid)=Morton_no
else
   {do ic^DB=1,2\}
      call get_Morton_number(tree%node%child(ic^D))
   {end do\}
end if

end subroutine get_Morton_number
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine amr_Morton_order
!=============================================================================
subroutine get_Morton_range
use mod_forest
include 'amrvacdef.f'

integer :: ipe
!-----------------------------------------------------------------------------
if (allocated(sfc_to_igrid)) deallocate(sfc_to_igrid)
do ipe=0,npe-1
   Morton_start(ipe)=1+ipe*int(nleafs/npe)+min(ipe,mod(nleafs,npe))
   Morton_stop(ipe)=(ipe+1)*int(nleafs/npe)+min(ipe+1,mod(nleafs,npe))
end do
if (Morton_stop(mype)>=Morton_start(mype)) then
   allocate(sfc_to_igrid(Morton_start(mype):Morton_stop(mype)))
end if

end subroutine get_Morton_range
!=============================================================================

!=============================================================================
subroutine get_Morton_range_active
use mod_forest
include 'amrvacdef.f'

! Cut the sfc based on weighted decision.  
! Oliver Porth, 02.02.2012

!!!Here you choose the weithts:!!
integer, parameter :: wa=3, wp=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! wp : Weight for passive block
! wa : Weight for active block
! wp=0 : balance load (active blocks) exactly and 
! don't care about memory imbalance
! wp=wa : balance memory exactly and don't care about load
! Maximum possible memory imbalance is X=wa/wp.

! If you run into memory issues, decrease this ratio.
! Scaling should be better if you allow for higher ratio.  
! Note also that passive cells still do regridding and boundary swap.  
! I have best results with a ratio 2:1, but it is problem dependent.
! It can't do magic though...
! Best to make sure that the sfc is properly aligned with your problem. 

integer :: ipe, Morton_no
integer :: Mtot, Mstop, Mcurr
! For debugging: 
integer :: nactive(0:npe-1),npassive(0:npe-1)
!-----------------------------------------------------------------------------
if (allocated(sfc_to_igrid)) deallocate(sfc_to_igrid)

Mtot  = nleafs_active*wa+(nleafs-nleafs_active)*wp
ipe = 0 
Mcurr = 0

nactive=0
npassive=0

Morton_start(0) = 1
do Morton_no=1,nleafs
   ! This is where we ideally would like to make the cuts:
   Mstop  = (ipe+1)*int(Mtot/npe)+min(ipe+1,mod(Mtot,npe))
   ! Build up mass:
   Mcurr = Mcurr + (wa*sfc(3,Morton_no)+wp*(1-sfc(3,Morton_no)))

   if (sfc(3,Morton_no)==1) then 
      nactive(ipe) = nactive(ipe) +1
      else
         npassive(ipe) = npassive(ipe) +1
      end if

   if (Mcurr >= Mstop) then 
      Morton_stop(ipe) = Morton_no
      ipe = ipe +1
      if (ipe>=npe) exit
      Morton_start(ipe) = Morton_no + 1
   end if
end do

Xmemory=dble(maxval(npassive+nactive))/&
     dble(minval(npassive+nactive))
Xload=dble(maxval(nactive))/&
     dble(minval(nactive))

if (Morton_stop(mype)>=Morton_start(mype)) then
   allocate(sfc_to_igrid(Morton_start(mype):Morton_stop(mype)))
end if

end subroutine get_Morton_range_active
!=============================================================================
