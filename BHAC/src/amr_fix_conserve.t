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
subroutine init_comm_fix_conserve(idim^LIM)
use mod_fix_conserve
include 'amrvacdef.f'

integer, intent(in) :: idim^LIM

integer :: iigrid, igrid, idims, iside, i^D, nxCo^D
integer :: ic^D, inc^D, ipe_neighbor
integer :: recvsize, sendsize
{#IFDEF STAGGERED
integer :: recvsize_cc, sendsize_cc
integer :: pi^D,mi^D,ph^D,mh^D,idir}
!-----------------------------------------------------------------------------
nsend=0
nrecv=0
recvsize=0
sendsize=0
{#IFDEF STAGGERED
! Special communication for diagonal 'coarse corners'
! nrecv/send_cc (for 'coarse corners' is a dim=ndim-1 array which
! stores the faces that must be communicated in each direction.
! nrecv/send_ct (for 'corners total' is the total number of
! necessary communications. These special cases have their own
! send and receive buffers (send/recvbuffer_cc), their tags, etc.
nsend_ct=0
nrecv_ct=0
recvsize_cc=0
sendsize_cc=0
}
do idims= idim^LIM
   select case (idims)
   {case (^D)
      nrecv=nrecv+nrecv_fc(^D)
      nsend=nsend+nsend_fc(^D)
      nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-dixB;
      isize(^D)={nxCo^DD*}*(nwflux)
      recvsize=recvsize+nrecv_fc(^D)*isize(^D)
      sendsize=sendsize+nsend_fc(^D)*isize(^D)
{#IFDEF STAGGERED
      ! This does not consider the 'coarse corner' case
      nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-dixB+1;
      isize_stg(^D)={nxCo^DD*}*(^ND-1)
      ! To avoid writing IFDEF STAGGERED everywhere,
      ! the whole size is used (cell centred and staggered)
      isize(^D)=isize(^D)+isize_stg(^D)      
      recvsize=recvsize+nrecv_fc(^D)*isize_stg(^D)
      sendsize=sendsize+nsend_fc(^D)*isize_stg(^D)
      ! Coarse corner case
      nrecv_ct=nrecv_ct+nrecv_cc(^D)
      nsend_ct=nsend_ct+nsend_cc(^D)
      recvsize_cc=recvsize_cc+nrecv_cc(^D)*isize_stg(^D)
      sendsize_cc=sendsize_cc+nsend_cc(^D)*isize_stg(^D)

}
   \}
   end select
end do

if (nrecv>0) then
   ! Receive for direct neighbors
   allocate(recvbuffer(recvsize),recvstatus(MPI_STATUS_SIZE,nrecv), &
            recvrequest(nrecv))
   recvrequest=MPI_REQUEST_NULL

   ibuf=1
   irecv=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_reflux_direct(igrid,idim^LIM)
   end do

end if

{#IFDEF STAGGERED
if (nrecv_ct.gt.0) then
   ! Receive corners
   allocate(recvbuffer_cc(recvsize_cc),recvstatus_stg(MPI_STATUS_SIZE,nrecv_ct), &
            recvrequest_stg(nrecv_ct))
   recvrequest_stg=MPI_REQUEST_NULL

   ibuf_cc=1
   irecv_cc=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_reflux_corners(igrid,idim^LIM)
   end do

end if
}
if (nsend>0) then
   allocate(sendbuffer(sendsize),sendstatus(MPI_STATUS_SIZE,nsend),sendrequest(nsend))
   sendrequest=MPI_REQUEST_NULL
   isend=0
   ibuf_send=1
end if
{#IFDEF STAGGERED
if (nsend_ct>0) then
   allocate(sendbuffer_cc(sendsize_cc),sendstatus_stg(MPI_STATUS_SIZE,nsend_ct),sendrequest_stg(nsend_ct))
   sendrequest_stg=MPI_REQUEST_NULL
   isend_cc=0
   ibuf_cc_send=1
end if
}

end subroutine init_comm_fix_conserve
!=============================================================================
subroutine make_task_reflux_direct(igrid,idim^LIM)
use mod_fix_conserve
include 'amrvacdef.f'

integer, intent(in) :: igrid,idim^LIM
integer             :: iside,idims,i^D,ic^D,inc^D
integer             :: ipe_neighbor
!-----------------------------------------------------------------------------

do idims= idim^LIM
   do iside=1,2
      i^D=kr(^D,idims)*(2*iside-3);
      {^IFPHI if (neighbor_pole(i^D,igrid)/=0) cycle}
      if (neighbor_type(i^D,igrid)/=4) cycle
      {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
         inc^DB=2*i^DB+ic^DB\}
         ipe_neighbor=neighbor_child(2,inc^D,igrid)
         if (ipe_neighbor/=mype) then
            irecv=irecv+1
            itag=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
            call MPI_IRECV(recvbuffer(ibuf),isize(idims), &
                           MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                           icomm,recvrequest(irecv),ierrmpi)
           !call add_task_to_list(itag,recvrequest(irecv),recvstatus(irecv)) 
            ibuf=ibuf+isize(idims)
         end if
      {end do\}
   end do
end do

end subroutine make_task_reflux_direct
!=============================================================================
{#IFDEF STAGGERED
subroutine make_task_reflux_corners(igrid,idim^LIM)
use mod_fix_conserve
include 'amrvacdef.f'

integer, intent(in) :: igrid,idim^LIM
integer :: iside,idims,i^D,ic^D,inc^D
integer :: ipe_neighbor
integer :: pi^D,mi^D,ph^D,mh^D,idir
!-----------------------------------------------------------------------------

do idims= idim^LIM
   do iside=1,2
     i^D=kr(^D,idims)*(2*iside-3);
     ! Check if there are special corners
     ! (Coarse block diagonal to a fine block)
     ! If there are, receive.
     ! Tags are calculated in the same way as for
     ! normal fluxes, but should not overlap because
     ! inc^D are different
     if (neighbor_type(i^D,igrid).eq.3) then
       do idir=idims+1,ndim
         pi^D=i^D+kr(idir,^D);
         mi^D=i^D-kr(idir,^D);
         ph^D=pi^D-kr(idims,^D)*(2*iside-3);
         mh^D=mi^D-kr(idims,^D)*(2*iside-3);

{^IFPHI  if (neighbor_pole(pi^D,igrid).eq.0) then}
         if (neighbor_type(pi^D,igrid).eq.4.and.&
             neighbor_type(ph^D,igrid).eq.3) then
          ! Loop on children (several in 3D)
          {do ic^DB=1+int((1-pi^DB)/2),2-int((1+pi^DB)/2)
              inc^DB=2*pi^DB+ic^DB\}
              ipe_neighbor=neighbor_child(2,inc^D,igrid)
              if (mype.ne.ipe_neighbor) then
                irecv_cc=irecv_cc+1
                itag_cc=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
                call MPI_IRECV(recvbuffer_cc(ibuf_cc),isize_stg(idims),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,recvrequest_stg(irecv_cc),ierrmpi)
               !call add_task_to_list(itag_cc,recvrequest_stg(irecv_cc),recvstatus_stg(irecv_cc))
                ibuf_cc=ibuf_cc+isize_stg(idims)
              end if
          {end do\}
         end if
{^IFPHI  end if}

{^IFPHI  if (neighbor_pole(mi^D,igrid).eq.0) then}
         if (neighbor_type(mi^D,igrid).eq.4.and.&
             neighbor_type(mh^D,igrid).eq.3) then
          ! Loop on children (several in 3D)
          {do ic^DB=1+int((1-mi^DB)/2),2-int((1+mi^DB)/2)
              inc^DB=2*mi^DB+ic^DB\}
              ipe_neighbor=neighbor_child(2,inc^D,igrid)
              if (mype.ne.ipe_neighbor) then
                irecv_cc=irecv_cc+1
                itag_cc=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
                call MPI_IRECV(recvbuffer_cc(ibuf_cc),isize_stg(idims),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,recvrequest_stg(irecv_cc),ierrmpi)
               !call add_task_to_list(itag_cc,recvrequest_stg(irecv_cc),recvstatus_stg(irecv_cc))
                ibuf_cc=ibuf_cc+isize_stg(idims)
              end if
          {end do\}
         end if
{^IFPHI  end if}
       end do
     end if
   end do
end do


end subroutine make_task_reflux_corners
!=============================================================================
}
subroutine allocateBflux
use mod_fix_conserve
include 'amrvacdef.f'

integer :: iigrid, igrid, iside, i^D, nx^D, nxCo^D
{#IFDEF STAGGERED
integer :: idir,idim,pi^D, mi^D, ph^D, mh^D ! To detect corners
}
!-----------------------------------------------------------------------------
nx^D=ixMhi^D-ixMlo^D+1;
nxCo^D=nx^D/2;

do iigrid=1,igridstail; igrid=igrids(iigrid);
   ! For every grid,
   ! arrays for the fluxes are allocated for every face direction(^D)
   ! and every side (1=left, 2=right)
   {do iside=1,2
      i^DD=kr(^DD,^D)*(2*iside-3);
      {^IFPHI if (neighbor_pole(i^DD,igrid)/=0) cycle}
      select case (neighbor_type(i^DD,igrid))
      case (4)
         allocate(pflux(iside,^D,igrid)%flux(1^D%1:nx^DD,1:nwflux){#IFDEF STAGGERED , pflux(iside,^D,igrid)%edge(1^D%0:nx^DD,1:^ND-1)})
      case (2)
         allocate(pflux(iside,^D,igrid)%flux(1^D%1:nxCo^DD,1:nwflux){#IFDEF STAGGERED , pflux(iside,^D,igrid)%edge(1^D%0:nxCo^DD,1:^ND-1)})
     {#IFDEF STAGGERED
      case(3)
      ! If there are staggered variables, it is necessary
      ! to store also some corners.
      ! Detect the corners that must be stored
      ! and store them
      ! This produces a warning when D=ndim, it might be improved.
      idim=^D
      do idir=idim+1,ndim
        pi^DD=i^DD+kr(idir,^DD);
        mi^DD=i^DD-kr(idir,^DD);
        ph^DD=pi^DD-kr(^D,^DD)*(2*iside-3);
        mh^DD=mi^DD-kr(^D,^DD)*(2*iside-3);
        if ((neighbor_type(pi^DD,igrid).eq.4&
            .and.neighbor_type(ph^DD,igrid).eq.3)&
            .or.&
            (neighbor_type(mi^DD,igrid).eq.4&
            .and.neighbor_type(mh^DD,igrid).eq.3)) then
          allocate(pflux(iside,^D,igrid)%edge(1^D%0:nx^DD,1:^ND-1))
          exit
        end if
      end do
     }
      end select
   end do\}
end do

end subroutine allocateBflux
!=============================================================================
subroutine deallocateBflux
use mod_fix_conserve
include 'amrvacdef.f'

integer :: iigrid, igrid, iside
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
   {do iside=1,2
      if (associated(pflux(iside,^D,igrid)%flux)) &
         deallocate(pflux(iside,^D,igrid)%flux)
!code test
         nullify(pflux(iside,^D,igrid)%flux)

{#IFDEF STAGGERED
      if (associated(pflux(iside,^D,igrid)%edge)) &
         deallocate(pflux(iside,^D,igrid)%edge)
!code test
         nullify(pflux(iside,^D,igrid)%edge)
}
   end do\}
end do

end subroutine deallocateBflux
!=============================================================================
subroutine fix_conserve(psuse,idim^LIM)
use mod_fix_conserve
include 'amrvacdef.f'

type(state) :: psuse(ngridshi)
integer, intent(in) :: idim^LIM

integer :: iigrid, igrid, idims, iside, iotherside, i^D, ic^D, inc^D, ix^L
integer :: nxCo^D, iw, ix, ipe_neighbor, ineighbor, ibufnext, nbuf
double precision :: CoFiratio
!-----------------------------------------------------------------------------
if (slab) then
   ! The flux is divided by volume of fine cell. We need, however,
   ! to divide by volume of coarse cell => muliply by volume ratio
   CoFiratio=one/dble(2**ndim)
end if

if (nrecv>0) then
   call MPI_WAITALL(nrecv,recvrequest,recvstatus,ierrmpi)
   ibuf=1
end if
nxCo^D=(ixMhi^D-ixMlo^D+1)/2;

!if (.false.) then 

! for all grids: perform flux update at Coarse-Fine interfaces
do iigrid=1,igridstail; igrid=igrids(iigrid);
   do idims= idim^LIM
      select case (idims)
      {case (^D)
         do iside=1,2
            i^DD=kr(^DD,^D)*(2*iside-3);
            {^IFPHI if (neighbor_pole(i^DD,igrid)/=0) cycle}
            if (neighbor_type(i^DD,igrid)/=4) cycle

! opedit: skip over active/passive interface since flux for passive ones is 
! not computed, keep the buffer counter up to date:
            if (.not.neighbor_active(i^DD,igrid).or.&
                .not.neighbor_active(0^DD&,igrid) ) then
               {do ic^DDB=1+int((1-i^DDB)/2),2-int((1+i^DDB)/2)
               inc^DDB=2*i^DDB+ic^DDB\}
               ipe_neighbor=neighbor_child(2,inc^DD,igrid)
               if (ipe_neighbor/=mype) then
                  ibufnext=ibuf+isize(^D)
                  ibuf=ibufnext
                  end if
               {end do\}
               cycle
            end if
!

            select case (iside)
            case (1)
               ix=ixMlo^D
            case (2)
               ix=ixMhi^D
            end select

            ! remove coarse flux
            if (slab) then
               psuse(igrid)%w%w(ix^D%ixM^T,1:nwflux) &
                  = psuse(igrid)%w%w(ix^D%ixM^T,1:nwflux) &
                   -pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux)
            else
               do iw=1,nwflux
                  psuse(igrid)%w%w(ix^D%ixM^T,iw)=psuse(igrid)%w%w(ix^D%ixM^T,iw)&
                     -pflux(iside,^D,igrid)%flux(1^D%:^DD&,iw) &
                     /pgeo(igrid)%dvolume(ix^D%ixM^T)
               end do
            end if


            ! add fine flux
            {do ic^DDB=1+int((1-i^DDB)/2),2-int((1+i^DDB)/2)
               inc^DDB=2*i^DDB+ic^DDB\}
               ineighbor=neighbor_child(1,inc^DD,igrid)
               ipe_neighbor=neighbor_child(2,inc^DD,igrid)
               ixmin^D=ix^D%ixmin^DD=ixMlo^DD+(ic^DD-1)*nxCo^DD;
               ixmax^D=ix^D%ixmax^DD=ixmin^DD-1+nxCo^DD;
               if (ipe_neighbor==mype) then
                  iotherside=3-iside
                  if (slab) then
                     psuse(igrid)%w%w(ix^S,1:nwflux) &
                       = psuse(igrid)%w%w(ix^S,1:nwflux) &
                       + pflux(iotherside,^D,ineighbor)%flux(:^DD&,1:nwflux)&
                       * CoFiratio
                  else
                     do iw=1,nwflux
                        psuse(igrid)%w%w(ix^S,iw)=psuse(igrid)%w%w(ix^S,iw) &
                            +pflux(iotherside,^D,ineighbor)%flux(:^DD&,iw) &
                            /pgeo(igrid)%dvolume(ix^S)
                     end do
                  end if
               else
                  if (slab) then
                     ibufnext=ibuf+isize(^D)
                     {#IFNDEF STAGGERED
                     psuse(igrid)%w%w(ix^S,1:nwflux) &
                         = psuse(igrid)%w%w(ix^S,1:nwflux)+CoFiratio*reshape(source=recvbuffer(ibuf:ibufnext-1), &
                                 shape=shape(psuse(igrid)%w%w(ix^S,1:nwflux)))
                     }
                     {#IFDEF STAGGERED
                     psuse(igrid)%w%w(ix^S,1:nwflux) &
                         = psuse(igrid)%w%w(ix^S,1:nwflux)+CoFiratio*reshape(source=recvbuffer(ibuf:ibufnext-isize_stg(^D)-1), &
                                 shape=shape(psuse(igrid)%w%w(ix^S,1:nwflux)))
                     }
                     ibuf=ibufnext
                  else
                     ibufnext=ibuf+isize(^D)
                     {#IFNDEF STAGGERED
                     nbuf=isize(^D)/nwflux
                     }
                     {#IFDEF STAGGERED
                     nbuf=(isize(^D)-isize_stg(^D))/nwflux
                     }
                     do iw=1,nwflux
                        psuse(igrid)%w%w(ix^S,iw)=psuse(igrid)%w%w(ix^S,iw) &
                           +reshape(source=recvbuffer(ibuf:ibuf+nbuf-1), &
                                    shape=shape(psuse(igrid)%w%w(ix^S,iw))) &
                           /pgeo(igrid)%dvolume(ix^S)
                        ibuf=ibuf+nbuf
                     end do
                     ibuf=ibufnext
                  end if
               end if
            {end do\}
         end do\}
      end select
   end do
end do

!end if

{#IFNDEF STAGGERED
if (nrecv>0) deallocate(recvbuffer,recvstatus,recvrequest)

if (nsend>0) then
   call MPI_WAITALL(nsend,sendrequest,sendstatus,ierrmpi)
   deallocate(sendbuffer,sendstatus,sendrequest)
end if
}

end subroutine fix_conserve
!=============================================================================
subroutine storeflux(igrid,fC,idim^LIM)
use mod_fix_conserve
include 'amrvacdef.f'

integer, intent(in)          :: igrid, idim^LIM
double precision, intent(in) :: fC(ixG^T,1:nwflux,1:ndim)

integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw
!-----------------------------------------------------------------------------
do idims = idim^LIM
   select case (idims)
   {case (^D)
      do iside=1,2
         i^DD=kr(^DD,^D)*(2*iside-3);
         {^IFPHI if (neighbor_pole(i^DD,igrid)/=0) cycle}
         select case (neighbor_type(i^DD,igrid))
         case (4)
            select case (iside)
            case (1)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux) = &
                  -fC(dixB^D%ixM^T,1:nwflux,^D)
            case (2)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux) = &
                  fC(ixMhi^D^D%ixM^T,1:nwflux,^D)
            end select
         case (2)
            nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-dixB;
            select case (iside)
            case (1)
               do iw=1,nwflux
                  {do ixCo^DDB=1,nxCo^DDB\}
                     ix^D=dixB^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
                     pflux(iside,^D,igrid)%flux(ixCo^DD,iw) &
                        = {^NOONEDsum}(fC(ix^D^D%ix^DD:ix^DD+1,iw,^D))
                  {end do\}
               end do
            case (2)
               do iw=1,nwflux
                  {do ixCo^DDB=1,nxCo^DDB\}
                     ix^D=ixMhi^D^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
                     pflux(iside,^D,igrid)%flux(ixCo^DD,iw) &
                        =-{^NOONEDsum}(fC(ix^D^D%ix^DD:ix^DD+1,iw,^D))
                  {end do\}
               end do
            end select
         end select
      end do\}
   end select
end do

end subroutine storeflux
!=============================================================================
subroutine sendflux(igrid,idim^LIM)
use mod_fix_conserve
include 'amrvacdef.f'

integer, intent(in) :: idim^LIM
integer             :: igrid, iigrid
integer             :: ibuf_send_next
integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw
integer :: ineighbor, ipe_neighbor
{#IFDEF STAGGERED
integer             :: idir,ibuf_cc_send_next,pi^D,ph^D,mi^D,mh^D
}
!----------------------------------------------------------------------------

do idims = idim^LIM
   select case (idims)
   {case (^D)
      do iside=1,2
         i^DD=kr(^DD,^D)*(2*iside-3);
         {^IFPHI if (neighbor_pole(i^DD,igrid)/=0) cycle}
         
         if (neighbor_type(i^DD,igrid)==2) then

            ineighbor=neighbor(1,i^DD,igrid)
            ipe_neighbor=neighbor(2,i^DD,igrid)
            if (ipe_neighbor/=mype) then
               ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
               inc^D=-2*i^D+ic^D^D%inc^DD=ic^DD;
               itag=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
               isend=isend+1

               ibuf_send_next=ibuf_send+isize(^D)
               {#IFNDEF STAGGERED
               sendbuffer(ibuf_send:ibuf_send_next-1)=&
               reshape(pflux(iside,^D,igrid)%flux,(/isize(^D)/))}

               {#IFDEF STAGGERED
               sendbuffer(ibuf_send:ibuf_send_next-isize_stg(^D)-1)=&
               reshape(pflux(iside,^D,igrid)%flux,(/isize(^D)-isize_stg(^D)/))

               sendbuffer(ibuf_send_next-isize_stg(^D):ibuf_send_next-1)=&
               reshape(pflux(iside,^D,igrid)%edge,(/isize_stg(^D)/))}

               call MPI_ISEND(sendbuffer(ibuf_send),isize(^D), &
                              MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                              icomm,sendrequest(isend),ierrmpi)

               ibuf_send=ibuf_send_next

            end if
            {#IFDEF STAGGERED
            ! If we are in a fine block surrounded by coarse blocks
            do idir=idims+1,ndim
              pi^DD=i^DD+kr(idir,^DD);
              mi^DD=i^DD-kr(idir,^DD);
              ph^DD=pi^DD-kr(idims,^DD)*(2*iside-3);
              mh^DD=mi^DD-kr(idims,^DD)*(2*iside-3);

              if (neighbor_type(pi^DD,igrid).eq.2.and.&
                  neighbor_type(ph^DD,igrid).eq.2.and.&
                  mype.ne.neighbor(2,pi^DD,igrid)) then
{^IFPHI         if (neighbor_pole(pi^DD,igrid).eq.0) then }
                ! Get relative position in the grid for tags
                ineighbor=neighbor(1,pi^DD,igrid)
                ipe_neighbor=neighbor(2,pi^DD,igrid)
                ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
                inc^DD=-2*pi^DD+ic^DD;
                itag_cc=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
                ! Reshape to buffer and send
                isend_cc=isend_cc+1
                ibuf_cc_send_next=ibuf_cc_send+isize_stg(^D)
                sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                reshape(pflux(iside,^D,igrid)%edge,shape=(/isize_stg(^D)/))
                call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(^D),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,sendrequest_stg(isend_cc),ierrmpi)
                ibuf_cc_send=ibuf_cc_send_next
{^IFPHI         end if }
              end if
 
              if (neighbor_type(mi^DD,igrid).eq.2.and.&
                  neighbor_type(mh^DD,igrid).eq.2.and.&
                  mype.ne.neighbor(2,mi^DD,igrid)) then
{^IFPHI         if (neighbor_pole(mi^DD,igrid).eq.0) then }
                ! Get relative position in the grid for tags
                ineighbor=neighbor(1,mi^DD,igrid)
                ipe_neighbor=neighbor(2,mi^DD,igrid)
                ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
                inc^DD=-2*pi^DD+ic^DD;
                inc^DD=-2*mi^DD+ic^DD;
                itag_cc=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
                ! Reshape to buffer and send
                isend_cc=isend_cc+1
                ibuf_cc_send_next=ibuf_cc_send+isize_stg(^D)
                sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                
                reshape(pflux(iside,^D,igrid)%edge,shape=(/isize_stg(^D)/))
                call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(^D),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,sendrequest_stg(isend_cc),ierrmpi)
                ibuf_cc_send=ibuf_cc_send_next
{^IFPHI         end if }
              end if
            end do
            }
         end if
      end do\}
   end select
end do

end subroutine sendflux
!=============================================================================
