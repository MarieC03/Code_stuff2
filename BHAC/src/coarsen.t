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
subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)

include 'amrvacdef.f'

integer, intent(in) :: igrid, ipe
integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
logical, intent(in) :: active

integer :: igridFi, ipeFi, ixCo^L, ixCoG^L, ixCoM^L, ic^D
{#IFDEF STAGGERED
integer :: idir
}
!-----------------------------------------------------------------------------
if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)

if (ipe==mype) call alloc_node(igrid)  !get sqrtgamma dV dS

! New passive cell, coarsen from initial condition:
if (.not. active) then

   if (ipe == mype) then
      call initial_condition(igrid)
      {do ic^DB=1,2\}
      igridFi=child_igrid(ic^D)
      ipeFi=child_ipe(ic^D)
      if (ipeFi==mype) then
         ! remove solution space of child      
         call dealloc_node(igridFi)
      end if
      {end do\}
      
   end if

   return
end if


{do ic^DB=1,2\}
   igridFi=child_igrid(ic^D)
   ipeFi=child_ipe(ic^D)

   if (ipeFi==mype) then
      ^D&dxlevel(^D)=rnode(rpdx^D_,igridFi);
      if (ipe==mype) then
         ixComin^D=ixMlo^D+(ic^D-1)*(ixMhi^D-ixMlo^D+1)/2;
         ixComax^D=ixMhi^D+(ic^D-2)*(ixMhi^D-ixMlo^D+1)/2;

         call coarsen_grid(ps(igridFi),px(igridFi)%x,ixG^LL,ixM^LL,ps(igrid),px(igrid)%x,ixG^LL, &
                     ixCo^L,pgeo(igridFi),pgeo(igrid),restrictprimitive,.false.)

         ! remove solution space of child
         call dealloc_node(igridFi)
      else
         ixCoGmin^D=1;
        ixCoGmax^D=ixGhi^D/2+dixB;
         ixCoM^L=ixCoG^L^LSUBdixB;
         call coarsen_grid(ps(igridFi),px(igridFi)%x,ixG^LL,ixM^LL,psCoarse(igridFi),pxCoarse(igridFi)%x, &
                           ixCoG^L,ixCoM^L,pgeo(igridFi),pgeoCoarse(igridFi),&
                           restrictprimitive,.false.)

         itag=ipeFi*ngridshi+igridFi
         isend=isend+1
         call MPI_ISEND(pwCoarse(igridFi)%w,1,type_coarse_block,ipe,itag, &
                        icomm,sendrequest(isend),ierrmpi)
{#IFDEF STAGGERED
         do idir=1,ndim
         itag_stg=(npe+ipeFi+1)*ngridshi+igridFi*(^NC-1+idir)
         call MPI_ISEND(psCoarse(igridFi)%ws%w,1,type_coarse_block_stg(idir,ic^D),ipe,itag_stg, &
                        icomm,sendrequest_stg(isend),ierrmpi)
         end do
}
         
      end if
   else
      if (ipe==mype) then
         itag=ipeFi*ngridshi+igridFi
         irecv=irecv+1
         call MPI_IRECV(pw(igrid)%w,1,type_sub_block(ic^D),ipeFi,itag, &
                        icomm,recvrequest(irecv),ierrmpi)
{#IFDEF STAGGERED
         do idir=1,ndim
         itag_stg=(npe+ipeFi+1)*ngridshi+igridFi*(^NC-1+idir)
         call MPI_IRECV(ps(igrid)%ws%w,1,type_sub_block_stg(idir,ic^D),ipeFi,itag_stg, &
                        icomm,recvrequest_stg(irecv),ierrmpi)
         end do
}
         
      end if
   end if
{end do\}

if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
end subroutine coarsen_grid_siblings
!=============================================================================
subroutine coarsen_grid(sFi,xFi,ixFiG^L,ixFi^L,sCo,xCo,ixCoG^L,ixCo^L,&
                        pgeogrid,pgeoCoarsegrid,coarsenprim,keepFi)

! coarsen by 2 in every direction - conservatively
include 'amrvacdef.f'

integer, intent(in)             :: ixFiG^L, ixFi^L, ixCoG^L, ixCo^L
double precision, intent(inout) :: xFi(ixFiG^S,1:ndim)
double precision,intent(inout)  :: xCo(ixCoG^S,1:ndim)
type(state), intent(inout)      :: sFi, sCo
type(geoalloc)                  :: pgeogrid, pgeoCoarsegrid
logical, intent(in)             :: coarsenprim, keepFi

! .. local ..
integer :: ixCo^D, ixFi^D, iw
double precision :: CoFiratio
!-----------------------------------------------------------------------------
associate(wFi=>sFi%w%w(ixFiG^S,1:nw), wCo=>sCo%w%w(ixCoG^S,1:nw))

{#IFDEF STAGGERED
staggered : associate(wFis=>sFi%ws%w,wCos=>sCo%ws%w)
}


{#IFNDEF DY_SP
if (covariant) myM=>pgeogrid%m
}

! code test erase first
{#IFNDEF DY_SP
if (amrentropy) then
   call e_to_rhos(ixFiG^L,ixFi^L,wFi,xFi)
else if (coarsenprim) then
   call primitive(ixFiG^L,ixFi^L,wFi,xFi)
end if
}

if (slab) then
   CoFiratio=one/dble(2**ndim)
   do iw=1,nw
      {do ixCo^DB = ixCo^LIM^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
         wCo(ixCo^D,iw)=sum(wFi(ixFi^D:ixFi^D+1,iw))*CoFiratio
      {end do\}
   end do
else
   do iw=1,nw
      {do ixCo^DB = ixCo^LIM^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
         wCo(ixCo^D,iw)= &
             sum(pgeogrid%dvolume(ixFi^D:ixFi^D+1)*wFi(ixFi^D:ixFi^D+1,iw)) &
            /pgeoCoarsegrid%dvolume(ixCo^D)
      {end do\}
   end do

   
{#IFDEF STAGGERED
   do iw=1,nws

       select case(iw)
       {case(^D)
       ! Start one layer before
       {do ixCo^DDB = ixComin^DDB-kr(^DDB,^D),ixComax^DDB
       ixFi^DDB=2*(ixCo^DDB-ixComin^DDB+kr(^DDB,^D))+ixFimin^DDB-kr(^DDB,^D)\}
       
            ! This if statement catches the axis where surface is zero:
            if (pgeoCoarsegrid%surfaceC^D(ixCo^DD) .gt. 1.0d-9*pgeoCoarsegrid%dvolume(ixCo^DD)) then ! Normal case
              wCos(ixCo^DD,iw)= &
                sum(pgeogrid%surfaceC^D(ixFi^DD:ixFi^DD+1-kr(^D,^DD))*wFis(ixFi^DD:ixFi^DD+1-kr(^D,^DD),iw)) &
                 /pgeoCoarsegrid%surfaceC^D(ixCo^DD)
            else ! On axis
               wCos(ixCo^DD,iw)= zero
            end if

       {end do\}
       }
       end select

   end do
   ! average to fill cell-centred values
   call faces2centers(ixCo^L,sCo)
   }
   
end if


! code test erase first
{#IFNDEF DY_SP
if (amrentropy) then
   if (keepFi) then
      {#IFNDEF DY_SP
      if (covariant) myM=>pgeogrid%m
      }
      call rhos_to_e(ixFiG^L,ixFi^L,wFi,xFi)
   end if
   {#IFNDEF DY_SP
   if (covariant) myM=>pgeoCoarsegrid%m
   }
   call rhos_to_e(ixCoG^L,ixCo^L,wCo,xCo)
else if (coarsenprim) then
   if (keepFi) then
      {#IFNDEF DY_SP
      if (covariant) myM=>pgeogrid%m
      }
      call conserve(ixFiG^L,ixFi^L,wFi,xFi,patchfalse)
   end if
   {#IFNDEF DY_SP
   if (covariant) myM=>pgeoCoarsegrid%m
   }
   call conserve(ixCoG^L,ixCo^L,wCo,xCo,patchfalse)
end if
}

{#IFDEF STAGGERED
end associate staggered
}


end associate
end subroutine coarsen_grid
!=============================================================================
