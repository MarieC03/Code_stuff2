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
subroutine setdt

! setdt  - set dt for all levels between levmin and levmax. 
!          dtpar>0  --> use fixed dtpar for all level
!          dtpar<=0 --> determine CFL limited timestep 
!        - set dtimpl

include 'amrvacdef.f'

integer :: iigrid, igrid, ncycle, ncycle2, ifile
double precision :: dtnew, qdtnew, dtmin_mype, factor, dx^D, dxmin^D
double precision :: qdtimpl,dtimpl_mype

double precision :: dtmax, dxmin
integer,save :: stepflag
!----------------------------------------------------------------------------


if(it==0) stepflag = 0
{#IFDEF TCRKL2
if(conduction) then
   dtimpl_mype=bigdouble
!$OMP PARALLEL DO PRIVATE(igrid,qdtnew,qdtimpl,&
!$OMP& dx^D) REDUCTION(min:dtimpl_mype)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      qdtimpl=bigdouble
      dx^D=rnode(rpdx^D_,igrid);
      saveigrid = igrid
      if (.not.slab) mygeo => pgeo(igrid)
      if (B0field) myB0 => pB0_cell(igrid)
      {#IFNDEF DY_SP
        if (covariant) myM => mygeo%m
      }
      qdtnew=bigdouble
{^IFMHDPHYS
      call getdt_heatconduct_mhd(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,px(igrid)%x)
}
{^IFHDPHYS
      call getdt_heatconduct_hd(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,px(igrid)%x)
}
      qdtimpl=min(qdtimpl,qdtnew)
      dtimpl_mype=min(dtimpl_mype,qdtimpl)
  end do
!$OMP END PARALLEL DO
   call MPI_ALLREDUCE(dtimpl_mype,dtimpl,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                         icomm,ierrmpi)
endif
}
if(sourceimpl) then
   dtimpl_mype=bigdouble
!$OMP PARALLEL DO PRIVATE(igrid,qdtnew,qdtimpl,&
!$OMP& dx^D) REDUCTION(min:dtimpl_mype)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      qdtimpl=bigdouble
      dx^D=rnode(rpdx^D_,igrid);
      saveigrid = igrid
      if (.not.slab) mygeo => pgeo(igrid)
      if (B0field) myB0 => pB0_cell(igrid)
      {#IFNDEF DY_SP
        if (covariant) myM => mygeo%m
      }
      call getdt_impl(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,px(igrid)%x)
      qdtimpl=min(qdtimpl,qdtnew)
      dtimpl_mype=min(dtimpl_mype,qdtimpl)
   end do
!$OMP END PARALLEL DO
   call MPI_ALLREDUCE(dtimpl_mype,dtimpl,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                         icomm,ierrmpi)
endif

if (dtpar<=zero) then
   dtmin_mype=bigdouble
   cmax_mype = zero
!$OMP PARALLEL DO PRIVATE(igrid,qdtnew,dtnew,dx^D) REDUCTION(min:dtmin_mype)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dtnew=bigdouble
      dx^D=rnode(rpdx^D_,igrid);
      call set_tmpGlobals(igrid)

      call getdt_courant(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,px(igrid)%x)

      dtnew=min(dtnew,qdtnew)
      
      call getdt(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,px(igrid)%x)
      
      dtnew=min(dtnew,qdtnew)
      
      call getdt_special(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,px(igrid)%x)
      
      dtnew=min(dtnew,qdtnew)
      dtmin_mype=min(dtmin_mype,dtnew)
      
      dt_grid(igrid)=dtnew
      
   end do
!$OMP END PARALLEL DO
else
   dtmin_mype=dtpar
end if

if (dtmin_mype<dtmin) then
   write(unitterm,*)"Warning: Time step too small!", dtmin_mype
   write(unitterm,*)"on processor:", mype
   write(unitterm,*)"at time:", t," step:", it
   call mpistop("too small timestep")
end if

if (slowsteps>it-itmin+1) then
   factor=one-(one-dble(it-itmin+1)/dble(slowsteps))**2
   if (time_accurate) then
      dtmin_mype=dtmin_mype*factor
   else
!$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
         dt_grid(igrid)=dt_grid(igrid)*factor
      end do
!$OMP END PARALLEL DO
   end if
end if


if (time_accurate) then

   if( stepflag<1.and.mype==0) then
      if(any(dtsave(1:nfile)<dtmin_mype )) then
         write(unitterm,1001) dtmin_mype, dtsave(1:nfile)
         stepflag = 1     
      endif
   endif   
 
   if (tmaxexact) then
      dtmin_mype=min(dtmin_mype,tmax-t)
   end if

   if(any(dtsave(1:nfile)<bigdouble).or.any(tsave(isavet(1:nfile),1:nfile)<bigdouble))then
      dtmax = minval((int(t/dtsave(1:nfile))+1)*dtsave(1:nfile))-t
      do ifile=1,nfile
         dtmax = min(tsave(isavet(ifile),ifile)-t,dtmax)
      end do
      if(dtmax<dtmin_mype .and. dtmax > smalldouble)then 
        dtmin_mype=min(dtmin_mype,dtmax)
      end if      
   end if

   if (dtpar<=zero) then
      call MPI_ALLREDUCE(dtmin_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                         icomm,ierrmpi)
   else
      dt=dtmin_mype
   end if
   
   if(sourceimpl)then
   if(sourceimplcycle)then
     ncycle2=ceiling(dt/dtimpl)
     if(ncycle2<=1)then
     !  write(unitterm,*)'implicit time step not smaller than CFL limit'
     !  write(unitterm,*)'dt =',dt,' dtimpl=',dtimpl
     !  call mpistop("no need for implicit cycling")
     else
       if(sourceparasts)then
         ncycle=floor(dsqrt(dble(ncycle2)))
         !!if(ncycle*dtimpl>dt)call mpistop("increased timestep!")
         dt=ncycle*dtimpl 
       else
         ncycle=ncycle2
       endif
       if (ncycle>ncyclemax) then
          if(mype==0) &
           print *,'too many subcycles, reducing dt to',ncyclemax,'dt_impl!!'
          dt=ncyclemax*dtimpl
       endif
     endif
   endif
   endif
{#IFDEF TCRKL2
   if(conduction) then
     ncycle=ceiling(dt/dtimpl)
     if (ncycle>ncyclemax) then
       if(mype==0 .and. .true.) then
         write(*,*) 'CLF time step is too many times larger than conduction time step',ncycle
         write(*,*) 'reducing dt to',ncyclemax,'times of dt_impl!!'
       endif
       dt=ncyclemax*dtimpl
     endif
   endif
}


!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dt_grid(igrid)=dt
   end do
!$OMP END PARALLEL DO
     

! global Lax-Friedrich finite difference flux splitting needs fastest wave-speed
! so does GLM: 
call MPI_ALLREDUCE(cmax_mype,cmax_global,1,&
     MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)

end if


1001 format(' Warning: timesteps: ',1x,1pe12.5,' exceeding output intervals ',2(1x,1pe12.5))

end subroutine setdt
!=============================================================================
subroutine getdt_courant(w,ixG^L,ix^L,dtnew,dx^D,x)

! compute CFL limited dt (for variable time stepping)
use mod_cfc_parameters
include 'amrvacdef.f'
 
integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

integer :: idims
double precision :: courantmax, dxinv(1:ndim), courantmaxtot, courantmaxtots
double precision :: cmax(ixG^T), cmaxtot(ixG^T), tmp(ixG^T)
!-----------------------------------------------------------------------------
dtnew=bigdouble

courantmax=zero
courantmaxtot=zero
courantmaxtots=zero



^D&dxinv(^D)=one/dx^D;

cmaxtot(ix^S)=zero

do idims=1,ndim
   call getcmax(w,x,ixG^L,ix^L,idims,cmax)
   !call getcmax_old(w,x,ixG^L,ix^L,idims,cmax,tmp,.false.)
   cmax_mype = max(cmax_mype,maxval(cmax(ix^S)))
{#IFDEF DY_SP
   if (.not.slab .and. (coordinate .ne. Cartesian) ) then
}
{#IFNDEF DY_SP
   if (.not.slab) then
}
      tmp(ix^S)=cmax(ix^S)/mygeo%dx(ix^S,idims)
      cmaxtot(ix^S)=cmaxtot(ix^S)+tmp(ix^S)
      courantmax=max(courantmax,maxval(tmp(ix^S)))
   else
      cmaxtot(ix^S)=cmaxtot(ix^S)+cmax(ix^S)*dxinv(idims)
      courantmax=max(courantmax,maxval(cmax(ix^S)*dxinv(idims)))
   end if
   courantmaxtot=courantmaxtot+courantmax
end do


select case (typecourant)
case ('minimum')
   ! courantmax='max(c/dx)'
   if (courantmax>smalldouble)     dtnew=min(dtnew,courantpar/courantmax)
case ('summax')
   ! courantmaxtot='summed max(c/dx)'
   if (courantmaxtot>smalldouble)  dtnew=min(dtnew,courantpar/courantmaxtot)
case ('maxsum')
   ! courantmaxtots='max(summed c/dx)'
   courantmaxtots=max(courantmaxtots,maxval(cmaxtot(ix^S)))
   if (courantmaxtots>smalldouble) dtnew=min(dtnew,courantpar/courantmaxtots)
case default
   write(unitterm,*)'Unknown typecourant=',typecourant
   call mpistop("Error from getdt_courant: no such typecourant!")
end select


end subroutine getdt_courant
!=============================================================================
subroutine getresidual(iit)

! compute residual for steady state calculations

include 'amrvacdef.f'
 
integer, intent(in) :: iit

integer :: iigrid,igrid,iw
double precision :: wnrm2_send(1:nwflux),wnrm2_recv(1:nwflux)
double precision :: wnrm2localgrids(1:nwflux),residlocalgrids(1:nwflux)
double precision :: resid_send(1:nwflux),resid_recv(1:nwflux)
!----------------------------------------------------------------------------

select case(typeresid)
  case('relative')
    wnrm2localgrids(1:nwflux)=zero
!$OMP PARALLEL DO PRIVATE(igrid,iw) REDUCTION(+:wnrm2localgrids)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      do iw=1,nwflux
        wnrm2localgrids(iw)=wnrm2localgrids(iw)+sum(pw(igrid)%w(ixM^T,iw)**2)
      enddo
    end do 
!$OMP END PARALLEL DO

    wnrm2_send(1:nwflux)=wnrm2localgrids(1:nwflux)

    call MPI_ALLREDUCE(wnrm2_send,wnrm2_recv,nwflux,MPI_DOUBLE_PRECISION, &
               MPI_SUM,icomm,ierrmpi)

    do iw=1,nwflux
      if(wnrm2_recv(iw)<smalldouble) wnrm2_recv(iw)=one
    enddo
  case('absolute')
    wnrm2_recv(1:nwflux)=one
  case default
    call mpistop('no such typeresid')
end select

residlocalgrids(1:nwflux)=zero
!$OMP PARALLEL DO PRIVATE(igrid,iw) REDUCTION(+:residlocalgrids)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
    do iw=1,nwflux
      residlocalgrids(iw)=residlocalgrids(iw) &
           +sum(pwres(igrid)%w(ixM^T,iw)**2)/wnrm2_recv(iw)
    enddo
end do 
!$OMP END PARALLEL DO

resid_send(1:nwflux)=residlocalgrids(1:nwflux)

call MPI_ALLREDUCE(resid_send,resid_recv,nwflux,MPI_DOUBLE_PRECISION, &
            MPI_SUM,icomm,ierrmpi)


residual=sqrt(sum(resid_recv(1:nwflux))/(nwflux))


end subroutine getresidual
!=============================================================================
