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

{#IFDEF STAGGERED
!=============================================================================
subroutine recalculateB
use mod_interpolate
include 'amrvacdef.f'


! ... local ...
integer :: igrid,iigrid
{#IFDEF DY_SP
double precision                  :: psi6(ixG^T,1:igridstail)
double precision                  :: psi6_surfC^D(ixG^T,1:igridstail)
double precision                  :: xi(ixG^T,1:ndim,1:igridstail),wi(ixG^T,1:nw,1:igridstail)
integer                           :: idims, ixC^LL
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
        ! Before this, you need psi
        do idims= 1, ^ND
            xi(ixG^T,1:ndim,iigrid) = px(igrid)%x(ixG^T,1:ndim)
            xi(ixG^T,idims,iigrid)  = px(igrid)%x(ixG^T,idims) + 0.5d0* pgeo(igrid)%dx(ixG^T,idims)

            call metric_interpolation(ixG^LL,ixG^LL,idims,pw(igrid)%w,px(igrid)%x,&
                 wi(ixG^T,1:nw,iigrid),xi(ixG^T,1:ndim,iigrid))
            select case (idims)
            {case (^D)
                psi6_surfC^D(ixG^T,iigrid) = wi(ixG^T, psi_metric_, iigrid)**6
                pgeo(igrid)%surfaceC^D(ixG^T)  = pgeo(igrid)%surfaceC^D(ixG^T) * psi6_surfC^D(ixG^T,iigrid)
            \}
            end select
        enddo
        psi6(ixG^T,iigrid) = pw(igrid)%w(ixG^T,psi_metric_)**6
        pgeo(igrid)%dvolume(ixG^T) = pgeo(igrid)%dvolume(ixG^T) * psi6(ixG^T,iigrid)
end do
!$OMP END PARALLEL DO
}

call init_comm_fix_conserve(1,^ND)

!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   ! Make zero the magnetic fluxes
   ! Fake advance, storing electric fields at edges
   call fake_advance(igrid,1,^ND,ps(igrid))

end do
!$OMP END PARALLEL DO

! Do correction

do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call sendflux(igrid,1,^ND) ! OMP: not threadsafe!
end do

call fix_conserve(ps,1,^ND)

call fix_edges(ps,1,^ND)

{#IFDEF DY_SP
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
    pgeo(igrid)%dvolume(ixG^T)           = pgeo(igrid)%dvolume(ixG^T) / psi6(ixG^T,iigrid)
    {^D&  pgeo(igrid)%surfaceC^D(ixG^T)  = pgeo(igrid)%surfaceC^D(ixG^T) / psi6_surfC^D(ixG^T,iigrid) \}

    do idims=1,ndim
      !ixClo^D=ixGlo^D;
      !ixChi^D=ixGhi^D-kr(idims,^D);
      xi(ixG^T,1:ndim,iigrid) = px(igrid)%x(ixG^T,1:ndim)
      xi(ixG^T,idims,iigrid)=xi(ixG^T, idims,iigrid)+0.5d0*pgeo(igrid)%dx(ixG^T, idims)
      call metric_interpolation(ixG^LL,ixG^LL,idims,pw(igrid)%w,px(igrid)%x,&
                                wi(ixG^T,1:nw,iigrid),xi(ixG^T,1:ndim,iigrid))
      select case (idims)
      {case (^D)
          psi6_surfC^D(ixG^T,iigrid) = wi(ixG^T, psi_metric_, iigrid)**6
          pws(igrid)%w(ixG^T, idims) = pws(igrid)%w(ixG^T, idims) * psi6_surfC^D(ixG^T,iigrid)
      \}
      end select
    enddo
end do
!$OMP END PARALLEL DO
}


! Now we fill the centers for the staggered variables
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call faces2centers(ixG^LL,ps(igrid))
end do
!$OMP END PARALLEL DO

end subroutine recalculateB
!=============================================================================
subroutine fake_advance(igrid,idim^LIM,s)

include 'amrvacdef.f'

integer       :: igrid,idim^LIM
type(state)   :: s
! ... local ...
double precision             :: dx^D
double precision             :: fC(ixG^T,1:nwflux,1:ndim)
double precision             :: fE(ixG^T,1:^NC)
!-----------------------------------------------------------------------------

dx^D=rnode(rpdx^D_,igrid);
call set_tmpGlobals(igrid)

call fake_update(ixG^LL,s,fC,fE,dx^D)

call storeflux(igrid,fC,idim^LIM)
call storeedge(igrid,ixG^LL,fE,idim^LIM) 

end subroutine fake_advance
!=============================================================================
subroutine fake_update(ixI^L,s,fC,fE,dx^D)
! In reality the same as b_from_vectorpotential for staggered case

include 'amrvacdef.f'

integer       :: ixI^L
type(state)   :: s
double precision             :: fC(ixI^S,1:nwflux,1:ndim)
double precision             :: fE(ixG^T,1:^NC)
double precision             :: dx^D
! ... local ...
integer                            :: ixIs^L,ixO^L,idir
double precision                   :: xC(ixI^S,1:ndim), A(ixI^S,1:ndir)
double precision                   :: circ(ixI^S,1:ndim), dxidir
!-----------------------------------------------------------------------------
associate(ws=>s%ws%w,x=>s%x%x)

A(:^D&,:)=zero
ws(:^D&,:)=zero

ixIs^L=s%ws%ixG^L;
ixO^L=ixI^L^LSUBdixB;

call b_from_vectorpotentialA(ixIs^L, ixI^L, ixO^L, ws, x, A)

! This is important only in 3D

do idir=1,ndim
   fE(:^D&,idir) =-A(:^D&,idir)*dxlevel(idir)
end do

end associate
end subroutine fake_update
!=============================================================================
}{#IFNDEF STAGGERED
subroutine fct_average(ixI^L, ixO^L, fC)

! Performs the average for the flux constrained transport from Toth 2000. 
! His equation (25)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim)

double precision                   :: fB(ixI^S,1:ndir,1:ndim)
integer                            :: iwdim, iwdir, idim, idir
integer                            :: ixJp^L, ixKp^L, ixKm^L, ixC^L, ixJpKm^L
!-----------------------------------------------------------------------------


do idim = 1, ndim ! idim is the direction of the flux interface
   do idir = 1, ndim ! idir is the component of the field to consider
      if (idim.eq.idir) cycle
      iwdir = b0_+idir; iwdim = b0_+idim

! Assemble indices:
      {ixCmax^D=ixOmax^D+1;}
      {ixCmin^D=ixOmin^D-1-kr(idim,^D);} ! Extend range in flux direction
 
      {ixJp^L=ixC^L+kr(idim,^D);}
      {ixKp^L=ixC^L+kr(idir,^D);}
      {ixKm^L=ixC^L-kr(idir,^D);}
      {ixJpKm^L=ixJp^L-kr(idir,^D);}


! Perform flux average:
         fB(ixC^S,idir,idim) = &
              0.125d0 * ( 2.0d0* fC(ixC^S,iwdir,idim) &
              + fC(ixKp^S,iwdir,idim) &
              + fC(ixKm^S,iwdir,idim) - dxlevel(idir) &
              * (fC(ixC^S,iwdim,idir)  &
              + fC(ixJp^S,iwdim,idir)  &
              + fC(ixKm^S,iwdim,idir)  &
              + fC(ixJpKm^S,iwdim,idir)  &
              )/dxlevel(idim))

   end do ! idir
end do ! idim

! Overwrite the new flux entries:
do idim = 1, ndim ! idim is the direction of the flux interface
   do idir = 1, ndim ! idir is the component of the field to consider
      iwdir = b0_+idir
      if (idim.eq.idir) then
         ! To conserve divb to machine precision, this needs to be zero.
         ! However, since restriction and prolongation can introduce divb
         ! when using AMR, additional measures for divb-control must be taken. 
         ! These rely on the normal field flux component (e.g. as in GLM).
         ! Thus when more than one level is used, we dont set this zero 
         ! to allow additional divb-control to work.
         if (mxnest .eq. 1) fC(ixI^S,iwdir,idim) = zero
         cycle
      end if

      fC(ixI^S,iwdir,idim) = fB(ixI^S,idir,idim)

   end do ! idir
end do ! idim

end subroutine fct_average
}
!=============================================================================
{#IFDEF STAGGERED
subroutine faces2centers(ixO^L,s)

  include 'amrvacdef.f'

  ! Non-staggered interpolation range
  integer, intent(in)                :: ixO^L
  type(state)                        :: s

  ! --- local ---
  integer                            :: hxO^L, idim, ix^D
  !-----------------------------------------------------------------------------
  associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)

  do idim=1,ndim
     ! Displace index to the left
     ! Even if ixI^L is the full size of the w arrays, this is ok
     ! because the staggered arrays have an additional place to the left.

     hxO^L=ixO^L-kr(idim,^D);

     ! Interpolate to cell barycentre using arithmetic average
     ! This might be done better later, to make the method less diffusive.
     select case(idim)
        {case(^D)
        w(ixO^S,b0_+idim)=(half*mygeo%dx(ixO^S,idim)/mygeo%dvolume(ixO^S))*&
             (ws(ixO^S,bs0_+idim)*mygeo%surfaceC^D(ixO^S)&
             +ws(hxO^S,bs0_+idim)*mygeo%surfaceC^D(hxO^S))
        {#IFDEF DY_SP
          w(ixO^S,b0_+idim) = w(ixO^S,b0_+idim) / w(ixO^S,psi_metric_)**6
        }
        \}
     end select
  end do

  end associate

end subroutine faces2centers
!=============================================================================
subroutine faces2centers4(ixO^L,s)

  include 'amrvacdef.f'

  ! Non-staggered interpolation range
  ! Fourth order equation.
  integer, intent(in)                :: ixO^L
  type(state)                        :: s

  ! --- local ---
  integer                            :: gxO^L, hxO^L, jxO^L, idim
  !-----------------------------------------------------------------------------
  associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)

  call mpistop('Stop using faces2centers4')

  do idim=1,ndim

     gxO^L=ixO^L-2*kr(idim,^D);
     hxO^L=ixO^L-kr(idim,^D);
     jxO^L=ixO^L+kr(idim,^D);

     ! Interpolate to cell barycentre using fourth order central formula
     select case(idim)
        {case(^D)     
        w(ixO^S,b0_+idim)=(1.0d0/16.0d0/mygeo%surface^D(ixO^S)) &
             *( &
             -ws(gxO^S,bs0_+idim)*mygeo%surfaceC^D(gxO^S) &
             +9.0d0*ws(hxO^S,bs0_+idim)*mygeo%surfaceC^D(hxO^S) &
             +9.0d0*ws(ixO^S,bs0_+idim)*mygeo%surfaceC^D(ixO^S) &
             -ws(jxO^S,bs0_+idim)*mygeo%surfaceC^D(jxO^S) &
             )
        \}
     end select
  end do

  end associate

end subroutine faces2centers4
!=============================================================================
subroutine faces2centers6(ixO^L,s)

  include 'amrvacdef.f'

  ! Non-staggered interpolation range
  ! Sixth order equation.
  integer, intent(in)                :: ixO^L
  type(state)                        :: s

  ! --- local ---
  integer                            :: fxO^L, gxO^L, hxO^L, jxO^L, kxO^L, idim
  !-----------------------------------------------------------------------------
  associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)

  call mpistop('Stop using faces2centers6')
  do idim=1,ndim

     fxO^L=ixO^L-3*kr(idim,^D);
     gxO^L=ixO^L-2*kr(idim,^D);
     hxO^L=ixO^L-kr(idim,^D);
     jxO^L=ixO^L+kr(idim,^D);
     kxO^L=ixO^L+2*kr(idim,^D);

     ! Interpolate to cell barycentre using sixth order central formula
     select case(idim)
        {case(^D)     
        w(ixO^S,b0_+idim)=(1.0d0/256.0d0/mygeo%surface^D(ixO^S)) &
             *( &
             +3.0d0*ws(fxO^S,bs0_+idim)*mygeo%surfaceC^D(fxO^S) &
             -25.0d0*ws(gxO^S,bs0_+idim)*mygeo%surfaceC^D(gxO^S) &
             +150.0d0*ws(hxO^S,bs0_+idim)*mygeo%surfaceC^D(hxO^S) &
             +150.0d0*ws(ixO^S,bs0_+idim)*mygeo%surfaceC^D(ixO^S) &
             -25.0d0*ws(jxO^S,bs0_+idim)*mygeo%surfaceC^D(jxO^S) &
             +3.0d0*ws(kxO^S,bs0_+idim)*mygeo%surfaceC^D(kxO^S) &
             )
        \}
     end select
  end do

  end associate

end subroutine faces2centers6
!=============================================================================
subroutine show_div_staggered(ixO^L,s)
! For calculating the divergence of a face-allocated vector field.
include 'amrvacdef.f'

integer, intent(in)           :: ixO^L
type(state)                   :: s
! ... local ...
integer                       :: hxO^L,idim
double precision              :: out(ixG^T)
!----------------------------------------------------------------------------
print *, 'In div_stg'

print *, 'x1 min = ', s%x%x(ixMlo^D,1)
print *, 'x2 min = ', s%x%x(ixMlo^D,2)

associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)
out=zero

do idim=1,ndim
   ! Displace index to the left
   ! Even if ixI^L is the full size of the w arrays, this is ok
   ! because the staggered arrays have an additional place to the left.

   hxO^L=ixO^L-kr(idim,^D);
!   print *, 'ixO^L', ixO^L
!   print *, 'hxO^L', hxO^L
   ! Calculate divergence by taking differences of fluxes
   select case(idim)
   {case(^D)
   out(ixO^S)=out(ixO^S)+(ws(ixO^S,bs0_+idim)*mygeo%surfaceC^D(ixO^S)&
                        -ws(hxO^S,bs0_+idim)*mygeo%surfaceC^D(hxO^S))&
                        /mygeo%dvolume(ixO^S)
      {#IFDEF DY_SP
        out(ixO^S) = out(ixO^S) / w(ixO^S, psi_metric_)**6
      }
   \}
   end select
end do
end associate

call printarray(ixG^LL,ixO^L,out)

end subroutine show_div_staggered
!---------------------------------------------------------------------------
subroutine div_staggered(ixI^L,ixO^L,s,out)
! For calculating the divergence of a face-allocated vector field.
include 'amrvacdef.f'

integer, intent(in)           :: ixI^L, ixO^L
type(state)                   :: s
double precision              :: out(ixI^S)
! ... local ...
integer                       :: hxO^L,idim
!----------------------------------------------------------------------------

associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)
out(ixO^S)=zero

do idim=1,ndim
   ! Displace index to the left
   ! Even if ixI^L is the full size of the w arrays, this is ok
   ! because the staggered arrays have an additional place to the left.

   hxO^L=ixO^L-kr(idim,^D);
!   print *, 'ixO^L', ixO^L
!   print *, 'hxO^L', hxO^L
   ! Calculate divergence by taking differences of fluxes
   select case(idim)
   {case(^D)
   out(ixO^S)=out(ixO^S)+(ws(ixO^S,bs0_+idim)*mygeo%surfaceC^D(ixO^S)&
                        -ws(hxO^S,bs0_+idim)*mygeo%surfaceC^D(hxO^S))&
      {#IFNDEF DY_SP
                        /mygeo%dvolume(ixO^S)
      }
      {#IFDEF DY_SP
                        !/ (mygeo%dvolume(ixO^S))
                        / (mygeo%dvolume(ixO^S) *  w(ixO^S, psi_metric_)**6)
      }
   \}
   end select
end do
end associate

!call printarray(ixG^LL,ixO^L,out)

end subroutine div_staggered
!=============================================================================
subroutine updatefaces(ixI^L,ixO^L,qdt,fC,fE,s)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: qdt
type(state)                        :: s
!double precision, intent(in)       :: w(ixI^S,1:nw)
!double precision, intent(in)       :: x(ixI^S)
double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
double precision, intent(inout)    :: fE(ixI^S,1:^NC)
!double precision, intent(inout)    :: bfaces(ixIs^S,1:ndir)

! --- local ---
integer                            :: ix^L(1:ndim)
integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
integer                            :: idim1,idim2,idir,iwdim1,iwdim2,i,j,k
!double precision                   :: edge(ixI^S,1:ndir)
double precision                   :: circ(ixI^S,1:ndim)
!-----------------------------------------------------------------------------
!associate(mygeo=>s%geo,bfaces=>s%ws%w)
associate(bfaces=>s%ws%w,x=>s%x%x)

! Calculate contribution to FEM of each edge,
! that is, estimate value of line integral of
! electric field in the positive idir direction.

ixCmax^D=ixOmax^D;
ixCmin^D=ixOmin^D-1;

fE(ixI^S,1:ndir)=zero

do idim1=1,ndim 
   iwdim1 = b0_+idim1
   do idim2=1,ndim
      iwdim2 = b0_+idim2
      do idir=1,ndir ! Direction of line integral
         ! Allow only even permutations
         if (lvc(idim1,idim2,idir).eq.1) then
            ! Assemble indices
            jxC^L=ixC^L+kr(idim1,^D);
            ixCp^L=ixC^L+kr(idim2,^D);
            ! Interpolate to edges
            fE(ixC^S,idir)=qdt*quarter*((fC(ixC^S,iwdim1,idim2)&
            +fC(jxC^S,iwdim1,idim2))/dxlevel(idim1)&
            -(fC(ixC^S,iwdim2,idim1)+fC(ixCp^S,iwdim2,idim1))&
            /dxlevel(idim2))

            if (typeaxial.ne.'slab' .and. typeaxial.ne.'cartesian') then
{^IFZIN
{#IFNDEF HARDBC
            ! Catch axis and origin, where the line integral is performed on a
            ! zero lenght element and should be zero.
            ! Remember that staggered quantities are assigned to the index
            ! of the cell at their left.
            where((abs(x(ixC^S,z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
                  (abs(x(ixC^S,z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))
}{#IFDEF HARDBC
            where((abs(x(ixC^S,z_)+half*dxlevel(z_)-xprobmin^Z).lt.smalldouble).or.&
                  (abs(x(ixC^S,z_)+half*dxlevel(z_)-xprobmax^Z).lt.smalldouble))
}
              fE(ixC^S,^PHI)=zero
              fE(ixC^S,r_)=zero
            end where}
            where(abs(x(ixC^S,r_)+half*dxlevel(r_)).lt.1.0d-9)
              fE(ixC^S,idir)=zero
            end where
            end if
         end if
      end do
   end do
end do

circ(ixI^S,1:ndim)=zero

! Calculate circulation on each face

do idim1=1,ndim ! Coordinate perpendicular to face 
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
        ! Assemble indices
        hxC^L=ixC^L-kr(idim2,^D);
        ! Add line integrals in direction idir
        circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixC^S,idir)&
                         -fE(hxC^S,idir))
      end do
   end do
end do

! Decrease bottom limit by one
do idim1=1, ndim
   ixmax^D(idim1)=ixOmax^D;ixmin^D=ixOmin^D-1;!-kr(^D,idim1);
end do

! Divide by the area of the face to get dB/dt

do idim1=1,ndim
   ixC^L=ix^L(idim1);
   select case(idim1)
   {case(^D)
   where(mygeo%surfaceC^D(ixC^S) .gt. 1.0d-9*mygeo%dvolume(ixC^S))
      circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                 /mygeo%surfaceC^D(ixC^S)
   elsewhere
      circ(ixC^S,idim1)=zero
   end where
   \}
    end select

   ! Time update
   ! minus!
   bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
end do

end associate

end subroutine updatefaces
!=============================================================================
! UCT contact, cheaper than UCT hll, since without reconstruction, mainly using quantites at cell centers
!> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
! the fC here without surfaceC
subroutine updatefacescontact(ixI^L,ixO^L,qdt,vnorm,fC,fE,wprim,s)
  use mod_imhd_intermediate
  include 'amrvacdef.f'

  integer, intent(in)                :: ixI^L, ixO^L
  double precision, intent(in)       :: qdt
  type(state)                        :: s
  double precision, intent(in)       :: wprim(ixI^S,1:nw)
  double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
  double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)
  double precision, intent(in)       :: vnorm(ixI^S,1:ndim)

  ! v_hat^i at cell centers
  double precision                   :: v_hat(ixI^S,1:ndir)
  ! electric field at cell centers
  double precision                   :: ECC(ixI^S,7-2*ndim:3)
  ! gradient of E at left and right side of a cell face
  double precision                   :: EL(ixI^S),ER(ixI^S)
  ! gradient of E at left and right side of a cell corner
  double precision                   :: ELC(ixI^S),ERC(ixI^S)
  double precision                   :: lfac(ixI^S), dxidir
  double precision                   :: sqrtg(s%ws%ixG^S)
  double precision                   :: sqrt_gamma_hat(s%ws%ixG^S)
  double precision                   :: xC(s%ws%ixG^S,1:ndim),xCC(s%ws%ixG^S,1:ndim)
  integer                            :: hxC^L,ixC^L,jxC^L,ixA^L,ixB^L
  integer                            :: idim,idim1,idim2,idir,iwdim1,iwdim2

  associate(bfaces=>s%ws%w,x=>s%x%x)


  {#IFDEF DY_SP
  call imhd_get_intermediate_variables(ixI^L, ixI^L, wprim, x, lfac=lfac(ixI^S))

    ! cell center v_hat
    {^C& v_hat(ixI^S,^C) = wprim(ixI^S, alp_metric_) * wprim(ixI^S, u^C_) / lfac(ixI^S) &
                      - wprim(ixI^S, beta_metric^C_) \}
  }
  {#IFNDEF DY_SP
    lfac(ixI^S) = wprim(ixI^S,lfac_)
    {^C& v_hat(ixI^S,^C) = mygeo%m%alpha(ixI^S) * wprim(ixI^S, u^C_) / lfac(ixI^S) &
                      - mygeo%m%beta(^C)%elem(ixI^S) \}
  }

  ixGs^L=s%ws%ixG^L;
  
  ! extend one layer of cell center locations in xCC
  xCC=0.d0
  xCC(ixI^S,1:ndim)=x(ixI^S,1:ndim)
  {
  xCC(ixGsmin^D^D%ixI^S,1:ndim) = x(ixImin^D^D%ixI^S,1:ndim)
  xCC(ixGsmin^D^D%ixGs^S,^D) = x({ixImin^DD,},^D) - s%geo%dx({ixImin^DD,},^D)
  \}
  {^IFTHREED
  xCC(ixImin1:ixImax1,ixGsmin2,ixGsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,1)
  xCC(ixGsmin1,ixImin2:ixImax2,ixGsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,2)
  xCC(ixGsmin1,ixGsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,3)
  }


  {^NOONED
  ! Calculate electric field E_hat_i at cell centers,
  ! with sqrt_gamma_hat factored out
  do idir=7-2*ndim,3
    idim1=mod(idir,ndir)+1   ! 'Next' direction
    idim2=mod(idir+1,ndir)+1 ! 'Next Next' direction
    ECC(ixI^S,idir) = - v_hat(ixI^S,idim1) * wprim(ixI^S,b0_+idim2) &
                      + v_hat(ixI^S,idim2) * wprim(ixI^S,b0_+idim1)
    {#IFDEF DY_SP
    ! conformal factor
    ECC(ixI^S,idir) = wprim(ixI^S,psi_metric_)**6 * ECC(ixI^S,idir)
    }
  end do

  ! Calculate contribution to FEM of each edge,
  ! that is, estimate value of line integral of
  ! electric field in the positive idir direction.
  fE=0.0d0
  ! evaluate electric field along cell edges according to equation (41)
  do idir=7-2*ndim,3 ! Direction of line integral

    idim1=mod(idir,ndir)+1   ! 'Next' direction
    idim2=mod(idir+1,ndir)+1 ! 'Next Next' direction
    iwdim1 = b0_+idim1
    iwdim2 = b0_+idim2

    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D+kr(idir,^D)-1;
    ! Assemble indices
    jxC^L=ixC^L+kr(idim1,^D);
    hxC^L=ixC^L+kr(idim2,^D);

    ! Interpolate sqrt gamma to the edges
  {#IFDEF DY_SP
    sqrtg(ixC^S)=zero
    ! Get edge coordinates and sqrt_gamma_hat
    do idim=1,ndim
      if (idim/=idir) then
        xC(ixC^S,idim)=xCC(ixC^S,idim) + 0.5d0 * s%geo%dx(ixC^S,idim)
      else
        xC(ixC^S,idim)=xCC(ixC^S,idim)
      end if
    end do
    call get_sqrt_gamma_hat(xC(ixGs^S, 1:ndim), ixGs^L, ixC^L, sqrt_gamma_hat(ixGs^S))
  }
  
  {#IFNDEF DY_SP
    sqrtg(ixC^S)=zero
    select case(idim1)
  { case(^D)
      sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
        mygeo%mSurface^D %sqrtgamma(ixC^S)+&
        mygeo%mSurface^D %sqrtgamma(ixCp^S))
  \}
    end select
  
    select case(idim2)
  { case(^D)
      sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
        mygeo%mSurface^D %sqrtgamma(ixC^S)+&
        mygeo%mSurface^D %sqrtgamma(jxC^S))
  \}
    end select
  }

    ! average cell-face electric field to cell edges
    fE(ixC^S,idir)=0.25d0*&
    (fC(ixC^S,iwdim1,idim2)+fC(jxC^S,iwdim1,idim2)&
    -fC(ixC^S,iwdim2,idim1)-fC(hxC^S,iwdim2,idim1))

    ! add slope in idim2 direction from equation (50)
    ixAmin^D=ixCmin^D;
    ixAmax^D=ixCmax^D+kr(idim1,^D);
    EL(ixA^S)=fC(ixA^S,iwdim1,idim2)-ECC(ixA^S,idir)
    hxC^L=ixA^L+kr(idim2,^D);
    ER(ixA^S)=fC(ixA^S,iwdim1,idim2)-ECC(hxC^S,idir)
    where(vnorm(ixC^S,idim1)>0.d0)
      ELC(ixC^S)=EL(ixC^S)
    else where(vnorm(ixC^S,idim1)<0.d0)
      ELC(ixC^S)=EL(jxC^S)
    else where
      ELC(ixC^S)=0.5d0*(EL(ixC^S)+EL(jxC^S))
    end where

    hxC^L=ixC^L+kr(idim2,^D);
    where(vnorm(hxC^S,idim1)>0.d0)
      ERC(ixC^S)=ER(ixC^S)
    else where(vnorm(hxC^S,idim1)<0.d0)
      ERC(ixC^S)=ER(jxC^S)
    else where
      ERC(ixC^S)=0.5d0*(ER(ixC^S)+ER(jxC^S))
    end where
    fE(ixC^S,idir)=fE(ixC^S,idir)+0.25d0*(ELC(ixC^S)+ERC(ixC^S))

    ! add slope in idim1 direction from equation (50)
    jxC^L=ixC^L+kr(idim2,^D);
    ixAmin^D=ixCmin^D;
    ixAmax^D=ixCmax^D+kr(idim2,^D);
    EL(ixA^S)=-fC(ixA^S,iwdim2,idim1)-ECC(ixA^S,idir)
    hxC^L=ixA^L+kr(idim1,^D);
    ER(ixA^S)=-fC(ixA^S,iwdim2,idim1)-ECC(hxC^S,idir)

    where(vnorm(ixC^S,idim2)>0.d0)
      ELC(ixC^S)=EL(ixC^S)
    else where(vnorm(ixC^S,idim2)<0.d0)
      ELC(ixC^S)=EL(jxC^S)
    else where
      ELC(ixC^S)=0.5d0*(EL(ixC^S)+EL(jxC^S))
    end where

    hxC^L=ixC^L+kr(idim1,^D);
    where(vnorm(hxC^S,idim2)>0.d0)
      ERC(ixC^S)=ER(ixC^S)
    else where(vnorm(hxC^S,idim2)<0.d0)
      ERC(ixC^S)=ER(jxC^S)
    else where
      ERC(ixC^S)=0.5d0*(ER(ixC^S)+ER(jxC^S))
    end where

    fE(ixC^S,idir)=fE(ixC^S,idir)+0.25d0*(ELC(ixC^S)+ERC(ixC^S))

    ! Calculate elctric field
    if (idir .le. ndim) then
       dxidir = dxlevel(idir)
    else
       dxidir = 1.0d0
    end if

    ! times time step and edge length
    fE(ixC^S,idir)=qdt * dxidir * &
                   {#IFDEF DY_SP
                     sqrt_gamma_hat(ixC^S) * &
                   }
                   {#IFNDEF DY_SP
                     sqrtg(ixC^S) * &
                   }
                   fE(ixC^S,idir)

    ! -------------------------------------------------------------------!
    !if (typeaxial.ne.'slab') then
    if (typeaxial.ne.'slab' .and. typeaxial.ne.'cartesian') then
  {^IFZIN
  {#IFNDEF HARDBC
    ! Catch axis and origin, where the line integral is performed on a
    ! zero lenght element and should be zero.
    ! Remember that staggered quantities are assigned to the index
    ! of the cell at their left.
    where((abs(x(ixC^S,z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
          (abs(x(ixC^S,z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))
  }{#IFDEF HARDBC
    where((abs(x(ixC^S,z_)+half*dxlevel(z_)-xprobmin^Z).lt.smalldouble).or.&
          (abs(x(ixC^S,z_)+half*dxlevel(z_)-xprobmax^Z).lt.smalldouble))
  }
      fE(ixC^S,^PHI)=zero
      fE(ixC^S,r_)=zero
    end where}
    where(abs(x(ixC^S,r_)+half*dxlevel(r_)).lt.1.0d-9)
      fE(ixC^S,idir)=zero
    end where
    end if
  !  {^IFCOORDCKS
  !  ! Use cell center position to decide on e-field excision:
  !  associate(x=>s%x%x)
  !    where( ((x(ixC^S,1)**2+x(ixC^S,2)**2).lt.(eqpar(a_)+coordpar(rcut_))**2) &
  !         {#IFDEF D3 .and.(abs(x(ixC^S,3)).lt.coordpar(rcut_))} )
  !       fE(ixC^S,idir)=zero
  !    end where
  !  end associate
  !  }
  
  end do

  circ(ixI^S,1:ndim)=zero

  ! Calculate circulation on each face
  
  do idim1=1,ndim ! Coordinate perpendicular to face
     ixCmax^D=ixOmax^D;
     ixCmin^D=ixOmin^D-kr(idim1,^D);
     do idim2=1,ndim
        do idir=1,ndir ! Direction of line integral
  
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                           -fE(hxC^S,idir))
        end do
     end do
  end do
  
  ! Divide by the area of the face to get dB/dt
  do idim1=1,ndim
     ixCmax^D=ixOmax^D;
     ixCmin^D=ixOmin^D-kr(idim1,^D);
     select case(idim1)
     {case(^D)
     where(mygeo%surfaceC^D(ixC^S) .gt. 1.0d-9*mygeo%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                   /mygeo%surfaceC^D(ixC^S)
     elsewhere
        circ(ixC^S,idim1)=zero
     end where
     \}
      end select
  
     ! Time update
     ! minus!
     bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
  end do

  }

  end associate
end subroutine updatefacescontact

!=============================================================================
! fixme: only updatefaces_constact needs psi6 for electric fields
! Uct2 is called upwind constrained transport HLL
subroutine updatefacesuct2(ixI^L,ixO^L,qdt,vbarC,cbarmin,cbarmax,fE,s)


  {^IFCOORDCKS
  ! This is for excision, we can do this more generally later
  use mod_metric, only: coordpar !, rcut_
  }
include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: qdt
type(state)                        :: s
double precision, intent(in)       :: vbarC(ixI^S,1:ndir,2)
double precision, intent(in)       :: cbarmin(ixI^S,ndim)
double precision, intent(in)       :: cbarmax(ixI^S,ndim)
double precision, intent(inout)    :: fE(ixI^S,1:ndir)
! --- local ---
double precision                   :: vtilL(ixI^S,2)
double precision                   :: vtilR(ixI^S,2)
double precision                   :: btilL(s%ws%ixG^S,ndim)
double precision                   :: btilR(s%ws%ixG^S,ndim)
double precision                   :: sqrtg(s%ws%ixG^S)
double precision                   :: sqrt_gamma_hat(s%ws%ixG^S)
double precision                   :: xC(s%ws%ixG^S,1:ndim),xCC(s%ws%ixG^S,1:ndim)
double precision                   :: cp(ixI^S,2)
double precision                   :: cm(ixI^S,2)
integer                            :: ixGs^L
integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
integer                            :: idim1,idim2,idir,idim
double precision                   :: circ(ixI^S,1:ndim), dxidir
integer                            :: i^D
!-----------------------------------------------------------------------------

associate(bfaces=>s%ws%w,x=>s%x%x,mygeo=>s%geo)

! Calculate contribution to FEM of each edge,
! that is, estimate value of line integral of
! electric field in the positive idir direction.

ixGs^L=s%ws%ixG^L;

! extend one layer of cell center locations in xCC
xCC=0.d0
xCC(ixI^S,1:ndim)=x(ixI^S,1:ndim)
{
xCC(ixGsmin^D^D%ixI^S,1:ndim) = x(ixImin^D^D%ixI^S,1:ndim)
xCC(ixGsmin^D^D%ixGs^S,^D) = x({ixImin^DD,},^D) - s%geo%dx({ixImin^DD,},^D)
\}
{^IFTHREED
xCC(ixImin1:ixImax1,ixGsmin2,ixGsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,1)
xCC(ixGsmin1,ixImin2:ixImax2,ixGsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,2)
xCC(ixGsmin1,ixGsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,3)
}


! Loop over components of electric field

! idir: electric field component we need to calculate
! idim1: directions in which we already performed the reconstruction
! idim2: directions in which we perform the reconstruction

fE(ixI^S,1:ndir)=zero

do idir=7-2*ndim,ndir
  ! Indices
  ! idir: electric field component
  ! idim1: one surface
  ! idim2: the other surface
  ! cyclic permutation: idim1,idim2,idir=1,2,3
  ! Velocity components on the surface
  ! follow cyclic premutations:
  ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

  ixCmax^D=ixOmax^D;
  ixCmin^D=ixOmin^D-1+kr(idir,^D);

  ! Set indices and directions
  idim1=mod(idir,ndir)+1
  idim2=mod(idir+1,ndir)+1

  jxC^L=ixC^L+kr(idim1,^D);
  ixCp^L=ixC^L+kr(idim2,^D);

  ! Interpolate sqrt gamma to the edges
{#IFDEF DY_SP
  sqrtg(ixC^S)=zero
  ! Get edge coordinates and sqrt_gamma_hat
  do idim=1,ndim
    if (idim/=idir) then
      xC(ixC^S,idim)=xCC(ixC^S,idim) + 0.5d0 * s%geo%dx(ixC^S,idim)
    else
      xC(ixC^S,idim)=xCC(ixC^S,idim)
    end if
  end do
  call get_sqrt_gamma_hat(xC(ixGs^S, 1:ndim), ixGs^L, ixC^L, sqrt_gamma_hat(ixGs^S))
}

{#IFNDEF DY_SP
  sqrtg(ixC^S)=zero
  select case(idim1)
{ case(^D)
    sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
      mygeo%mSurface^D %sqrtgamma(ixC^S)+&
      mygeo%mSurface^D %sqrtgamma(ixCp^S))
\}
  end select
  
  select case(idim2)
{ case(^D)
    sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
      mygeo%mSurface^D %sqrtgamma(ixC^S)+&
      mygeo%mSurface^D %sqrtgamma(jxC^S))
\}
  end select
}
!checked DY_SP with correct calculation of sqrtg
!write(*,*) sqrtg(ixC^S)
!stop 'uct2'

  ! Reconstruct transverse transport velocities
  call rec(ixI^L,ixC^L,idim2,vbarC(ixI^S,idim1,1),&
           vtilL(ixI^S,2),vtilR(ixI^S,2))

  call rec(ixI^L,ixC^L,idim1,vbarC(ixI^S,idim2,2),&
           vtilL(ixI^S,1),vtilR(ixI^S,1))

  ! Reconstruct magnetic fields
  ! Eventhough the arrays are larger, rec works with
  ! the limits ixG.
  call rec(ixI^L,ixC^L,idim2,bfaces(ixI^S,idim1),&
           btilL(ixI^S,idim1),btilR(ixI^S,idim1))

  call rec(ixI^L,ixC^L,idim1,bfaces(ixI^S,idim2),&
           btilL(ixI^S,idim2),btilR(ixI^S,idim2))

  ! Take the maximum characteristic

  cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1))
  cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))

  cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2))
  cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))
 
  ! Calculate elctric field
  if (idir .le. ndim) then
     dxidir = dxlevel(idir)
  else
     dxidir = 1.0d0
  end if

  fE(ixC^S,idir)=qdt * dxidir * &
  {#IFDEF DY_SP
               sqrt_gamma_hat(ixC^S) * (&
  }
  {#IFNDEF DY_SP
               sqrtg(ixC^S) * (&
  }
               -(cp(ixC^S,1)*vtilL(ixC^S,1)*btilL(ixC^S,idim2) &
               + cm(ixC^S,1)*vtilR(ixC^S,1)*btilR(ixC^S,idim2) &
               -cp(ixC^S,1)*cm(ixC^S,1)*(btilR(ixC^S,idim2)-btilL(ixC^S,idim2)))&
               /(cp(ixC^S,1) + cm(ixC^S,1)) &
               +(cp(ixC^S,2)*vtilL(ixC^S,2)*btilL(ixC^S,idim1) &
               + cm(ixC^S,2)*vtilR(ixC^S,2)*btilR(ixC^S,idim1) &
               -cp(ixC^S,2)*cm(ixC^S,2)*(btilR(ixC^S,idim1)-btilL(ixC^S,idim1)))&
               /(cp(ixC^S,2) + cm(ixC^S,2)) &
               )

  !if (typeaxial.ne.'slab') then
  if (typeaxial.ne.'slab' .and. typeaxial.ne.'cartesian') then
{^IFZIN
{#IFNDEF HARDBC
  ! Catch axis and origin, where the line integral is performed on a
  ! zero lenght element and should be zero.
  ! Remember that staggered quantities are assigned to the index
  ! of the cell at their left.
  where((abs(x(ixC^S,z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
        (abs(x(ixC^S,z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))
}{#IFDEF HARDBC
  where((abs(x(ixC^S,z_)+half*dxlevel(z_)-xprobmin^Z).lt.smalldouble).or.&
        (abs(x(ixC^S,z_)+half*dxlevel(z_)-xprobmax^Z).lt.smalldouble))
}
    fE(ixC^S,^PHI)=zero
    fE(ixC^S,r_)=zero
  end where}
  where(abs(x(ixC^S,r_)+half*dxlevel(r_)).lt.1.0d-9)
    fE(ixC^S,idir)=zero
  end where
  end if
!  {^IFCOORDCKS
!  ! Use cell center position to decide on e-field excision:
!  associate(x=>s%x%x)
!    where( ((x(ixC^S,1)**2+x(ixC^S,2)**2).lt.(eqpar(a_)+coordpar(rcut_))**2) &
!         {#IFDEF D3 .and.(abs(x(ixC^S,3)).lt.coordpar(rcut_))} )
!       fE(ixC^S,idir)=zero
!    end where
!  end associate
!  }

end do
!--------------------------------------------

circ(ixI^S,1:ndim)=zero

! Calculate circulation on each face

do idim1=1,ndim ! Coordinate perpendicular to face 
   ixCmax^D=ixOmax^D;
   ixCmin^D=ixOmin^D-kr(idim1,^D);
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
         
        ! Assemble indices
        hxC^L=ixC^L-kr(idim2,^D);
        ! Add line integrals in direction idir
        circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixC^S,idir)&
                         -fE(hxC^S,idir))
      end do
   end do
end do

! Divide by the area of the face to get dB/dt
do idim1=1,ndim
   ixCmax^D=ixOmax^D;
   ixCmin^D=ixOmin^D-kr(idim1,^D);
   select case(idim1)
   {case(^D)
   where(mygeo%surfaceC^D(ixC^S) .gt. 1.0d-9*mygeo%dvolume(ixC^S))
      circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                 /mygeo%surfaceC^D(ixC^S)
   elsewhere
      circ(ixC^S,idim1)=zero
   end where
   \}
    end select

   ! Time update
   ! minus!
   bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
end do

end associate
end subroutine updatefacesuct2
!=============================================================================
subroutine updatefacesuct2av(ixI^L,ixO^L,qdt,vbarC,cbarmin,cbarmax,fC,fE,s)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: qdt
type(state)                        :: s
double precision, intent(in)       :: vbarC(ixI^S,1:ndir,2)
double precision, intent(in)       :: cbarmin(ixI^S,ndim)
double precision, intent(in)       :: cbarmax(ixI^S,ndim)
double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
double precision, intent(inout)    :: fE(ixI^S,1:ndir)
! --- local ---
double precision                   :: vtilL(ixI^S,2)
double precision                   :: vtilR(ixI^S,2)
double precision                   :: btilL(s%ws%ixG^S,ndim)
double precision                   :: btilR(s%ws%ixG^S,ndim)
double precision                   :: sqrtg(s%ws%ixG^S)
double precision                   :: sqrt_gamma_hat(s%ws%ixG^S)
double precision                   :: xC(s%ws%ixG^S,1:ndim),xCC(s%ws%ixG^S,1:ndim)
double precision                   :: cp(ixI^S,2)
double precision                   :: cm(ixI^S,2)
double precision                   :: fEc(ixI^S)
double precision                   :: fEvar(ixI^S)
double precision                   :: fEmin(ixI^S)
double precision                   :: fEmax(ixI^S)
integer                            :: ixGs^L
integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,jxCp^L,ixCm^L
integer                            :: idim1,idim2,idir,iwdim1,iwdim2, idim
double precision                   :: circ(ixI^S,1:ndim), dxidir
integer                            :: i^D
!-----------------------------------------------------------------------------

associate(bfaces=>s%ws%w,w=>s%w%w,x=>s%x%x,mygeo=>s%geo)
{#IFNDEF DY_SP
call mpistop('DY_SP cannot use uct2av')
}
! Calculate contribution to FEM of each edge,
! that is, estimate value of line integral of
! electric field in the positive idir direction.

ixGs^L=s%ws%ixG^L;

! extend one layer of cell center locations in xCC
xCC=0.d0
xCC(ixI^S,1:ndim)=x(ixI^S,1:ndim)
{
xCC(ixGsmin^D^D%ixI^S,1:ndim) = x(ixImin^D^D%ixI^S,1:ndim)
xCC(ixGsmin^D^D%ixGs^S,^D) = x({ixImin^DD,},^D) - s%geo%dx({ixImin^DD,},^D)
\}
{^IFTHREED
xCC(ixImin1:ixImax1,ixGsmin2,ixGsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,1)
xCC(ixGsmin1,ixImin2:ixImax2,ixGsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,2)
xCC(ixGsmin1,ixGsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,3)
}


! Loop over components of electric field

! idir: electric field component we need to calculate
! idim1: directions in which we already performed the reconstruction
! idim2: directions in which we perform the reconstruction

fE(ixI^S,1:ndir)=zero

do idir=7-2*ndim,ndir
  ! Indices
  ! idir: electric field component
  ! idim1: one surface
  ! idim2: the other surface
  ! cyclic permutation: idim1,idim2,idir=1,2,3
  ! Velocity components on the surface
  ! follow cyclic premutations:
  ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

  ixCmax^D=ixOmax^D;
  ixCmin^D=ixOmin^D-1+kr(idir,^D);

  ! Set indices and directions
  idim1=mod(idir,ndir)+1
  idim2=mod(idir+1,ndir)+1

  jxC^L=ixC^L+kr(idim1,^D);
  ixCp^L=ixC^L+kr(idim2,^D);
  jxCp^L=jxC^L+kr(idim2,^D);

  ! Interpolate sqrt gamma to the edges

{#IFDEF DY_SP
  sqrtg(ixC^S)=zero
  ! Get edge coordinates and sqrt_gamma_hat
  do idim=1,ndim
    if (idim/=idir) then
      xC(ixC^S,idim)=xCC(ixC^S,idim) + 0.5d0 * s%geo%dx(ixC^S,idim)
    else
      xC(ixC^S,idim)=xCC(ixC^S,idim)
    end if
  end do
  call get_sqrt_gamma_hat(xC(ixGs^S, 1:ndim), ixGs^L, ixC^L, sqrt_gamma_hat(ixGs^S))
}

{#IFNDEF DY_SP
  sqrtg(ixC^S)=zero
  select case(idim1)
{ case(^D)
    sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
      mygeo%mSurface^D %sqrtgamma(ixC^S)+&
      mygeo%mSurface^D %sqrtgamma(ixCp^S))
\}
  end select
  
  select case(idim2)
{ case(^D)
    sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
      mygeo%mSurface^D %sqrtgamma(ixC^S)+&
      mygeo%mSurface^D %sqrtgamma(jxC^S))
\}
  end select
}

  ! Reconstruct transverse transport velocities
  call rec(ixI^L,ixC^L,idim2,vbarC(ixI^S,idim1,1),&
           vtilL(ixI^S,2),vtilR(ixI^S,2))

  call rec(ixI^L,ixC^L,idim1,vbarC(ixI^S,idim2,2),&
           vtilL(ixI^S,1),vtilR(ixI^S,1))

  ! Reconstruct magnetic fields
  ! Eventhough the arrays are larger, rec works with
  ! the limits ixG.
  call rec(ixI^L,ixC^L,idim2,bfaces(ixI^S,idim1),&
           btilL(ixI^S,idim1),btilR(ixI^S,idim1))

  call rec(ixI^L,ixC^L,idim1,bfaces(ixI^S,idim2),&
           btilL(ixI^S,idim2),btilR(ixI^S,idim2))

  ! Take the maximum characteristic

  cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1))
  cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))

  cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2))
  cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))
 
  ! Calculate elctric field
  if (idir .le. ndim) then
     dxidir = dxlevel(idir)
  else
     dxidir = 1.0d0
  end if

  fE(ixC^S,idir)=qdt * dxidir * &
  {#IFDEF DY_SP
               sqrt_gamma_hat(ixC^S) * (&
  }
  {#IFNDEF DY_SP
               sqrtg(ixC^S) * (&
  }
               -(cp(ixC^S,1)*vtilL(ixC^S,1)*btilL(ixC^S,idim2) &
               + cm(ixC^S,1)*vtilR(ixC^S,1)*btilR(ixC^S,idim2) &
               -cp(ixC^S,1)*cm(ixC^S,1)*(btilR(ixC^S,idim2)-btilL(ixC^S,idim2)))&
               /(cp(ixC^S,1) + cm(ixC^S,1)) &
               +(cp(ixC^S,2)*vtilL(ixC^S,2)*btilL(ixC^S,idim1) &
               + cm(ixC^S,2)*vtilR(ixC^S,2)*btilR(ixC^S,idim1) &
               -cp(ixC^S,2)*cm(ixC^S,2)*(btilR(ixC^S,idim1)-btilL(ixC^S,idim1)))&
               /(cp(ixC^S,2) + cm(ixC^S,2)) &
               )

  
  !if (typeaxial.ne.'slab') then
  if (typeaxial.ne.'slab' .and. typeaxial.ne.'cartesian') then
{^IFZIN
{#IFNDEF HARDBC
  ! Catch axis and origin, where the line integral is performed on a
  ! zero lenght element and should be zero.
  ! Remember that staggered quantities are assigned to the index
  ! of the cell at their left.
  where((abs(x(ixC^S,z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
        (abs(x(ixC^S,z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))
}{#IFDEF HARDBC
  where((abs(x(ixC^S,z_)+half*dxlevel(z_)-xprobmin^Z).lt.smalldouble).or.&
        (abs(x(ixC^S,z_)+half*dxlevel(z_)-xprobmax^Z).lt.smalldouble))
}
    fE(ixC^S,^PHI)=zero
    fE(ixC^S,r_)=zero
  end where}
  where(abs(x(ixC^S,r_)+half*dxlevel(r_)).lt.1.0d-9)
    fE(ixC^S,idir)=zero
  end where
  end if

  ! Allowed range for the electric field in that direction.
  ! For simplicity, use twice cell-centred values, to allow
  ! some overshooting due to reconstruction
  ! Recall: 1,2,3 --> idir, idim1, idim2

  fEc(ixI^S) = myM%beta(idim2)%elem(ixI^S)*w(ixI^S,b0_+idim1) &
             - myM%beta(idim1)%elem(ixI^S)*w(ixI^S,b0_+idim2)

  fEvar(ixI^S) = 2.0*max(abs(w(ixI^S,b0_+idim2) - w(ixI^S,b0_+idim1)),&
                     abs(w(ixI^S,b0_+idim1)),&
                     abs(w(ixI^S,b0_+idim2)),1e-6)

  fEmax(ixI^S) = fEc(ixI^S) + myM%alpha(ixI^S) &
             *fEvar(ixI^S)

  fEmin(ixI^S) = fEc(ixI^S) - myM%alpha(ixI^S) &
             *fEvar(ixI^S)

!  if (any(fEmax(ixI^S).ne.fEmax(ixI^S))) print *,'fEmax is NaN','it=',it
!  if (any(abs(fEmax(ixI^S)).gt.bigdouble)) print *,'fEmax is inf','it=',it
!
!  if (any(fEmin(ixI^S).ne.fEmin(ixI^S))) print *,'fEmin is NaN','it=',it
!  if (any(abs(fEmin(ixI^S)).gt.bigdouble)) print *,'fEmin is inf','it=',it



  ! For the edge-allocated electric fields, take the maximum and
  ! minimum of the surrounding cells

{#IFNDEF D1

  fEmax(ixC^S) = qdt * dxidir * sqrtg(ixC^S) * max(fEmax(ixC^S),fEmax(jxC^S),fEmax(ixCp^S),fEmax(jxCp^S)) 

  fEmin(ixC^S) = qdt * dxidir * sqrtg(ixC^S) * min(fEmin(ixC^S),fEmin(jxC^S),fEmin(ixCp^S),fEmin(jxCp^S)) 

}
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Fall back to Balsara-Spicer averaging
  ! in case of problems
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (any(.not.(fE(ixC^S,idir).lt.fEmax(ixC^S).and.fE(ixC^S,idir).gt.fEmin(ixC^S)))) then
!     write(*,*) mype, 'WARNING: Had to fall back to Balsara-Spicer','it=',it
     iwdim1 = b0_+idim1; iwdim2 = b0_+idim2
     where (.not.(fE(ixC^S,idir).lt.fEmax(ixC^S).and.fE(ixC^S,idir).gt.fEmin(ixC^S)))
        fE(ixC^S,idir)=qdt*quarter*((fC(ixC^S,iwdim1,idim2)&
             +fC(jxC^S,iwdim1,idim2))/dxlevel(idim1)&
             -(fC(ixC^S,iwdim2,idim1)+fC(ixCp^S,iwdim2,idim1))&
             /dxlevel(idim2))
     end where
  end if

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Before clipping the electric field,
  ! forcefully remove NaNs that could be there
  ! Brutalistic (probably never helps
  ! as NaNs will also be in other variables then)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (any(fE(ixC^S,idir) .ne. fE(ixC^S,idir))) then
!     write(*,*) mype, 'WARNING: Had to forecefully set E=0','it=',it
     where (fE(ixC^S,idir) .ne. fE(ixC^S,idir))
        fE(ixC^S,idir) = zero
     end where
  end if
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! If even Balsara-Spicer yields electric
  ! field out of bounds, reset to a fraction
  ! of the allowed field
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (any(.not.(fE(ixC^S,idir).lt.fEmax(ixC^S).and.fE(ixC^S,idir).gt.fEmin(ixC^S)))) then
!     write(*,*) mype, 'WARNING: Clipping the electric field','it=',it
     where (fE(ixC^S,idir).gt.fEmax(ixC^S))
        fE(ixC^S,idir)=fEmin(ixC^S) + 0.8 * (fEmax(ixC^S) - fEmin(ixC^S))
     end where

     where (fE(ixC^S,idir).lt.fEmin(ixC^S))
        fE(ixC^S,idir)=fEmin(ixC^S) + 0.2 * (fEmax(ixC^S) - fEmin(ixC^S))
     end where

  end if

!  if (any(fE(ixC^S,idir).ne.fE(ixC^S,idir))) print *,'fE is NaN','it=',it
!  if (any(abs(fE(ixC^S,idir)).gt.bigdouble)) print *,'fE is inf','it=',it




end do
!--------------------------------------------

circ(ixI^S,1:ndim)=zero

! Calculate circulation on each face

do idim1=1,ndim ! Coordinate perpendicular to face 
   ixCmax^D=ixOmax^D;
   ixCmin^D=ixOmin^D-kr(idim1,^D);
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
         
        ! Assemble indices
        hxC^L=ixC^L-kr(idim2,^D);
        ! Add line integrals in direction idir
        circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixC^S,idir)&
                         -fE(hxC^S,idir))
      end do
   end do
end do

! Divide by the area of the face to get dB/dt
do idim1=1,ndim
   ixCmax^D=ixOmax^D;
   ixCmin^D=ixOmin^D-kr(idim1,^D);
   select case(idim1)
   {case(^D)
   where(mygeo%surfaceC^D(ixC^S) .gt. 1.0d-9*mygeo%dvolume(ixC^S))
      circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                 /mygeo%surfaceC^D(ixC^S)
   elsewhere
      circ(ixC^S,idim1)=zero
   end where
   \}
    end select

   ! Time update
   ! minus!
    bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)

!  if (any(bfaces(ixC^S,idim1) .ne. bfaces(ixC^S,idim1))) &
!       write(*,*) mype, 'WARNING: Magnetic field is NaN in uct2av'

    
end do

end associate
end subroutine updatefacesuct2av
!=============================================================================
subroutine updatefacesuct1(ixI^L,ixO^L,qdt,vbarRC,vbarLC,cbarmin,cbarmax,fE,s)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: qdt
type(state)                        :: s
double precision, intent(in)       :: vbarRC(ixI^S,1:ndir,2)
double precision, intent(in)       :: vbarLC(ixI^S,1:ndir,2)
double precision, intent(in)       :: cbarmin(ixI^S,ndim)
double precision, intent(in)       :: cbarmax(ixI^S,ndim)
double precision, intent(inout)    :: fE(ixI^S,1:ndir)
! --- local ---
double precision                   :: vtilRL(ixI^S,2)
double precision                   :: vtilLL(ixI^S,2)
double precision                   :: vtilRR(ixI^S,2)
double precision                   :: vtilLR(ixI^S,2)
double precision                   :: btilL(s%ws%ixG^S,ndim)
double precision                   :: btilR(s%ws%ixG^S,ndim)
double precision                   :: sqrtg(s%ws%ixG^S)
double precision                   :: sqrt_gamma_hat(s%ws%ixG^S)
double precision                   :: xC(s%ws%ixG^S,1:ndim),xCC(s%ws%ixG^S,1:ndim)
double precision                   :: cp(ixI^S,2)
double precision                   :: cm(ixI^S,2)
double precision                   :: ELL(ixI^S)
double precision                   :: ELR(ixI^S)
double precision                   :: ERL(ixI^S)
double precision                   :: ERR(ixI^S)
integer                            :: ixGs^L
integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
integer                            :: idim1,idim2,idir, idim
double precision                   :: circ(ixI^S,1:ndim), dxidir

!-----------------------------------------------------------------------------

associate(bfaces=>s%ws%w,x=>s%x%x,mygeo=>s%geo)

! Calculate contribution to FEM of each edge,
! that is, estimate value of line integral of
! electric field in the positive idir direction.

ixGs^L=s%ws%ixG^L;

! extend one layer of cell center locations in xCC
xCC=0.d0
xCC(ixI^S,1:ndim)=x(ixI^S,1:ndim)
{
xCC(ixGsmin^D^D%ixI^S,1:ndim) = x(ixImin^D^D%ixI^S,1:ndim)
xCC(ixGsmin^D^D%ixGs^S,^D) = x({ixImin^DD,},^D) - s%geo%dx({ixImin^DD,},^D)
\}
{^IFTHREED
xCC(ixImin1:ixImax1,ixGsmin2,ixGsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,1)
xCC(ixGsmin1,ixImin2:ixImax2,ixGsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,2)
xCC(ixGsmin1,ixGsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,3)
}


! Loop over components of electric field

! idir: electric field component we need to calculate
! idim1: directions in which we already performed the reconstruction
! idim2: directions in which we perform the reconstruction

fE(ixI^S,1:ndir)=zero

do idir=7-2*ndim,ndir
  ! Indices
  ! idir: electric field component
  ! idim1: where the first reconstruction was done
  ! idim2: where we will do the second reconstruction
  ! cyclic permutation: idim1,idim2,idir=1,2,3

  ixCmax^D=ixOmax^D;
  ixCmin^D=ixOmin^D-1+kr(idir,^D);

 ! Set indices and directions
  idim1=mod(idir,ndir)+1
  idim2=mod(idir+1,ndir)+1

  jxC^L=ixC^L+kr(idim1,^D);
  ixCp^L=ixC^L+kr(idim2,^D);

  ! Interpolate sqrt gamma to the edges

{#IFDEF DY_SP
  sqrtg(ixC^S)=zero
  ! Get edge coordinates and sqrt_gamma_hat
  do idim=1,ndim
    if (idim/=idir) then
      xC(ixC^S,idim)=xCC(ixC^S,idim) + 0.5d0 * s%geo%dx(ixC^S,idim)
    else
      xC(ixC^S,idim)=xCC(ixC^S,idim)
    end if
  end do
  call get_sqrt_gamma_hat(xC(ixGs^S, 1:ndim), ixGs^L, ixC^L, sqrt_gamma_hat(ixGs^S))

}

{#IFNDEF DY_SP
  sqrtg(ixC^S)=zero
  select case(idim1)
{ case(^D)
    sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
      mygeo%mSurface^D %sqrtgamma(ixC^S)+&
      mygeo%mSurface^D %sqrtgamma(ixCp^S))
\}
  end select
  
  select case(idim2)
{ case(^D)
    sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
      mygeo%mSurface^D %sqrtgamma(ixC^S)+&
      mygeo%mSurface^D %sqrtgamma(jxC^S))
\}
  end select
}
! any DY_SP or Non DY_SP, Nan in uct1 sometimes
!write(*,*) sqrtg(ixC^S)
!stop 'uct1'

  ! Reconstruct velocities
  call rec(ixI^L,ixC^L,idim2,vbarLC(ixI^S,idir,1),&
           vtilLL(ixI^S,1),vtilLR(ixI^S,1)) 
  call rec(ixI^L,ixC^L,idim2,vbarRC(ixI^S,idir,1),&
           vtilRL(ixI^S,1),vtilRR(ixI^S,1))

  call rec(ixI^L,ixC^L,idim2,vbarLC(ixI^S,idir,2),&
           vtilLL(ixI^S,2),vtilLR(ixI^S,2)) 
  call rec(ixI^L,ixC^L,idim2,vbarRC(ixI^S,idir,2),&
           vtilRL(ixI^S,2),vtilRR(ixI^S,2))

  ! Reconstruct magnetic fields
  ! Eventhough the arrays are larger, rec works with
  ! the limits ixG.
  call rec(ixI^L,ixC^L,idim2,bfaces(ixI^S,idim1),&
           btilL(ixI^S,idim1),btilR(ixI^S,idim1))
  call rec(ixI^L,ixC^L,idim1,bfaces(ixI^S,idim2),&
           btilL(ixI^S,idim2),btilR(ixI^S,idim2))

  ! Take the maximum characteristic
  cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1)) 
  cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))

  cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2)) 
  cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))

  ! Calculate component idir of vxB partial
  ELL(ixC^S)= vtilLL(ixC^S,2)*btilL(ixC^S,idim1)&
       -vtilLL(ixC^S,1)*btilL(ixC^S,idim2)

  ELR(ixC^S)= vtilLR(ixC^S,2)*btilR(ixC^S,idim1)&
       -vtilLR(ixC^S,1)*btilL(ixC^S,idim2)

  ERL(ixC^S)= vtilRL(ixC^S,2)*btilL(ixC^S,idim1)&
       -vtilRL(ixC^S,1)*btilR(ixC^S,idim2)

  ERR(ixC^S)= vtilRR(ixC^S,2)*btilR(ixC^S,idim1)&
       -vtilRR(ixC^S,1)*btilR(ixC^S,idim2)

  ! For 3D, use interval
  if (idir .le. ndim) then
     dxidir = dxlevel(idir)
  else
     dxidir = 1.0d0
  end if
  
  ! Calculate elctric field
  fE(ixC^S,idir)=qdt * dxidir * &
  {#IFDEF DY_SP
   sqrt_gamma_hat(ixC^S) * (&
  }
  {#IFNDEF DY_SP
   sqrtg(ixC^S) * (&
  }
   (cp(ixC^S,1)*cp(ixC^S,2)*ELL(ixC^S)&
   +cp(ixC^S,1)*cm(ixC^S,2)*ELR(ixC^S)&
   +cm(ixC^S,1)*cp(ixC^S,2)*ERL(ixC^S)&
   +cm(ixC^S,1)*cm(ixC^S,2)*ERR(ixC^S)&
   )/&
    ((cp(ixC^S,1)+cm(ixC^S,1))*&
     (cp(ixC^S,2)+cm(ixC^S,2)))&
   +(cp(ixC^S,1)*cm(ixC^S,1)*(btilR(ixC^S,idim2)-btilL(ixC^S,idim2)))/&
     (cp(ixC^S,1) + cm(ixC^S,1))&
   -(cp(ixC^S,2)*cm(ixC^S,2)*(btilR(ixC^S,idim1)-btilL(ixC^S,idim1)))/&
     (cp(ixC^S,2) + cm(ixC^S,2)))

  !if (typeaxial.ne.'slab') then
  if (typeaxial.ne.'slab' .and. typeaxial.ne.'cartesian') then
{^IFZIN
  ! Catch axis and origin, where the line integral is performed on a
  ! zero lenght element and should be zero.
  ! Remember that staggered quantities are assigned to the index
  ! of the cell at their left.
  where((abs(x(ixC^S,z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
        (abs(x(ixC^S,z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))
    fE(ixC^S,^PHI)=zero
    fE(ixC^S,r_)=zero
  end where}
  where(abs(x(ixC^S,r_)+half*dxlevel(r_)).lt.1.0d-9)
    fE(ixC^S,idir)=zero
  end where
  end if

end do


!--------------------------------------------

circ(ixI^S,1:ndim)=zero

! Calculate circulation on each face

do idim1=1,ndim ! Coordinate perpendicular to face 
   ixCmax^D=ixOmax^D;
   ixCmin^D=ixOmin^D-kr(idim1,^D);
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
                  
        ! Assemble indices
        hxC^L=ixC^L-kr(idim2,^D);
        ! Add line integrals in direction idir
        circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixC^S,idir)&
                         -fE(hxC^S,idir))
      end do
   end do
end do

! Divide by the area of the face to get dB/dt
do idim1=1,ndim
   ixCmax^D=ixOmax^D;
   ixCmin^D=ixOmin^D-kr(idim1,^D);
   select case(idim1)
   {case(^D)
   where(mygeo%surfaceC^D(ixC^S) .gt. 1.0d-9*mygeo%dvolume(ixC^S))
      circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                 /mygeo%surfaceC^D(ixC^S)
   elsewhere
      circ(ixC^S,idim1)=zero
   end where
   \}
    end select

   ! Time update
   ! minus!
   bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
end do

end associate
end subroutine updatefacesuct1
!=============================================================================
subroutine storeedge(igrid,ixI^L,fE,idim^LIM)
!use mod_fix_conserve
include 'amrvacdef.f'

integer, intent(in)          :: igrid, ixI^L, idim^LIM
double precision, intent(in) :: fE(ixI^S,1:^NC)

integer :: idims, idir, iside, i^D
integer :: pi^D, mi^D, ph^D, mh^D ! To detect corners
integer :: ixMc^L
!-----------------------------------------------------------------------------

do idims = idim^LIM  !loop over face directions
  !! Loop over block faces
  do iside=1,2 
    i^D=kr(^D,idims)*(2*iside-3);
    {^IFPHI if (neighbor_pole(i^D,igrid)/=0) cycle}
    select case (neighbor_type(i^D,igrid))
    case (4)
      ! The neighbour is finer
      ! Face direction, side (left or right), restrict required?, fE
      call fluxtoedge(igrid,ixI^L,idims,iside,.false.,fE)
    case(2)
      ! The neighbour is coarser
      call fluxtoedge(igrid,ixI^L,idims,iside,.true.,fE)
    case(3)
      ! If the neighbour is at the same level,
      ! check if there are corners
      ! If there is any corner, store the fluxes from that side
      do idir=idims+1,ndim
        pi^D=i^D+kr(idir,^D);
        mi^D=i^D-kr(idir,^D);
        ph^D=pi^D-kr(idims,^D)*(2*iside-3);
        mh^D=mi^D-kr(idims,^D)*(2*iside-3);
        if (neighbor_type(pi^D,igrid).eq.4&
          .and.neighbor_type(ph^D,igrid).eq.3) then
          call fluxtoedge(igrid,ixI^L,idims,iside,.false.,fE)
        end if
        if (neighbor_type(mi^D,igrid).eq.4&
          .and.neighbor_type(mh^D,igrid).eq.3) then
          call fluxtoedge(igrid,ixI^L,idims,iside,.false.,fE)
        end if
      end do
    end select
  end do
end do

end subroutine storeedge
!=============================================================================
subroutine fluxtoedge(igrid,ixI^L,idims,iside,restrict,fE)
use mod_fix_conserve
include 'amrvacdef.f'

integer                      :: igrid,ixI^L,idims,iside
logical                      :: restrict
double precision, intent(in) :: fE(ixI^S,1:^NC)
! ... local ...
integer                      :: idir1,idir2
integer                      :: ixE^L,ixF^L{^IFTHREED, jxF^L,}, nx^D,nxCo^D
!-----------------------------------------------------------------------------

nx^D=ixMhi^D-ixMlo^D+1;
nxCo^D=nx^D/2;

! ixE are the indices on the 'edge' array.
! ixF are the indices on the 'fE' array
! jxF are indices advanced to perform the flux restriction (sum) in 3D
! A line integral of the electric field on the coarse side
! lies over two edges on the fine side. So, in 3D we restrict by summing
! over two cells on the fine side.

do idir1=1,ndim-1
 {^IFTHREED ! 3D: rotate indices among 1 and 2 to save space 
  idir2=mod(idir1+idims-1,3)+1}
 {^IFTWOD ! Assign the only field component (3) to the only index (1)
  idir2=3}

  if (restrict) then
    ! Set up indices for restriction
    ixFmin^D=ixMlo^D-1+kr(^D,idir2);
    ixFmax^D=ixMhi^D-kr(^D,idir2);
    {^IFTHREED
    jxF^L=ixF^L+kr(^D,idir2);}

    ixEmin^D=0+kr(^D,idir2);
    ixEmax^D=nxCo^D;
    select case(idims)
   {case(^D)
      ixEmin^D=1;ixEmax^D=1;
      select case(iside)
      case(1)
        ixFmax^D=ixFmin^D
        {^IFTHREEDjxFmax^D=ixFmin^D}
      case(2)
        ixFmin^D=ixFmax^D
        {^IFTHREEDjxFmin^D=ixFmax^D}
      end select
   \}
    end select

  pflux(iside,idims,igrid)%edge(ixE^S,idir1)=&
    fE(ixFmin^D:ixFmax^D:2,idir2){^IFTHREED +&
    fE(jxFmin^D:jxFmax^D:2,idir2)};

  else
    ! Set up indices for copying 
    ixFmin^D=ixMlo^D-1+kr(^D,idir2);
    ixFmax^D=ixMhi^D;

    ixEmin^D=0+kr(^D,idir2);
    ixEmax^D=nx^D;

    select case(idims)
   {case(^D)
      ixEmin^D=1;ixEmax^D=1;
      select case(iside)
      case(1)
        ixFmax^D=ixFmin^D
      case(2)
        ixFmin^D=ixFmax^D
      end select
   \}
    end select

    pflux(iside,idims,igrid)%edge(ixE^S,idir1)=&
    fE(ixF^S,idir2)

  end if

end do

end subroutine fluxtoedge
!=============================================================================
subroutine fix_edges(psuse,idim^LIM)
use mod_fix_conserve
include 'amrvacdef.f'

type(state) :: psuse(ngridshi)
integer, intent(in) :: idim^LIM

integer :: iigrid, igrid, idims, iside, iotherside, i^D, ic^D, inc^D, ixMc^L
integer :: nbuf, ibufnext
integer :: ibufnext_cc
integer :: pi^D, mi^D, ph^D, mh^D ! To detect corners
integer :: ixE^L(1:ndir), ixtE^L, ixF^L(1:ndim), ixfE^L(1:ndir)
integer :: nx^D, idir, ix, ipe_neighbor, ineighbor
logical :: pcorner(1:ndim),mcorner(1:ndim)
!-----------------------------------------------------------------------------
if (nrecv_ct>0) then
   call MPI_WAITALL(nrecv_ct,recvrequest_stg,recvstatus_stg,ierrmpi)
end if

! Initialise buffer counter again
ibuf=1
ibuf_cc=1

do iigrid=1,igridstail; igrid=igrids(iigrid);
  do idims= idim^LIM
    do iside=1,2
      i^D=kr(^D,idims)*(2*iside-3);
      {^IFPHI if (neighbor_pole(i^D,igrid)/=0) cycle}
      select case(neighbor_type(i^D,igrid))
      case(4)
      ! The first neighbour is finer
      if (.not.neighbor_active(i^D,igrid).or.&
          .not.neighbor_active(0^D&,igrid) ) then
         {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
         inc^DB=2*i^DB+ic^DB\}
         ipe_neighbor=neighbor_child(2,inc^D,igrid)
         !! When the neighbour is in a different process
         if (ipe_neighbor/=mype) then
            ibufnext=ibuf+isize(idims)
            ibuf=ibufnext
            end if
         {end do\}
         cycle
      end if

      ! Check if there are corners
      pcorner=.false.
      mcorner=.false.
      do idir=1,ndim
        pi^D=i^D+kr(idir,^D);
        mi^D=i^D-kr(idir,^D);
        ph^D=pi^D-kr(idims,^D)*(2*iside-3);
        mh^D=mi^D-kr(idims,^D)*(2*iside-3);
        if (neighbor_type(ph^D,igrid).eq.4) pcorner(idir)=.true.
        if (neighbor_type(mh^D,igrid).eq.4) mcorner(idir)=.true.
      end do

      ! Calculate indices range
      call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.false.,.false.,0^D&,pcorner,mcorner)

      ! Remove coarse part of circulation
      call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iside,idims,igrid)%edge,idims,iside,.false.,psuse(igrid))

      ! Add fine part of the circulation
      {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
         inc^DB=2*i^DB+ic^DB\}
         ineighbor=neighbor_child(1,inc^D,igrid)
         ipe_neighbor=neighbor_child(2,inc^D,igrid)
         iotherside=3-iside
         nx^D=(ixMhi^D-ixMlo^D+1)/2;

         call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.true.,.false.,inc^D,pcorner,mcorner)
         if (ipe_neighbor.eq.mype) then
         call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iotherside,idims,ineighbor)%edge,idims,iside,.true.,psuse(igrid))

         else
         ibufnext=ibuf+isize(idims)
         call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,&
         reshape(source=recvbuffer(ibufnext-isize_stg(idims):ibufnext-1),&
         shape=(/ ixtEmax^D-ixtEmin^D+1 ,^ND-1 /)),&
         idims,iside,.true.,psuse(igrid))

         ibuf=ibufnext

         end if
      {end do\}

      case(3)
      ! The first neighbour is at the same level
      ! Check if there are corners
      do idir=idims+1,ndim
        pcorner=.false.
        mcorner=.false.
        pi^D=i^D+kr(idir,^D);
        mi^D=i^D-kr(idir,^D);
        ph^D=pi^D-kr(idims,^D)*(2*iside-3);
        mh^D=mi^D-kr(idims,^D)*(2*iside-3);
        if (neighbor_type(pi^D,igrid).eq.4&
          .and.neighbor_type(ph^D,igrid).eq.3&
          .and.neighbor_pole(pi^D,igrid).eq.0) then
          pcorner(idir)=.true.
        ! Remove coarse part
        ! Set indices
          call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.false.,.true.,0^D&,pcorner,mcorner)
        ! Remove
          call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iside,idims,igrid)%edge,idims,iside,.false.,psuse(igrid))

        ! Add fine part
        ! Find relative position of finer grid
         {^IFTHREED do ix=1,2}

          inc^D=kr(idims,^D)*3*(iside-1)+3*kr(idir,^D){^IFTHREED+kr(6-idir-idims,^D)*ix};
          ineighbor=neighbor_child(1,inc^D,igrid)
          ipe_neighbor=neighbor_child(2,inc^D,igrid)
          iotherside=3-iside

        ! Set indices
          call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.true.,.true.,inc^D,pcorner,mcorner)

        ! add

          if (ipe_neighbor.eq.mype) then

          call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iotherside,idims,ineighbor)%edge,idims,iside,.true.,psuse(igrid))

          else


          ibufnext_cc=ibuf_cc+isize_stg(idims)
          call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,&
          reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
          shape=(/ ixtEmax^D-ixtEmin^D+1 ,^ND-1 /)),&
          idims,iside,.true.,psuse(igrid))

          ibuf_cc=ibufnext_cc

          end if

         {^IFTHREED end do}
        ! Set CoCorner to false again for next step

          pcorner(idir)=.false.

          !call div_staggered(ixI^L,psuse(igrid))
        end if

        if (neighbor_type(mi^D,igrid).eq.4&
          .and.neighbor_type(mh^D,igrid).eq.3&
          .and.neighbor_pole(mi^D,igrid).eq.0) then
          mcorner(idir)=.true.
        ! Remove coarse part
        ! Set indices
          call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.false.,.true.,0^D&,pcorner,mcorner)
        ! Remove
          call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iside,idims,igrid)%edge,idims,iside,.false.,psuse(igrid))

        ! Add fine part
        ! Find relative position of finer grid
         {^IFTHREED do ix=1,2}

          inc^D=kr(idims,^D)*3*(iside-1){^IFTHREED+kr(6-idir-idims,^D)*ix};
          ineighbor=neighbor_child(1,inc^D,igrid)
          ipe_neighbor=neighbor_child(2,inc^D,igrid)
          iotherside=3-iside

        ! Set indices
          call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.true.,.true.,inc^D,pcorner,mcorner)

        ! add

          if (ipe_neighbor.eq.mype) then

          call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iotherside,idims,ineighbor)%edge,idims,iside,.true.,psuse(igrid))

          else

          ibufnext_cc=ibuf_cc+isize_stg(idims)
          call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,&
          reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
          shape=(/ ixtEmax^D-ixtEmin^D+1 ,^ND-1 /)),&
          idims,iside,.true.,psuse(igrid))

          ibuf_cc=ibufnext_cc

          end if

         {^IFTHREED end do}
        ! Set CoCorner to false again for next step

         mcorner(idir)=.false.
        end if
      end do
      end select
    end do
  end do
 ! Average to centers
 call faces2centers(ixG^LL,psuse(igrid)) 
end do

if (nrecv>0) deallocate(recvbuffer,recvstatus,recvrequest)
if (nrecv_ct>0) deallocate(recvbuffer_cc,recvstatus_stg,recvrequest_stg)

if (nsend>0) then
   call MPI_WAITALL(nsend,sendrequest,sendstatus,ierrmpi)
   deallocate(sendbuffer,sendstatus,sendrequest)
end if

if (nsend_ct>0) then
   call MPI_WAITALL(nsend_ct,sendrequest_stg,sendstatus_stg,ierrmpi)
   deallocate(sendbuffer_cc,sendstatus_stg,sendrequest_stg)
end if

end subroutine fix_edges
!=============================================================================
subroutine set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,add,CoCorner,inc^D,pcorner,mcorner)
include 'amrvacdef.f'
! This routine sets the indexes for the correction
! of the circulation according to several different
! cases, as grids located in different cpus,
! presence of corners, and different relative locations
! of the fine grid respect to the coarse one

integer,intent(in)    :: igrid,idims,iside,inc^D
logical,intent(in)    :: add,CoCorner
logical,intent(inout) :: pcorner(1:ndim),mcorner(1:ndim)
integer,intent(out)   :: ixF^L(1:ndim),ixtE^L,ixE^L(1:ndir),ixfE^L(1:ndir) ! Indices for faces and edges
! ... local ...
integer               :: icor^D,idim1,idir,nx^D,middle^D
integer               :: ixtfE^L
!-----------------------------------------------------------------------------
! ixF -> Indices for the _F_aces, and
! depends on the field component
! ixtE -> are the _t_otal range of the 'edge' array
! ixE -> are the ranges of the edge array,
! depending on the component
! ixfE -> are the ranges of the fE array (3D),
! and also depend on the component

! ... General ...
! Assign indices for the size of the E field array

ixtfEmin^D=ixMlo^D-1;
ixtfEmax^D=ixMhi^D;

if (add) then
nx^D=(ixMhi^D-ixMlo^D+1)/2;
else
nx^D=ixMhi^D-ixMlo^D+1;
end if

do idim1=1,ndim
  ixtEmin^D=0;
  ixtEmax^D=nx^D;
  select case(idims)
  {case(^D)
    ixtEmin^D=1;ixtEmax^D=1;
    if (iside.eq.1) ixtfEmax^D=ixtfEmin^D;
    if (iside.eq.2) ixtfEmin^D=ixtfEmax^D;
  \}
  end select
end do

! Assign indices, considering only the face
! (idims and iside)
do idim1=1,ndim
  ixFmin^D(idim1)=ixMlo^D-kr(idim1,^D);
  ixFmax^D(idim1)=ixMhi^D;
  select case(idims)
  {case(^D)
     select case(iside)
     case(1)
     ixFmax^D(idim1)=ixFmin^D(idim1)
     case(2)
     ixFmin^D(idim1)=ixFmax^D(idim1)
     end select
  \}
  end select
end do

! ... Relative position ...
! Restrict range using relative position
if (add) then
middle^D=(ixMhi^D+ixMlo^D)/2;
{
if (inc^D.eq.1) then
ixFmax^D(:)=middle^D
ixtfEmax^D=middle^D
end if
if (inc^D.eq.2) then
ixFmin^D(:)=middle^D+1
ixtfEmin^D=middle^D
end if
\}
end if

! ... Adjust ranges of edges according to direction ...

do idim1=1,ndir
  ixfEmax^D(idim1)=ixtfEmax^D;
  ixEmax^D(idim1)=ixtEmax^D;
  ixfEmin^D(idim1)=ixtfEmin^D+kr(idim1,^D);
  ixEmin^D(idim1)=ixtEmin^D+kr(idim1,^D);
end do



! ... Corners ...
! 'Coarse' corners
if (CoCorner) then
  do idim1=idims+1,ndim
    if (pcorner(idim1)) then
      do idir=1,ndir !Index arrays have size ndim
        if (idir.eq.6-idim1-idims) then
         !!! Something here has to change
         !!! Array ixfE must have size ndir, while
         !!! ixE must have size ndim
         {if (^D.eq.idim1) then
            ixfEmin^D(idir)=ixfEmax^D(idir)
            if (add) then
              ixEmax^D(idir) =ixEmin^D(idir)
            else
              ixEmin^D(idir) =ixEmax^D(idir)
            end if
          end if\}
        else
          ixEmin^D(idir)=1;
          ixEmax^D(idir)=0;
          ixfEmin^D(idir)=1;
          ixfEmax^D(idir)=0;
        end if
      end do
    end if
    if (mcorner(idim1)) then
      do idir=1,ndir
        if (idir.eq.6-idim1-idims) then
         {if (^D.eq.idim1) then
            ixfEmax^D(idir)=ixfEmin^D(idir)
            if (add) then
              ixEmin^D(idir) =ixEmax^D(idir)
            else
              ixEmax^D(idir) =ixEmin^D(idir)
            end if
          end if\}
        else
          ixEmin^D(idir)=1;
          ixEmax^D(idir)=0;
          ixfEmin^D(idir)=1;
          ixfEmax^D(idir)=0;
        end if
      end do
    end if
  end do
else
! Other kinds of corners
!! Crop ranges to account for corners
!! When the fine fluxes are added, we consider 
!! whether they come from the same cpu or from
!! a different one, in order to mimimise the 
!! amount of communication
!
!!!!Case for different processors still not implemented!!!
  
   {if((idims.gt.^D).and.pcorner(^D)) then
      if((.not.add).or.(inc^D.eq.2)) then
        !ixFmax^DD(:)=ixFmax^DD(:)-kr(^D,^DD);
        do idir=1,ndir
          if ((idir.eq.idims).or.(idir.eq.^D)) cycle
            ixfEmax^D(idir)=ixfEmax^D(idir)-1
            ixEmax^D(idir)=ixEmax^D(idir)-1
        end do
      end if
    end if\}
   {if((idims.gt.^D).and.mcorner(^D)) then
      if((.not.add).or.(inc^D.eq.1)) then
        !ixFmin^DD(:)=ixFmin^DD(:)+kr(^D,^DD);
        do idir=1,ndir
          if ((idir.eq.idims).or.(idir.eq.^D)) cycle
            ixfEmin^D(idir)=ixfEmin^D(idir)+1
            ixEmin^D(idir)=ixEmin^D(idir)+1
        end do
      end if
    end if\}
end if

end subroutine set_ix_circ
!=============================================================================
subroutine add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,edge,idims,iside,add,s)
use mod_fix_conserve
include 'amrvacdef.f'

type(state)        :: s
integer,intent(in) :: idims,iside
integer            :: ixF^L(1:ndim),ixtE^L,ixE^L(1:ndir),ixfE^L(1:ndir)
double precision   :: edge(ixtE^S,1:ndim-1)
logical,intent(in) :: add
! ... local ...
integer            :: idim1,idim2,idir,middle^D
integer            :: ixfEC^L,ixEC^L
double precision   :: fE(ixG^T,1:^NC) !!!!!!!!
double precision   :: circ(ixG^T,1:ndim) !!!!!!!!
integer            :: ix^L,hx^L,ixC^L,hxC^L ! Indices for edges
!-----------------------------------------------------------------------------
! ixF -> Indices for the faces, depends on the field component
! ixE -> Total range for the edges
! ixfE -> Edges in fE (3D) array
associate(mygeo=>s%geo,bfaces=>s%ws%w)
! ix,hx,ixC,hxC -> Auxiliary indices

! Assign quantities stored ad edges to make it as similar as 
! possible to the routine updatefaces.

fE(:^D&,:)=zero

do idim1=1,ndim-1
  {^IFTHREED ! 3D: rotate indices (see routine fluxtoedge)
  idir=mod(idim1+idims-1,3)+1}
  {^IFTWOD   ! 2D: move E back to directon 3
  idir=3}
  ixfEC^L=ixfE^L(idir);
  ixEC^L=ixE^L(idir);
  fE(ixfEC^S,idir)=edge(ixEC^S,idim1)

end do

! Calculate part of circulation needed
circ=zero
do idim1=1,ndim
   do idim2=1,ndim
      do idir=1,ndir
        if (lvc(idim1,idim2,idir).eq.0) cycle
        ! Assemble indices
        ixC^L=ixF^L(idim1);
        hxC^L=ixC^L-kr(idim2,^D);
        if (idim1.eq.idims) then
        circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixC^S,idir)&
                         -fE(hxC^S,idir))

        else
         
          select case(iside)
          case(2)
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *fE(ixC^S,idir)

          case(1)
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                         -lvc(idim1,idim2,idir)&
                         *fE(hxC^S,idir)

          end select
        end if
      end do
   end do
end do

! Divide circulation by surface and add
do idim1=1,ndim
   ixC^L=ixF^L(idim1);
   select case(idim1)
   {case(^D)
   where(mygeo%surfaceC^D(ixC^S) .gt. 1.0d-9*mygeo%dvolume(ixC^S))
      circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                 /mygeo%surfaceC^D(ixC^S)
   elsewhere
      circ(ixC^S,idim1)=zero
   end where
   \}
    end select
   ! Add/subtract to field at face

   if (add) then
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
   else
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)+circ(ixC^S,idim1)
   end if

end do

end associate

end subroutine add_sub_circ
}
!=============================================================================
{#IFNDEF STAGGERED
subroutine b_from_vectorpotential(ixI^L, ixO^L, w, x)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)

integer                            :: ixC^L, ixCp^L, ixCm^L, ixOm^L, hxO^L, idim, idir, j, k
double precision                   :: xC(ixI^S,1:ndim), A(ixI^S,1:ndir), dxC(ixI^S,1:ndim)
double precision                   :: B(ixI^S,1:^ND)
!-----------------------------------------------------------------------------

A(:^D&,:)=zero

{ixCmax^D=ixOmax^D;}
{ixCmin^D=ixOmin^D-1;} ! Extend range by one

do idir=7-2*ndim,ndir
 do idim=1,ndim
  ! Get corner coordinates
  if (idim.ne.idir) then
   {ixCp^L=ixC^L+kr(idim,^D);}
   {ixCm^L=ixC^L-kr(idim,^D);}
   xC(ixC^S,idim) = half * (x(ixCp^S,idim) + x(ixC^S,idim))
   dxC(ixC^S,idim) = xC(ixC^S,idim) - xC(ixCm^S,idim)
  else
   xC(ixC^S,idim)=x(ixC^S,idim)
  end if
 end do
 ! Initialise vector potencial at the corners

 call initvecpot_usr(ixI^L, ixC^L, xC, A(ixI^S,idir), idir)

end do

! ---
! take the curl of the vectorpotential: 

B(:^D&,:) = zero
do idir = 1, ^ND
   do j = 1, ndim
      if (idir .eq. j) cycle
      {ixCmin^D=ixOmin^D-kr(idir,^D);}
      {ixCm^L=ixC^L-kr(j,^D);}
      do k = 1,ndir
         ! Field on the faces
         select case(idir)
         {case(^D)
         {^IFTWOD
         where (mygeo%surfaceC^D(ixC^S) .ne. zero)
            B(ixC^S,idir) = B(ixC^S,idir) &
                 + lvc(idir,j,k) &
                 * (A(ixC^S,k)-A(ixCm^S,k))
         elsewhere
            B(ixC^S,idir) = zero
         end where
         }
         {^IFTHREED
         where (mygeo%surfaceC^D(ixC^S) .ne. zero)
            B(ixC^S,idir) = B(ixC^S,idir) &
                 + lvc(idir,j,k) * dxlevel(k) &
                 * (A(ixC^S,k)-A(ixCm^S,k)) 
         elsewhere
            B(ixC^S,idir) = zero
         end where
         }
         \}
         end select
      end do
   end do
end do



! Average to the cell centers and fill solution array:

do idir = 1, ^ND
   {ixOm^L=ixO^L-kr(idir,^D);}
   w(ixO^S,b0_+idir) = half * dxC(ixO^S,idir)*(B(ixO^S,idir) + B(ixOm^S,idir))/mygeo%dvolume(ixO^S)
end do

end subroutine b_from_vectorpotential
!=============================================================================
}
{#IFDEF STAGGERED
subroutine b_from_vectorpotential(ixIs^L, ixI^L, ixO^L, ws, x)

include 'amrvacdef.f'

integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
double precision, intent(inout)    :: ws(ixIs^S,1:nws)
double precision, intent(in)       :: x(ixI^S,1:ndim)
! ... local ...
double precision                   :: Adummy(ixI^S,1:ndir)
!----------------------------------------------------------------

call b_from_vectorpotentialA(ixIs^L, ixI^L, ixO^L, ws, x, Adummy)

end subroutine b_from_vectorpotential
!=============================================================================
subroutine b_from_vectorpotentialA(ixIs^L, ixI^L, ixO^L, ws, x, A)

include 'amrvacdef.f'

integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
double precision, intent(inout)    :: ws(ixIs^S,1:nws),A(ixI^S,1:ndir)
double precision, intent(in)       :: x(ixI^S,1:ndim)
! .. local ..
integer                            :: ixC^L, hxC^L, ixCp^L, ixCm^L, hxO^L, idim, idim1, idim2, idir
double precision                   :: xC(ixI^S,1:ndim)
double precision                   :: circ(ixI^S,1:ndim), dxidir
!-----------------------------------------------------------------------------

A(:^D&,:)=zero
ws(:^D&,:)=zero

{ixCmax^D=ixOmax^D;}
{ixCmin^D=ixOmin^D-1;} ! Extend range by one
do idir=7-2*ndim,ndir
 do idim=1,ndim
  ! Get edge coordinates
  if (idim.ne.idir) then
   {ixCp^L=ixC^L+kr(idim,^D);}
   {ixCm^L=ixC^L-kr(idim,^D);}
   xC(ixC^S,idim) = half * (x(ixCp^S,idim) + x(ixC^S,idim))
  else
   xC(ixC^S,idim)=x(ixC^S,idim)
  end if
end do

 ! Initialise vector potencial at the edge
 call initvecpot_usr(ixI^L, ixC^L, xC, A(ixI^S,idir), idir)


end do

! Set NaN to zero (can happen e.g. on axis):
where(A(ixG^T,1:ndir).ne.A(ixG^T,1:ndir))
   A(ixG^T,1:ndir)=zero
end where


! --------------------------------------------------
! Take the curl of the vector potential 
! --------------------------------------------------
circ(:^D&,:) = zero

! Calculate circulation on each face
do idim1=1,ndim ! Coordinate perpendicular to face 
   ixCmax^D=ixOmax^D;
   ixCmin^D=ixOmin^D-kr(idim1,^D);
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
        
        ! Assemble indices
        hxC^L=ixC^L-kr(idim2,^D);
        ! Add line integrals in direction idir
        if (idir .le. ndim) then
                dxidir = dxlevel(idir)
        else
                dxidir = 1.0d0
        end if
        circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                         +lvc(idim1,idim2,idir)*dxidir &
                         *(A(ixC^S,idir)&
                         -A(hxC^S,idir))
      end do
   end do
end do

! Set NaN to zero (should not occur)
where(circ.ne.circ)
   circ=zero
end where


! --------------------------------------------------
! Divide by the area of the face to get B
! --------------------------------------------------
do idim1=1,ndim
   ixCmax^D=ixOmax^D;
   ixCmin^D=ixOmin^D-kr(idim1,^D);
   select case(idim1)
   {case(^D)
     where(mygeo%surfaceC^D(ixC^S) .gt. 1.0d-9*mygeo%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                   /mygeo%surfaceC^D(ixC^S)
     elsewhere
        circ(ixC^S,idim1)=zero
     end where
   \}
   end select
   ws(ixC^S,idim1) = circ(ixC^S,idim1)
end do

end subroutine b_from_vectorpotentialA
!=============================================================================
}
subroutine printarray(ixI^L,ixO^L,array)
! Routine for easily printing an array, just for debugging, erase later

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: array(ixO^S)
!double precision, intent(in)  :: array(:,:)
integer                       :: i,j,k
!----------------------------------------------------------------------------

print *, array

{^IFTWOD
do j=ixOmax2,ixOmin2,-1
}
   do i=ixOmin1,ixOmax1

{^IFONED
      write (*,'(I2,E17.9,$)') i,array(i)
}
{^IFTWOD
      write (*,'(I2,I2,E17.9,$)') i,j,array(i,j)
}
{^IFTHREED
      do k=ixOmin3,ixOmax3
        write (*,'(I2,I2,I2,E17.9,$)') i,j,k,array(i,j,k)
      end do
      print *, ' '
      }
   end do

   print *, '---'
{^IFTWOD
end do
}

write (*,'(A)') '------------------------------------------------------------'
end subroutine printarray
!=============================================================================
