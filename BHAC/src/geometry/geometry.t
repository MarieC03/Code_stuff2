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
subroutine set_pole

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
select case (typeaxial)
case ("spherical") {^IFTHREED
  ! For spherical grid, check whether phi-direction is periodic
  if(periodB(ndim)) then
    if(^PHI/=3) call mpistop("set setamrvac -phi=3 and recompile")
    if(mod(ng3(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin2)<smalldouble) then
      if(mype==0) write(unitterm,*) &
        "Will apply pi-periodic conditions at northpole!"
      poleB(1,2)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no northpole!"
    end if
    if(abs(xprobmax2-dpi)<smalldouble) then
      if(mype==0) write(unitterm,*) &
        "Will apply pi-periodic conditions at southpole!"
      poleB(2,2)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no southpole!"
    end if
  end if}
case ("cylindrical")
  {if (^D==^PHI.and.periodB(^D)) then
    if(mod(ng^D(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin1)<smalldouble) then
      if(mype==0) write(unitterm,*) "Will apply pi-periodic conditions at r=0"
      poleB(1,1)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no cylindrical axis!"
    end if
  end if\}
end select

end subroutine set_pole
!=============================================================================
subroutine getgridgeo(igrid)

use mod_metric
include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: ix^L, ixCoG^L, ixCoM^L, ixCo^L, ixCoCoG^L, ixGext^L
double precision :: xmin^D, dx^D
!-----------------------------------------------------------------------------
!!!!!!!!!!!!!!!! For staggered constrained transport, surfaces must be defined
!!!!!!!!!!!!!!!! also on the faces with index 0. This at the moment is done only
!!!!!!!!!!!!!!!! in fillgeo_covariant
!ix^L=ixM^LL^LADD1;

ix^L=ixG^LL^LSUB1;
if (2*int(dixB/2)==dixB) then
   ixGext^L=ixG^LL;
else
   ixGext^L=ixG^LL^LADD1;
end if
allocate(pgeo(igrid)%surfaceC^D(ixGlo^D-1:ixGhi^D^D%ixG^T), &
         pgeo(igrid)%surface^D(ixmin^D-1:ixmax^D^D%ix^S), &
         pgeo(igrid)%dvolume(ixGext^S), &
         pgeo(igrid)%dx(ixGext^S,1:ndim), &
         pgeo(igrid)%xbar(ixGext^S,1:ndim))
! Covariant coordinates rely on the metric datastructure
if (covariant) then
     call allocate_metric(pgeo(igrid)%m,ixGext^L)
     {^D& call allocate_metric(pgeo(igrid)%mSurface^D,{ixGextmin^DD-1},ixGextmax^DD,need_derivs=.false.)\}
end if

dx^D=rnode(rpdx^D_,igrid);
xmin^D=rnode(rpxmin^D_,igrid);
call fillgeo(pgeo(igrid),pw(igrid),ixG^LL,ixGext^L,xmin^D,dx^D,.false.)

{#IFDEF DY_SP
  call deallocate_metric(pgeo(igrid)%m)
  {^D& call deallocate_metric(pgeo(igrid)%mSurface^D) \}
}

if (errorestimate==1) then
call mpistop('stop using errorestimate==1')
   ixCoGmin^D=1; ixCoGmax^D=ixGhi^D/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGext^L=ixCoG^L;
   else
      ixGext^L=ixCoG^L^LADD1;
   end if
   ixCoM^L=ixCoG^L^LSUBdixB;
   ixCo^L=ixCoM^L^LADD1;
   ixCo^L=ixCoG^L^LSUB1;

   allocate(pgeoCoarse(igrid)%surfaceC^D(ixComin^D-1:ixComax^D^D%ixCo^S), &
        pgeoCoarse(igrid)%surface^D(ixComin^D-1:ixComax^D^D%ixCo^S), &
        pgeoCoarse(igrid)%dvolume(ixGext^S), &
        pgeoCoarse(igrid)%dx(ixGext^S,1:ndim), &
        pgeoCoarse(igrid)%xbar(ixGext^S,1:ndim))

   if (covariant) then
      call allocate_metric(pgeoCoarse(igrid)%m,ixGext^L)
      {^D& call allocate_metric(pgeoCoarse(igrid)%mSurface^D,{ixGextmin^DD-1},ixGextmax^DD,need_derivs=.false.)\}
   end if

   dx^D=two*rnode(rpdx^D_,igrid);
   ! Warning: if you do not use pw(igrid) instead of pwCoarse, bc_gc_r will copy a wrong pole
   call fillgeo(pgeoCoarse(igrid),pw(igrid),ixCoG^L,ixGext^L,xmin^D,dx^D,.false.)  

   ixCoCoGmin^D=1; ixCoCoGmax^D=ixCoGmax^D/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGext^L=ixCoCoG^L;
   else
      ixGext^L=ixCoCoG^L^LADD1;
   end if


   allocate(pgeoCoCo(igrid)%dvolume(ixGext^S), &
        pgeoCoCo(igrid)%xbar(ixGext^S,1:ndim))
   if (covariant) &
        call allocate_metric(pgeoCoCo(igrid)%m,ixGext^L)

   dx^D=4.0d0*rnode(rpdx^D_,igrid);
   call fillgeo(pgeoCoCo(igrid),pwCoCo(igrid),ixCoCoG^L,ixGext^L,xmin^D,dx^D,.true.)
else
   ixCoGmin^D=1; ixCoGmax^D=ixGhi^D/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGext^L=ixCoG^L;
   else
      ixGext^L=ixCoG^L^LADD1;
   end if
{#IFNDEF STAGGERED
   allocate(pgeoCoarse(igrid)%dvolume(ixGext^S), &
        pgeoCoarse(igrid)%xbar(ixGext^S,1:ndim))
   if (covariant) &
        call allocate_metric(pgeoCoarse(igrid)%m,ixGext^L)

   dx^D=two*rnode(rpdx^D_,igrid);
   ! Warning: if you do not use pw(igrid) instead of pwCoarse, bc_gc_r will copy a wrong pole
   call fillgeo(pgeoCoarse(igrid),pw(igrid),ixCoG^L,ixGext^L,xmin^D,dx^D,.true.)
}
{#IFDEF STAGGERED
   ixCoM^L=ixCoG^L^LSUBdixB;
   ixCo^L=ixCoM^L^LADD1;
   ixCo^L=ixCoG^L^LSUB1;

   allocate(pgeoCoarse(igrid)%surfaceC^D(ixCoGmin^D-1:ixCoGmax^D^D%ixCoG^S),&
        pgeoCoarse(igrid)%surface^D(ixComin^D-1:ixComax^D^D%ixCo^S),&
        pgeoCoarse(igrid)%dvolume(ixGext^S),&
        pgeoCoarse(igrid)%dx(ixGext^S,1:ndim),&
        pgeoCoarse(igrid)%xbar(ixGext^S,1:ndim))

   if (covariant) then
      call allocate_metric(pgeoCoarse(igrid)%m,ixGext^L)
      {^D& call allocate_metric(pgeoCoarse(igrid)%mSurface^D,{ixGextmin^DD-1},ixGextmax^DD,need_derivs=.false.)\}
   end if

   dx^D=two*rnode(rpdx^D_,igrid);
   ! Warning: if you do not use pw(igrid) instead of pwCoarse, bc_gc_r will copy a wrong pole
   call fillgeo(pgeoCoarse(igrid),pw(igrid),ixCoG^L,ixGext^L,xmin^D,dx^D,.false.)
}
{#IFDEF DY_SP
  call deallocate_metric(pgeoCoarse(igrid)%m)
  {#IFDEF STAGGERED
  {^D& call deallocate_metric(pgeoCoarse(igrid)%mSurface^D) \}
  }
}

end if

end subroutine getgridgeo
!=============================================================================
subroutine dealloc_myM(igrid)
  use mod_metric
  include 'amrvacdef.f'
integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (covariant) then
{#IFDEF DY_SP
   call deallocate_metric(pgeo(igrid)%m)
   call deallocate_metric(pgeoCoarse(igrid)%m)
   {^D& call deallocate_metric(pgeo(igrid)%mSurface^D)\}
   
   {#IFDEF STAGGERED
   {^D& call deallocate_metric(pgeoCoarse(igrid)%mSurface^D)\}
   }
}
{#IFNDEF DY_SP
  call mpistop('not define DY_SP, should not dealloc_myM')
}
endif 
end subroutine dealloc_myM

!=============================================================================
subroutine putgridgeo(igrid)

  use mod_metric
  include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
! Warning: AMR refine or coarse will call here. DY_SP should not deallocate metric again
{#IFNDEF DY_SP
if (covariant) then
   call deallocate_metric(pgeo(igrid)%m)
   call deallocate_metric(pgeoCoarse(igrid)%m)
   {^D& call deallocate_metric(pgeo(igrid)%mSurface^D)\}
end if
}

if (errorestimate==1) then
   deallocate(pgeo(igrid)%surfaceC^D,pgeo(igrid)%surface^D,&
        pgeo(igrid)%dvolume,pgeo(igrid)%dx, &
        pgeo(igrid)%xbar, &
        pgeoCoarse(igrid)%surfaceC^D,pgeoCoarse(igrid)%surface^D,&
        pgeoCoarse(igrid)%dvolume,pgeoCoarse(igrid)%dx,&
        pgeoCoarse(igrid)%xbar,&
        pgeoCoCo(igrid)%dvolume,&
        pgeoCoCo(igrid)%xbar)
   
   if (covariant) then
      call deallocate_metric(pgeoCoCo(igrid)%m)
      {^D& call deallocate_metric(pgeoCoarse(igrid)%mSurface^D)\}
   end if
   
else
   deallocate(pgeo(igrid)%surfaceC^D,pgeo(igrid)%surface^D,&
        pgeo(igrid)%dvolume,pgeo(igrid)%dx,pgeo(igrid)%xbar,&
        pgeoCoarse(igrid)%dvolume,pgeoCoarse(igrid)%xbar)
{#IFDEF STAGGERED
   deallocate(pgeoCoarse(igrid)%surfaceC^D,pgeoCoarse(igrid)%surface^D,&
        pgeoCoarse(igrid)%dx)
   {#IFNDEF DY_SP
   if (covariant) then
      {^D& call deallocate_metric(pgeoCoarse(igrid)%mSurface^D)\}
   end if
   }
}
end if

end subroutine putgridgeo
!=============================================================================
subroutine fillgeo(pgeogrid,pwgrid,ixG^L,ixGext^L,xmin^D,dx^D,need_only_volume)
!subroutine fillgeo(pgeogrid,psgrid,ixG^L,ixGext^L,xmin^D,dx^D,need_only_volume)

use mod_metric
include 'amrvacdef.f'

type(geoalloc) :: pgeogrid
type(walloc)    :: pwgrid
!type(state)    :: psgrid
integer, intent(in) :: ixG^L, ixGext^L
double precision, intent(in) :: xmin^D, dx^D
logical, intent(in) :: need_only_volume

integer :: idims, ix, ixM^L, ix^L, ixC^L
double precision :: x(ixGext^S,ndim)
!-----------------------------------------------------------------------------
!associate(x=>psgrid%x%x,w=>psgrid%w%w{#IFDEF STAGGERED ,ws=>psgrid%ws%w})

ixM^L=ixG^L^LSUBdixB;
!ix^L=ixM^L^LADD1;
ix^L=ixG^L^LSUB1;
! write(*,*) 'inside fillgeo' 
if (covariant) then
   ! call for the general geometry
   call fillgeo_covariant(pgeogrid,pwgrid,ixG^L,ixGext^L,xmin^D,dx^D,need_only_volume)
   !call fillgeo_covariant(pgeogrid,psgrid,ixG^L,ixGext^L,xmin^D,dx^D,need_only_volume)
   return
end if


call mpistop('Dont pass the below')
select case (typecoord)
case ("slabtest")

   pgeogrid%dvolume(ixGext^S) = ^D&dx^D*

   if (need_only_volume) return

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   pgeogrid%surfaceC1(ixC^S)={^IFONED one}{^NOONED dx2}{^IFTHREED*dx3}
   pgeogrid%surface1(ixC^S) ={^IFONED one}{^NOONED dx2}{^IFTHREED*dx3}
   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   pgeogrid%surfaceC2(ixC^S)=dx1}{^IFTHREED*dx3}
   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   pgeogrid%surface2(ixC^S)=dx1}{^IFTHREED*dx3}
   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   pgeogrid%surfaceC3(ixC^S)=dx1*dx2
   pgeogrid%surface3(ixC^S)=dx1*dx2}

   ^D&pgeogrid%dx(ixGext^S,^D)=dx^D;
   ! Not yet correct:
   ^D&pgeogrid%xbar(ixGext^S,^D)=0.0d0;


case ("spherical")

   do idims=1,min(ndim,2)
      select case(idims)
      {case(^D)
         do ix = ixGext^LIM^D
            x(ix^D%ixGext^S,^D)=xmin^D+(dble(ix-dixB)-half)*dx^D
         end do\}
      end select
   end do

   if(typespherical==0) then
     pgeogrid%dvolume(ixGext^S)=(x(ixGext^S,1)**2+dx1**2/12.0d0)*dx1 {^NOONED &
              *two*dabs(dsin(x(ixGext^S,2)))*dsin(half*dx2)}{^IFTHREED*dx3}
   else
     pgeogrid%dvolume(ixGext^S)=(x(ixGext^S,1)**2)*dx1 {^NOONED &
              *dabs(dsin(x(ixGext^S,2)))*dx2}{^IFTHREED*dx3}
   endif

   if (need_only_volume) return

   ! Not yet correct:
   ^D&pgeogrid%xbar(ixGext^S,^D)=0.0d0;

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   if(typespherical==0) then
       pgeogrid%surfaceC1(ixC^S)=(x(ixC^S,1)+half*dx1)**2 {^NOONED &
              *two*dsin(x(ixC^S,2))*dsin(half*dx2)}{^IFTHREED*dx3}
   else
       pgeogrid%surfaceC1(ixC^S)=(x(ixC^S,1)+half*dx1)**2 {^NOONED &
              *dsin(x(ixC^S,2))*dx2}{^IFTHREED*dx3}
   endif

   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   pgeogrid%surfaceC2(ixC^S)=x(ixC^S,1)*dx1 &
              *dsin(x(ixC^S,2)+half*dx2)}{^IFTHREED*dx3}

   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   pgeogrid%surfaceC3(ixC^S)=x(ixC^S,1)*dx1*dx2}

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   if(typespherical==0) then
       pgeogrid%surface1(ixC^S)=x(ixC^S,1)**2 {^NOONED &
              *two*dsin(x(ixC^S,2))*dsin(half*dx2)}{^IFTHREED*dx3}
   else
      pgeogrid%surface1(ixC^S)=x(ixC^S,1)**2 {^NOONED &
              *dsin(x(ixC^S,2))*dx2}{^IFTHREED*dx3}
   endif

   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   pgeogrid%surface2(ixC^S)=x(ixC^S,1)*dx1 &
              *dsin(x(ixC^S,2))}{^IFTHREED*dx3}

   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   pgeogrid%surface3(ixC^S)=x(ixC^S,1)*dx1*dx2}

   pgeogrid%dx(ixGext^S,1)=dx1
   {^NOONED pgeogrid%dx(ixGext^S,2)=x(ixGext^S,1)*dx2}
   {^IFTHREED pgeogrid%dx(ixGext^S,3)=x(ixGext^S,1)*dsin(x(ixGext^S,2))*dx3}

case ("cylindrical")

   do ix = ixGext^LIM1
      x(ix,ixGext^SE,1)=xmin1+(dble(ix-dixB)-half)*dx1
   end do

   pgeogrid%dvolume(ixGext^S)=dabs(x(ixGext^S,1))*{^D&dx^D*}
   pgeogrid%dvolume(ixGext^S)=dabs(half*((x(ixGext^S,1)+half*dx1)**2-(x(ixGext^S,1)-half*dx1)**2)){^DE&*dx^DE }

   if (need_only_volume) return

   ! Not yet correct:
   ^D&pgeogrid%xbar(ixGext^S,^D)=0.0d0;

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   !!pgeogrid%surfaceC1(ixC^S)=(x(ixC^S,1)+half*dx1){^DE&*dx^DE }
   pgeogrid%surfaceC1(ixC^S)=dabs((x(ixC^S,1)+half*dx1)){^DE&*dx^DE }
   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   if (^Z==2) pgeogrid%surfaceC2(ixC^S)=x(ixC^S,1)*dx1{^IFTHREED*dx3}
   if (^PHI==2) pgeogrid%surfaceC2(ixC^S)=dx1{^IFTHREED*dx3}}
   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   if (^Z==3) pgeogrid%surfaceC3(ixC^S)=x(ixC^S,1)*dx1*dx2
   if (^PHI==3) pgeogrid%surfaceC3(ixC^S)=dx1*dx2}

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   !!pgeogrid%surface1(ixC^S)=x(ixC^S,1){^DE&*dx^DE }
   pgeogrid%surface1(ixC^S)=dabs(x(ixC^S,1)){^DE&*dx^DE }
   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   if (^Z==2) pgeogrid%surface2(ixC^S)=x(ixC^S,1)*dx1{^IFTHREED*dx3}
   if (^PHI==2) pgeogrid%surface2(ixC^S)=dx1{^IFTHREED*dx3}}
   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   if (^Z==3) pgeogrid%surface3(ixC^S)=x(ixC^S,1)*dx1*dx2
   if (^PHI==3) pgeogrid%surface3(ixC^S)=dx1*dx2}


   pgeogrid%dx(ixGext^S,1)=dx1
   {^IFZ {^DE&if (^DE==^Z) pgeogrid%dx(ixGext^S,^DE)=dx^DE\}}
   {^IFPHI {if (^DE==^PHI) pgeogrid%dx(ixGext^S,^DE)=x(ixGext^S,1)*dx^DE\}}

case default

   call mpistop("Sorry, typecoord unknown")
   
end select
!end associate
end subroutine fillgeo
!=============================================================================
subroutine gradient(q,ix^L,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir

include 'amrvacdef.f'

integer :: ix^L, idir
double precision :: q(ixG^T), gradq(ixG^T)

double precision :: qC(ixG^T),invdx
integer :: jx^L, hx^L, ixC^L, jxC^L {#IFDEF FOURTHORDER , lx^L, kx^L}

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
if (slab) then
{#IFNDEF FOURTHORDER
   jx^L=ix^L+kr(idir,^D);
   hx^L=ix^L-kr(idir,^D);
   gradq(ix^S) = half*(q(jx^S)-q(hx^S))*invdx
}{#IFDEF FOURTHORDER
   lx^L=ix^L+2*kr(idir,^D);
   jx^L=ix^L+kr(idir,^D);
   hx^L=ix^L-kr(idir,^D);
   kx^L=ix^L-2*kr(idir,^D);
   gradq(ix^S) = (-q(lx^S) + 8.0d0 * q(jx^S) - 8.0d0 * q(hx^S) + q(kx^S)) &
        /(12.0d0 * dxlevel(idir))
}
else
   hx^L=ix^L-kr(idir,^D);
   ixCmin^D=hxmin^D;ixCmax^D=ixmax^D;
   jxC^L=ixC^L+kr(idir,^D);
   select case(idir)
   {case(^D)
      qC(ixC^S)=mygeo%surfaceC^D(ixC^S)*half*(q(ixC^S)+q(jxC^S))
      gradq(ix^S)=(qC(ix^S)-qC(hx^S))/mygeo%dvolume(ix^S)
      ! Substract difference divergence and gradient
      gradq(ix^S)=gradq(ix^S)-q(ix^S) &
                     *(mygeo%surfaceC^D(ix^S)-mygeo%surfaceC^D(hx^S)) &
                    /mygeo%dvolume(ix^S) \}
   end select
end if

end subroutine gradient
!=============================================================================
subroutine upwindGradientS(ixI^L,ixO^L,idir,w,x,var,gradient)
!dont know why you need this, fixme

  use mod_limiter
  include 'amrvacdef.f'
  
  integer, intent(in)             :: ixI^L, ixO^L, idir
  double precision, intent(in)    :: w(ixI^S,1:nw)
  double precision, intent(in)    :: x(ixI^S,1:ndim)
  double precision, intent(in)    :: var(ixI^S)
  double precision, intent(out)   :: gradient(ixI^S)
  ! .. local ..
  integer                                 :: hxO^L, ixC^L, jxC^L, gxC^L, hxC^L
  double precision, dimension(ixG^T,1:ndim)             ::  xi
  double precision, dimension(ixG^T,1:nw) :: wLC, wRC, wLp, wRp, wprim
  double precision, dimension(ixG^T)      :: cmaxC, cmaxRC, cmaxLC
  double precision, dimension(ixG^T)      :: cminC, cminRC, cminLC
  integer, dimension(ixG^T)               :: patchf
  double precision,dimension(ixI^S)       :: qL, qR, q, ldq, rdq, dqC
  logical          :: fattening = .False.

  !-----------------------------------------------------------------------------
  ! Always primitive reconstruction:
  wprim(ixI^S,1:nw) = w(ixI^S,1:nw)
call mpistop('not changed yet for upwindGradientS')
  call primitive(ixI^L,ixI^L,wprim,x) ! also applies floors to wprim


  ! Index ranges:
  hxO^L=ixO^L-kr(idir,^D);
  ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
  jxC^L=ixC^L+kr(idir,^D);
  gxCmin^D=ixCmin^D-kr(idir,^D);gxCmax^D=jxCmax^D;
  hxC^L=gxC^L+kr(idir,^D);
  
  !==================================================
  ! left and right states:
  !==================================================
  wRp(ixC^S,1:nwflux)=wprim(jxC^S,1:nwflux)
  wLp(ixC^S,1:nwflux)=wprim(ixC^S,1:nwflux)
  xi(ixC^S,idir) = half * ( x(ixC^S,idir)+x(jxC^S,idir) )

  call upwindLR(ixI^L,ixC^L,ixC^L,idir,wprim,wprim,wLC,wRC,wLp,wRp,x,dxlevel(idir))
  ! get auxilaries for L and R states
  !if (nwaux>0) then
  !   call getaux(.true.,wLC,xi,ixG^LL,ixC^L,'upwindGradientS-L')
  !   call getaux(.true.,wRC,xi,ixG^LL,ixC^L,'upwindGradientS-R')
  !end if
  call getcmax(wLC,xi,ixG^LL,ixC^L,idir,cmaxLC,cminLC,.true.)
  call getcmax(wRC,xi,ixG^LL,ixC^L,idir,cmaxRC,cminRC,.true.)
  ! now take the maximum of left and right states
  cmaxC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))
  cminC(ixC^S)=min(cminRC(ixC^S),cminLC(ixC^S))

  ! select cases based on fasted waves:
  patchf(ixC^S) =  1
  where(cminC(ixC^S) >= zero)
     patchf(ixC^S) = -2
  elsewhere(cmaxC(ixC^S) <= zero)
     patchf(ixC^S) =  2
  endwhere
  !==================================================
  !==================================================
  
  !==================================================
  ! reconstruct the variable itself, from left and right
  !==================================================
  qR(gxC^S) = var(hxC^S)
  qL(gxC^S) = var(gxC^S)

  select case (typegradlimiter)
     
  case (limiter_ppm)
     fattening = .True. 
     call PPMlimiter(ixI^L,ixM^LL,idir,var,qL,qR,PPM_extrema)
     !call PPMlimiter(ixI^L,ixM^LL,idir,var,var,qL,qR,fattening)
    
     !call PPMlimitervar(ixI^L,ixM^LL,idir,var,var,qL,qR)

  case default

     dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)

     call dwlimiter2(dqC,ixI^L,gxC^L,idir,typegradlimiter,ldq,rdq)

     qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
     qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)

  end select
  !==================================================


  ! Sl>=0, use left interface value:
  where (patchf(ixC^S) .eq. -2)
     q(ixC^S) = qL(ixC^S)     
  ! Sr<=0, use right interface value:
  elsewhere  (patchf(ixC^S) .eq. 2)
     q(ixC^S) = qR(ixC^S)
  ! Outgoing in both directions, use average value:
  elsewhere
     q(ixC^S) = half * (qL(ixC^S) + qR(ixC^S))
  end where

  if (slab) then
     gradient(ixO^S)=half*(q(ixO^S)-q(hxO^S))/dxlevel(idir)
  else
     select case(idir)
        {case(^D)
        gradient(ixO^S)=(q(ixO^S)-q(hxO^S))/mygeo%dx(ixO^S,idir) \}
     end select
  end if

end subroutine upwindGradientS
!=============================================================================
subroutine gradientS(q,ix^L,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir
! first use limiter to go from cell center to edge

use mod_limiter
include 'amrvacdef.f'

integer :: ix^L, idir
double precision :: q(ixG^T), gradq(ixG^T)
double precision :: dxdim

double precision :: qC(ixG^T)
double precision,dimension(ixG^T):: qL,qR,dqC,ldq,rdq
integer                          :: hx^L,ixC^L,jxC^L,gxC^L,hxC^L
logical          :: fattening = .False.

!-----------------------------------------------------------------------------


hx^L=ix^L-kr(idir,^D);
ixCmin^D=hxmin^D;ixCmax^D=ixmax^D;
jxC^L=ixC^L+kr(idir,^D);
gxCmin^D=ixCmin^D-kr(idir,^D);gxCmax^D=jxCmax^D;
hxC^L=gxC^L+kr(idir,^D);


!==================================================
! reconstruct the variable itself, from left and right
!==================================================
qR(gxC^S) = q(hxC^S)
qL(gxC^S) = q(gxC^S)

select case (typegradlimiter)

case (limiter_ppm)
   fattening = .True. 
   call PPMlimiter(ixG^LL,ixM^LL,idir,q,qL,qR,PPM_extrema)
   !call PPMlimiter(ixG^LL,ixM^LL,idir,q,q,qL,qR,fattening)
   !call PPMlimitervar(ixG^LL,ixM^LL,idir,q,q,qL,qR)

case default

   dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)

   call dwlimiter2(dqC,ixG^LL,gxC^L,idir,typegradlimiter,ldq,rdq)

   qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
   qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)

end select
!==================================================


if (slab) then
   gradq(ix^S)=half*(qR(ix^S)-qL(hx^S))/dxlevel(idir)
else
   select case(idir)
   {case(^D)
    gradq(ix^S)=(qR(ix^S)-qL(hx^S))/mygeo%dx(ix^S,idir) \}
   end select
end if

end subroutine gradientS
!=============================================================================
subroutine divvector(qvec,ixI^L,ixO^L,divq)

! Calculate divergence of a vector qvec within ixL

include 'amrvacdef.f'

integer :: ixI^L,ixO^L
double precision :: qvec(ixG^T,1:ndir), divq(ixG^T)

double precision :: qC(ixG^T), invdx(1:ndim)
integer :: jxO^L, hxO^L, ixC^L, jxC^L, idims, ix^L {#IFDEF FOURTHORDER , gxO^L, kxO^L}
!-----------------------------------------------------------------------------
{#IFNDEF FOURTHORDER
ix^L=ixO^L^LADD1;
}{#IFDEF FOURTHORDER
ix^L=ixO^L^LADD2;
}
if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in divvector: Non-conforming input limits")
invdx=1.d0/dxlevel
divq(ixO^S)=zero
do idims=1,ndim
   if (slab) then
{#IFNDEF FOURTHORDER
     jxO^L=ixO^L+kr(idims,^D);
     hxO^L=ixO^L-kr(idims,^D);
     divq(ixO^S)=divq(ixO^S)+half*(qvec(jxO^S,idims)-qvec(hxO^S,idims))*invdx(idims)
}{#IFDEF FOURTHORDER
      kxO^L=ixO^L+2*kr(idims,^D);
      jxO^L=ixO^L+kr(idims,^D);
      hxO^L=ixO^L-kr(idims,^D);
      gxO^L=ixO^L-2*kr(idims,^D);
      divq(ixO^S)=divq(ixO^S)+&
           (-qvec(kxO^S,idims) + 8.0d0 * qvec(jxO^S,idims) - 8.0d0 * qvec(hxO^S,idims) + qvec(gxO^S,idims))/(12.0d0 * dxlevel(idims))
}
   else
     hxO^L=ixO^L-kr(idims,^D);
     ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
     jxC^L=ixC^L+kr(idims,^D);
     select case(idims)
     {case(^D)
        qC(ixC^S)=mygeo%surfaceC^D(ixC^S)*half*(qvec(ixC^S,idims)+qvec(jxC^S,idims))
        divq(ixO^S)=divq(ixO^S)+(qC(ixO^S)-qC(hxO^S))/mygeo%dvolume(ixO^S) \}
      end select
   end if
end do


end subroutine divvector 
!=============================================================================
subroutine curl3(qvec,ixI^L,ixO^L,curlvec)
! To calculate the curl 'curlvec' of a covariant 3-vector 'qvec'
! in the range ixO^L. The total extent ixI^L must be at least 1
! cell bigger than ixO^L in all directions.
! The curl is an array of three components.
! A pointer myM must be associated to the metric.
! Derivatives are calculated with a 2nd order approximation and
! without staggering.
! This routine works only for vectors in coordinate basis.

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: qvec(ixI^S,1:3)
double precision, intent(out) :: curlvec(ixI^S,1:3)

! ... local ...

integer :: idir1,idir2,idir3
integer :: jxO^L,hxO^L,ixOp^L,ixOm^L
double precision :: invdx(1:^ND)
!-----------------------------------------------------------------------------
{#IFDEF DY_SP
call mpistop('DY_SP enters curl3')
}

invdx=1.d0/dxlevel

curlvec=zero

do idir1=1,3
  ! Directions
  idir2=modulo(idir1,3)+1
  idir3=modulo(idir2,3)+1

  ! Indices
  hxO^L=ixO^L-kr(idir2,^D);
  jxO^L=ixO^L+kr(idir2,^D);
  ixOm^L=ixO^L-kr(idir3,^D);
  ixOp^L=ixO^L+kr(idir3,^D);

  if (idir2.le.^ND) &
    curlvec(ixO^S,idir1) = curlvec(ixO^S,idir1) &
     + invdx(idir2)*half*(qvec(jxO^S,idir3)-qvec(hxO^S,idir3))

  if (idir3.le.^ND) &
    curlvec(ixO^S,idir1) = curlvec(ixO^S,idir1) &
     - invdx(idir3)*half*(qvec(ixOp^S,idir2)-qvec(ixOm^S,idir2))

  !! When covariant, divide between metric element
  if (covariant) curlvec(ixO^S,idir1)=&
    curlvec(ixO^S,idir1)/myM%sqrtgamma(ixO^S)

end do

end subroutine curl3
!=============================================================================
subroutine curlvector(qvec,ixI^L,ixO^L,curlvec,idirmin,idirmin0,ndir0)

! Calculate curl of a vector qvec within ixL

include 'amrvacdef.f'

integer :: ixI^L,ixO^L,idirmin,ix^L,idir,jdir,kdir,hxO^L,jxO^L,ndir0,idirmin0{#IFDEF FOURTHORDER , gxO^L, kxO^L}
double precision :: qvec(ixG^T,1:ndir0),curlvec(ixG^T,idirmin0:3), invdx(1:ndim)
double precision :: tmp(ixG^T),tmp2(ixG^T),surface(ixG^T),mydx(ixG^T)
!-----------------------------------------------------------------------------
{#IFNDEF FOURTHORDER
ix^L=ixO^L^LADD1;
}{#IFDEF FOURTHORDER
ix^L=ixO^L^LADD2;
}
if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in curl: Non-conforming input limits")

! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
! Curl can have components (idirmin0:3)
! Determine exact value of idirmin while doing the loop.

invdx=1.d0/dxlevel
idirmin=4
curlvec(ixO^S,idirmin0:3)=zero

do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ix^S)=qvec(ix^S,kdir)
      hxO^L=ixO^L-kr(jdir,^D);
      jxO^L=ixO^L+kr(jdir,^D);
      if(slab)then
{#IFDEF FOURTHORDER
      kxO^L=ixO^L+2*kr(jdir,^D);
      gxO^L=ixO^L-2*kr(jdir,^D);
      tmp2(ixO^S)=(-tmp(kxO^S) + 8.0d0 * tmp(jxO^S) - 8.0d0 * tmp(hxO^S) + tmp(gxO^S)) &
           / (12.0d0 * dxlevel(jdir))
}
{#IFNDEF FOURTHORDER
         tmp2(ixO^S)=half*(tmp(jxO^S)-tmp(hxO^S))*invdx(jdir)
}
      else
         ! approximate formula, reduces to slab case
         ! and avoids staggering

         if (kdir .le. ndim) then 
            mydx(ix^S)=mygeo%dx(ix^S,kdir)
         else 
            mydx(ix^S)=one
         end if

         select case(idir)
           {case(^D)
             surface(ixO^S)=mygeo%surface^D(ixO^S)
             tmp2(ixO^S)=half*(mydx(jxO^S)*tmp(jxO^S) &
                              -mydx(hxO^S)*tmp(hxO^S)) &
                     /surface(ixO^S) \}
          end select
      endif
      if(lvc(idir,jdir,kdir)==1)then
         curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
      else
         curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
      endif
      if(idir<idirmin)idirmin=idir
   endif
enddo; enddo; enddo;

end subroutine curlvector 
!=============================================================================
subroutine divvectorS(qvec,ixI^L,ixO^L,divq)

! Calculate divergence of a vector qvec within ixL
! using limited extrapolation to cell edges

use mod_limiter
include 'amrvacdef.f'

integer :: ixI^L,ixO^L
double precision :: qvec(ixI^S,1:ndir), divq(ixI^S)

double precision,dimension(ixI^S):: qL,qR,dqC,ldq,rdq
double precision :: dxdim, invdx(1:ndim)

integer :: hxO^L,ixC^L,jxC^L,idims,ix^L,gxC^L,hxC^L,idummy
character*79, save :: savetypelimiter,savetypegradlimiter,save2typelimiter
logical          :: fattening = .False.
!-----------------------------------------------------------------------------
ix^L=ixO^L^LADD2;

if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in divvectorS: Non-conforming input limits")

idummy=0
invdx=1.d0/dxlevel
divq(ixO^S)=zero
do idims=1,ndim
   hxO^L=ixO^L-kr(idims,^D);
   ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
   jxC^L=ixC^L+kr(idims,^D);
   gxCmin^D=ixCmin^D-kr(idims,^D);gxCmax^D=jxCmax^D;
   hxC^L=gxC^L+kr(idims,^D);

   qR(gxC^S) = qvec(hxC^S,idims)
   qL(gxC^S) = qvec(gxC^S,idims)

  !==================================================
  ! reconstruct the variable itself, from left and right
  !==================================================
  select case (typegradlimiter)
     
  case (limiter_ppm)
     
     dqC(ixI^S)=qvec(ixI^S,idims)
     fattening = .True. 
     call PPMlimiter(ixI^L,ixO^L,idims,dqC,qL,qR,PPM_extrema)
     !call PPMlimiter(ixI^L,ixO^L,idims,dqC,dqC,qL,qR,fattening)
     !call PPMlimitervar(ixI^L,ixO^L,idims,dqC,dqC,qL,qR)

  case default

     dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)

     call dwlimiter2(dqC,ixI^L,gxC^L,idims,typegradlimiter,ldq,rdq)

     qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
     qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)

  end select
  !==================================================

   if (slab) then
     divq(ixO^S)=divq(ixO^S)+half*(qR(ixO^S)-qL(hxO^S))*invdx(idims)
   else
     select case(idims)
     {case(^D)
        qR(ixC^S)=mygeo%surfaceC^D(ixC^S)*qR(ixC^S)
        qL(ixC^S)=mygeo%surfaceC^D(ixC^S)*qL(ixC^S)
        divq(ixO^S)=divq(ixO^S)+(qR(ixO^S)-qL(hxO^S))/mygeo%dvolume(ixO^S) \}
      end select
   end if
end do

end subroutine divvectorS
!=============================================================================
subroutine rec(ixI^L,ixC^L,idir,q,qL,qR)

  ! Reconstruct scalar q within ixO^L to 1/2 dx in direction idir
  ! Return both left and right reconstructed values 

use mod_limiter
include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixC^L, idir
double precision, intent(in)       :: q(ixI^S)
double precision, intent(out)      :: qL(ixI^S), qR(ixI^S)

double precision                   :: qC(ixI^S)
double precision,dimension(ixI^S)  :: dqC,ldq,rdq
integer                            :: ixO^L,jxC^L,gxC^L,hxC^L
logical          :: fattening = .False.
!-----------------------------------------------------------------------------

jxC^L=ixC^L+kr(idir,^D);
gxCmin^D=ixCmin^D-kr(idir,^D);gxCmax^D=jxCmax^D;
hxC^L=gxC^L+kr(idir,^D);

qR(gxC^S) = q(hxC^S)
qL(gxC^S) = q(gxC^S)

select case (typelimiter)
   
case (limiter_ppm)
   ! the ordinary grid-index:
   ixOmin^D=ixCmin^D+kr(idir,^D);
   ixOmax^D=ixCmax^D;
   !fattening = .True. 
   call PPMlimiter(ixI^L,ixC^L,idir,q,qL,qR,PPM_extrema)
   !call PPMlimiter(ixI^L,ixC^L,idir,q,q,qL,qR,fattening)
   !call PPMlimitervar(ixI^L,ixO^L,idir,q,q,qL,qR)
   
case (limiter_mp5)
   call MP5limiter(ixI^L,ixC^L,idir,q,qL,qR)
   !call MP5limitervar(ixI^L,ixC^L,idir,q,qL,qR)
case (limiter_weno5)
   call WENO5limitervar(ixI^L,ixC^L,idir,q,qL,qR)
case (limiter_wenoZ5)
   call WENOZ5limitervar(ixI^L,ixC^L,idir,dxlevel(idir),q,qL,qR)
case (limiter_wenoZP)
   call WENOZPlimitervar(ixI^L,ixC^L,idir,dxlevel(idir),q,qL,qR)
case default

   dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
   
   
   call dwlimiter2(dqC,ixI^L,gxC^L,idir,typelimiter,ldq,rdq)
   
   qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
   qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)
   
end select


end subroutine rec

  !> get the 3-metric in flat space gamma_ij_hat
  subroutine get_gamma_ij_hat(x, ixI^L, ixO^L, gamma)
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(out)   :: gamma(ixI^S, 1:3, 1:3)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    gamma = 0.0d0

    select case (coordinate)
    case (cartesian)
      gamma(ixO^S,1,1) = 1.0d0
      gamma(ixO^S,2,2) = 1.0d0
      gamma(ixO^S,3,3) = 1.0d0
    case (cylindrical)
      ! note: r_ = 1,z_ = 2, phi_ = 3.
      gamma(ixO^S,1,1) = 1.0d0
      gamma(ixO^S,2,2) = 1.0d0
      gamma(ixO^S,3,3) = x(ixO^S, 1)**2
    case (spherical)
      ! note: r_ = 1,theta_=2, phi_ = 3.
      gamma(ixO^S,1,1) = 1.0d0
      gamma(ixO^S,2,2) = x(ixO^S, 1)**2
      gamma(ixO^S,3,3) = gamma(ixO^S,2,2) {^NOONED * dsin(x(ixO^S, 2))**2 }
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_gamma_ij_hat
  !> get the 3-metric in flat space gammainvij_hat
  subroutine get_gammainvij_hat(x, ixI^L, ixO^L, gammainv)
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(out)   :: gammainv(ixI^S, 1:3, 1:3)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    gammainv = 0.0d0

    select case (coordinate)
    case (cartesian)
      gammainv(ixO^S,1,1) = 1.0d0
      gammainv(ixO^S,2,2) = 1.0d0
      gammainv(ixO^S,3,3) = 1.0d0
    case (cylindrical)
      ! note: r_ = 1,z_ = 2, phi_ = 3.
      gammainv(ixO^S,1,1) = 1.0d0
      gammainv(ixO^S,2,2) = 1.0d0
      gammainv(ixO^S,3,3) = 1.0d0 / x(ixO^S, 1)**2
    case (spherical)
      ! note: r_ = 1,theta_=2, phi_ = 3.
      gammainv(ixO^S,1,1) = 1.0d0
      gammainv(ixO^S,2,2) = 1.0d0 / x(ixO^S, 1)**2
      gammainv(ixO^S,3,3) = 1.0d0 / (x(ixO^S, 1)**2  {^NOONED * dsin(x(ixO^S, 2))**2 } )
    case default
      call mpistop("Sorry, coordinate unknown in getting gammainv")
    end select
  end subroutine get_gammainvij_hat

  !> Calculate dervitives of a scalar q within ixL in direction idir
  ! 2pt stencil
  subroutine partial_d(q,ixI^L,ixO^L,x,idir,gradq)
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S)
    double precision, intent(inout) :: gradq(ixI^S)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer                         :: jxO^L, hxO^L

    hxO^L=ixO^L-kr(idir,^D);
    jxO^L=ixO^L+kr(idir,^D);

    select case(coordinate)
    case(cartesian)
    !  gradq(ixO^S)=0.5d0*(q(jxO^S)-q(hxO^S))/dxlevel(idir)
      gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,1)-x(hxO^S,1)))
        {^NOONED
      case(^Z)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,2)-x(hxO^S,2))
        }
        {^IFTHREED
      case(^PHI)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,3)-x(hxO^S,3))
        }
      end select
    case(cylindrical)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
    case default
      call mpistop('Unknown geometry')
    end select
  end subroutine partial_d

  ! 5 pts stencil
  subroutine partial_d_5pts(q,ixI^L,ixO^L,x,idir,gradq)
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S)
    double precision, intent(inout) :: gradq(ixI^S)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer                         :: jxO^L, hxO^L
    integer                         :: j2xO^L, h2xO^L

    hxO^L=ixO^L-kr(idir,^D);
    jxO^L=ixO^L+kr(idir,^D);
    h2xO^L=ixO^L-kr(idir,^D)-kr(idir,^D);
    j2xO^L=ixO^L+kr(idir,^D)+kr(idir,^D);
    select case(coordinate)
    case(Cartesian)
        gradq(ixO^S)=(-q(j2xO^S)+8.0d0*q(jxO^S)-8.0d0*q(hxO^S)+q(h2xO^S))/&
                     (6.0d0*(x(jxO^S,idir)-x(hxO^S,idir)))
!                     (12.0d0*dxlevel(idir))
    case(Cartesian_stretched)
        gradq(ixO^S)=(-q(j2xO^S)+8.0d0*q(jxO^S)-8.0d0*q(hxO^S)+q(h2xO^S))/&
                     (6.0d0*(x(jxO^S,idir)-x(hxO^S,idir)))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(-q(j2xO^S)+8.0d0*q(jxO^S)-8.0d0*q(hxO^S)+q(h2xO^S))/&
                     (6.0d0*(x(jxO^S,1)-x(hxO^S,1)))
        {^NOONED
      case(2)
        gradq(ixO^S)=(-q(j2xO^S)+8.0d0*q(jxO^S)-8.0d0*q(hxO^S)+q(h2xO^S))/&
                     (6.0d0*(x(jxO^S,2)-x(hxO^S,2)))
        }
        {^IFTHREED
      case(3)
        gradq(ixO^S)=(-q(j2xO^S)+8.0d0*q(jxO^S)-8.0d0*q(hxO^S)+q(h2xO^S))/&
                     (6.0d0*(x(jxO^S,3)-x(hxO^S,3)))
        }
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixO^S)=(-q(j2xO^S)+8.0d0*q(jxO^S)-8.0d0*q(hxO^S)+q(h2xO^S))/&
                     (6.0d0*(x(jxO^S,phi_)-x(hxO^S,phi_)))
      else
        gradq(ixO^S)=(-q(j2xO^S)+8.0d0*q(jxO^S)-8.0d0*q(hxO^S)+q(h2xO^S))/&
                     (6.0d0*(x(jxO^S,idir)-x(hxO^S,idir)))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

  end subroutine partial_d_5pts

  !> Calculate dervitives of gammaij_hat in all coordinate
  subroutine partial_d_gamma_ij_hat(x,ixI^L,ixO^L,d_gamma_ij_hat)
    ! as some numerical metric derivatives cannot calulate by complex_de,
    ! so we use a two pt stencil derivative for metric components
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(inout) :: d_gamma_ij_hat(ixI^S, 1:3, 1:3, 1:3)


    d_gamma_ij_hat = 0.0d0

    select case (coordinate)
    case (cartesian)
      d_gamma_ij_hat(ixO^S,:,:,:) = 0.0d0
    case (cylindrical)
      ! note: r_ = 1,z_ = 2, phi_ = 3.
      d_gamma_ij_hat(ixO^S,3,3,1) = 2.0d0 * x(ixO^S, 1)
    case (spherical)
      ! note: r_ = 1,theta_=2, phi_ = 3.
      d_gamma_ij_hat(ixO^S,2,2,1) = 2.0d0 * x(ixO^S, 1)

      d_gamma_ij_hat(ixO^S,3,3,1) = d_gamma_ij_hat(ixO^S,2,2,1) {^NOONED * dsin(x(ixO^S, 2))**2 }
      d_gamma_ij_hat(ixO^S,3,3,2) = x(ixO^S, 1)**2 * 2.0d0 {^NOONED * dcos(x(ixO^S, 2)) * dsin(x(ixO^S, 2)) }
    case default
      call mpistop("Sorry, coordinate unknown in partial_d_gamma_ij_hat")
    end select
  end subroutine partial_d_gamma_ij_hat

  !> get the 3-metric in flat space sqrt_gamma_hat
  subroutine get_sqrt_gamma_hat(x, ixI^L, ixO^L, sqrt_gamma)
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(out) :: sqrt_gamma(ixI^S)
    double precision, intent(in)  :: x(ixI^S, 1:^ND)

    select case (coordinate)
    case (cartesian)
      sqrt_gamma(ixO^S) = 1.0d0
    case (cylindrical)
      sqrt_gamma(ixO^S) = dabs( x(ixO^S, 1) )
    case (spherical)
      sqrt_gamma(ixO^S) = dabs( x(ixO^S, 1)**2 {^NOONED * dsin(x(ixO^S, 2)) } )
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_sqrt_gamma_hat

  !> calculate christoffel symbols of cells after calculating the surface area and volume
  subroutine get_christoffel(ixI^L,ixO^L,x,christoffel)
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(inout) :: christoffel(ixI^S,1:3,1:3,1:3)
    double precision :: drs(ixI^S), dx2(ixI^S), dx3(ixI^S)
    integer          :: h1x^L{^NOONED, h2x^L}


    christoffel=0.0d0
    select case (coordinate)
    case (cartesian)
       return ! nothing to do here
    case (spherical)
      h1x^L=ixI^L-kr(1,^D); {^NOONED h2x^L=ixI^L-kr(2,^D);}
      drs(ixO^S)=mygeo%dx(ixO^S,1)
      {^NOONED
      dx2(ixO^S)=mygeo%dx(ixO^S,2)
      }
      {^IFTHREED
      dx3(ixO^S)=mygeo%dx(ixO^S,3)
      }

      christoffel(ixO^S,2,1,2)=0.5d0 * (mygeo%surfaceC1(ixO^S) - mygeo%surfaceC1(h1x^S)) / mygeo%dvolume(ixO^S)
      christoffel(ixO^S,2,2,1)=christoffel(ixO^S,2,1,2)
      christoffel(ixO^S,3,3,1)=christoffel(ixO^S,2,1,2)
      christoffel(ixO^S,3,1,3)=christoffel(ixO^S,2,1,2)

      christoffel(ixO^S,1,2,2)=-0.25d0 * ((x(ixO^S,1)+0.5d0*drs(ixO^S))**2*mygeo%surfaceC1(ixO^S) &
                                          - (x(ixO^S,1)-0.5d0*drs(ixO^S))**2*mygeo%surfaceC1(h1x^S)) / mygeo%dvolume(ixO^S)

      christoffel(ixO^S,1,3,3)=-0.25d0 / mygeo%dvolume(ixO^S) &
                           * ((x(ixO^S,1)+0.5d0*drs(ixO^S))**4 - (x(ixO^S,1)-0.5d0*drs(ixO^S))**4){^NOONED &
           *( (dcos(x(ixO^S,2)+0.5d0*dx2(ixO^S))**3/3.0d0 - dcos(x(ixO^S,2)+0.5d0*dx2(ixO^S))) &
             -(dcos(x(ixO^S,2)-0.5d0*dx2(ixO^S))**3/3.0d0 - dcos(x(ixO^S,2)-0.5d0*dx2(ixO^S))) )}{^IFTHREED*dx3(ixO^S)}

      {^NOONED
      christoffel(ixO^S,2,3,3)= - ( &
            (dsin(x(ixO^S,2)+0.5d0*dx2(ixO^S))**2)*mygeo%surfaceC2(ixO^S) &
          - (dsin(x(ixO^S,2)-0.5d0*dx2(ixO^S))**2)*mygeo%surfaceC2(h2x^S) &
              ) / 3.0d0 / mygeo%dvolume(ixO^S)

      christoffel(ixO^S,3,2,3)=dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
      christoffel(ixO^S,3,3,2)=christoffel(ixO^S,3,2,3)
      }
    case (cylindrical)
      call mpistop("Sorry, cylindertrial, will add later")
   !   drs(ixO^S)=mygeo%dx(ixO^S,1)
   !   {^NOONED
   !   dx2(ixO^S)=mygeo%dx(ixO^S,2)}
   !   {^IFTHREED
   !   dx3(ixO^S)=mygeo%dx(ixO^S,3)}

   !   if ( z_ == 2 ) then
   !      christoffel(ixO^S,1,3,3)= - (x(ixO^S,1)**2 +drs(ixO^S)**2/12.0d0) / x(ixO^S,1)
   !      christoffel(ixO^S,3,1,3)= 1.0d0 / x(ixO^S,1)
   !      christoffel(ixO^S,3,3,1) = christoffel(ixO^S,3,1,3)
   !   else
   !      christoffel(ixO^S,1,2,2)= - (x(ixO^S,1)**2 +drs(ixO^S)**2/12.0d0) / x(ixO^S,1)
   !      christoffel(ixO^S,2,1,2)= 1.0d0 / x(ixO^S,1)
   !      christoffel(ixO^S,2,2,1) = christoffel(ixO^S,2,1,2)
   !   end if
    case default
      call mpistop("Sorry, coordinate unknown")
    end select

  end subroutine get_christoffel

  !> transform vectors between natural and orthonomal basis
  subroutine get_natural2orthonormal(x, ixI^L, ixO^L, N2R)
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(out)   :: N2R(ixI^S, 1:3)

    N2R(ixO^S,1:3) = 1.0d0

    select case (coordinate)
    case (Cartesian)
      ! nothing to do here
    case (cylindrical)
      !if (^Z == 2) then
         N2R(ixO^S,3) = x(ixO^S, 1)
      !else
      !   N2R(ixO^S,2) = x(ixO^S, 1)
      !end if
    case (spherical)
      ! note: r_ = 1,theta_=2, phi_ = 3.
      N2R(ixO^S,2) = x(ixO^S, 1)
      N2R(ixO^S,3) = x(ixO^S, 1) {^NOONED * dsin(x(ixO^S, 2)) }
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_natural2orthonormal




!=============================================================================
