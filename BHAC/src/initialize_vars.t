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
subroutine initialize_vars
  use mod_forest
  use mod_metric
include 'amrvacdef.f'
integer :: igrid, level, ipe, ig^D, i, j, k
logical :: ok
!-----------------------------------------------------------------------------

! Initialize Kronecker delta, and Levi-Civita tensor
do i=0,3
   do j=0,3
      if(i==j)then
         kr(i,j)=1
      else
         kr(i,j)=0
      endif
      if (i .gt. 0 .and. j .gt. 0) then
         do k=1,3
            if(i==j.or.j==k.or.k==i)then
               lvc(i,j,k)=0
            else if(i+1==j.or.i-2==j)then
               lvc(i,j,k)=1
            else
               lvc(i,j,k)=-1
            endif
         enddo
      endif
   enddo
enddo

! set time, time counter
if(.not.treset)t=zero
if(.not.itreset)it=0
dt=zero
dtimpl=zero
itmin=0

if(.not.time_accurate.or.residmin>smalldouble) then
  residual=one
endif 

! set all dt to zero
dt_grid(1:ngridshi)=zero

! check resolution
if ({mod(ixGhi^D,2)/=0|.or.}) then
   call mpistop("mesh widths must give even number grid points")
end if
ixM^LL=ixG^LL^LSUBdixB;
if (errorestimate==1) then
   if ({mod(ixMhi^D-ixMlo^D+1,4)/=0|.or.}) then
      call mpistop("mesh widths must be divisible by 4 for Richardson")
   end if
end if

if (nbufferx^D>(ixMhi^D-ixMlo^D+1)|.or.) then
   write(unitterm,*) 'nbufferx^D bigger than mesh size makes no sense.'
   write(unitterm,*) 'Decrease nbufferx or increase mesh size'
   call mpistop("")
end if

! initialize dx arrays on finer (>1) levels
do level=2,mxnest
   {dx(^D,level) = dx(^D,level-1) * half\}  ! refine ratio 2
end do

! domain decomposition
! physical extent of a grid block at level 1, per dimension
^D&dg^D(1)=dx(^D,1)*dble(ixGhi^D-2*dixB)\

! number of grid blocks at level 1 in simulation domain, per dimension
^D&ng^D(1)=nint((xprobmax^D-xprobmin^D)/dg^D(1))\

! total number of grid blocks at level 1
nglev1={ng^D(1)*}

do level=2,mxnest
   dg^D(level)=half*dg^D(level-1);
   ng^D(level)=ng^D(level-1)*2;
end do

! check that specified stepsize correctly divides domain

!as changed smalldouble to be 1.0d-16
ok=({(abs(dble(ng^D(1))*dg^D(1)-(xprobmax^D-xprobmin^D))<=1.0d-12)|.and.})
!ok=({(abs(dble(ng^D(1))*dg^D(1)-(xprobmax^D-xprobmin^D))<=smalldouble)|.and.})

if (.not.ok) then
   write(unitterm,*)"domain cannot be divided by meshes of given gridsize"
   call mpistop("domain cannot be divided by meshes of given gridsize")
end if


poleB=.false.
if (.not.slab) call set_pole

do igrid=1,ngridshi
! All nullification already at the declaration stage.   
! Associate the state structure with the pw arrays.  

   psold(igrid) = state(igrid=igrid,w=pwold(igrid),x=px(igrid),geo=pgeo(igrid) &
        {#IFDEF STAGGERED , ws=pwsold(igrid)} &
        )

   ps(igrid)        = state(igrid=igrid,w=pw(igrid),x=px(igrid),geo=pgeo(igrid) &
        {#IFDEF STAGGERED , ws=pws(igrid)} &
        )

   ps1(igrid)       = state(igrid=igrid,w=pw1(igrid),x=px(igrid),geo=pgeo(igrid) &
        {#IFDEF STAGGERED , ws=pws1(igrid)} &
        )

   psCoarse(igrid)  = state(igrid=igrid,w=pwCoarse(igrid),x=pxCoarse(igrid),geo=pgeoCoarse(igrid), &
        is_coarse=.true. &
        {#IFDEF STAGGERED , ws=pwsCoarse(igrid)} &
        )
   if (errorestimate==1) then
     psCoCo(igrid)    = state(igrid=igrid,w=pwCoCo(igrid),geo=pgeoCoCo(igrid), &
          is_coarse=.true. &
          {#IFDEF STAGGERED , ws=pwsCoCo(igrid)} &
          )
   endif

   {#IFDEF M1
   psm1(igrid)       = m1_impl_state(igrid=igrid,pwrad=pwm1(igrid))
   }

   if (use_gw_br .and. gw_br_use_I3_ij) then
     ! method I does not need this
     ps_t_dot_old(igrid)  = state(igrid=igrid,w=pw_t_dot_old(igrid))
     ps_t_dot_old2(igrid) = state(igrid=igrid,w=pw_t_dot_old2(igrid))
   endif


   select case(typeadvance)
   case("threestep","fourstep","ImEx122","ImEx222","FullyImplicit")
      ps2(igrid)    = state(igrid=igrid,w=pw2(igrid),x=px(igrid),geo=pgeo(igrid) &
        {#IFDEF STAGGERED , ws=pws2(igrid)} &
        )
   case("ssprk43","rk4")
      ps2(igrid)    = state(igrid=igrid,w=pw2(igrid),x=px(igrid),geo=pgeo(igrid) &
        {#IFDEF STAGGERED , ws=pws2(igrid)} &
        )
      ps3(igrid)    = state(igrid=igrid,w=pw3(igrid),x=px(igrid),geo=pgeo(igrid) &
        {#IFDEF STAGGERED , ws=pws3(igrid)} &
        )
   case("ssprk54")
      ps4(igrid)    = state(igrid=igrid,w=pw4(igrid),x=px(igrid),geo=pgeo(igrid) &
        {#IFDEF STAGGERED , ws=pws4(igrid)} &
        )
   case default
   end select

!   if (nstep>2) then
!      ps2(igrid)    = state(igrid=igrid,w=pw2(igrid),x=px(igrid),geo=pgeo(igrid) &
!        {#IFDEF STAGGERED , ws=pws2(igrid)} &
!        )
!   end if
!
!   if (nstep>3) then
!      ps3(igrid)    = state(igrid=igrid,w=pw3(igrid),x=px(igrid),geo=pgeo(igrid) &
!        {#IFDEF STAGGERED , ws=pws3(igrid)} &
!        )
!   end if
!
!   if (nstep>4) then
!      ps4(igrid)    = state(igrid=igrid,w=pw4(igrid),x=px(igrid),geo=pgeo(igrid) &
!        {#IFDEF STAGGERED , ws=pws4(igrid)} &
!        )
!   end if

   if (residmin>smalldouble) then
      psres(igrid)  = state(igrid=igrid,w=pwres(igrid),x=px(igrid),geo=pgeo(igrid) &
        {#IFDEF STAGGERED , ws=pwsres(igrid)} &
        )
   end if

end do

! on each processor, create for later use a default patch array
allocate(patchfalse(ixG^T))
patchfalse(ixG^T)=.false.

! initialize connectivity data
igridstail=0

! allocate memory for forest data structures
allocate(level_head(mxnest),level_tail(mxnest))
do level=1,mxnest
   nullify(level_head(level)%node,level_tail(level)%node)
end do

allocate(igrid_to_node(ngridshi,0:npe-1))
do ipe=0,npe-1
   do igrid=1,ngridshi
      nullify(igrid_to_node(igrid,ipe)%node)
   end do
end do

allocate(sfc(1:3,ngridshi*npe))

allocate(igrid_to_sfc(ngridshi))

sfc=0
allocate(Morton_start(0:npe-1),Morton_stop(0:npe-1))
allocate(Morton_sub_start(0:npe-1),Morton_sub_stop(0:npe-1))

allocate(nleafs_level(1:nlevelshi))

allocate(coarsen(ngridshi,0:npe-1),refine(ngridshi,0:npe-1))
coarsen=.false.
refine=.false.
if (nbufferx^D/=0|.or.) then
   allocate(buffer(ngridshi,0:npe-1))
   buffer=.false.
end if
allocate(igrid_inuse(ngridshi,0:npe-1))
igrid_inuse=.false.

allocate(tree_root(1:ng^D(1)))
{do ig^DB=1,ng^DB(1)\}
   nullify(tree_root(ig^D)%node)
{end do\}


! default the physical scaling parameters:
UNIT_LENGTH   = ONE
UNIT_DENSITY  = ONE
UNIT_VELOCITY = ONE

! initialize the meta-info for the metric datastructure:
if (covariant) then
   call init_metric
end if


end subroutine initialize_vars
!=============================================================================
subroutine set_tmpGlobals(igrid)

  include 'amrvacdef.f'

  integer, intent(in)              :: igrid
  !-----------------------------------------------------------------------------

  saveigrid       = igrid
  typelimiter     = type_limiter(node(plevel_,igrid))
  typegradlimiter = type_gradient_limiter(node(plevel_,igrid))

  ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  
  if (.not.slab) mygeo => pgeo(igrid)

  {#IFNDEF DY_SP
    if(covariant)myM => mygeo%m
  }
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     myB0      => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if

end subroutine set_tmpGlobals
!=============================================================================
