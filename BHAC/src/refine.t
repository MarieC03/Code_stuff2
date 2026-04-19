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
subroutine refine_grid(child_igrid,child_ipe,igrid,ipe,active)

include 'amrvacdef.f'

integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe
logical, intent(in) :: active

integer :: ic^D
!-----------------------------------------------------------------------------

! allocate solution space for new children
{do ic^DB=1,2\}
   call alloc_node(child_igrid(ic^D))
{end do\}

if ((time_advance .and. active).or.convert.or.firstprocess) then
   ! prolong igrid to new children
   call prolong_grid(child_igrid,child_ipe,igrid,ipe)
else
   ! Fill new created children with initial condition
   {do ic^DB=1,2\}
      call set_tmpGlobals(child_igrid(ic^D))
      call initial_condition(child_igrid(ic^D))
   {end do\}
end if

! remove solution space of original igrid, becoz now we use child_grid
call dealloc_node(igrid)
end subroutine refine_grid
!=============================================================================
subroutine prolong_grid(child_igrid,child_ipe,igrid,ipe)
{#IFDEF STAGGERED
use mod_amr_fct
}
include 'amrvacdef.f'

integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe

integer :: ix^L, ichild, ixCo^L, ic^D
double precision :: dxCo^D, xComin^D, dxFi^D, xFimin^D
!-----------------------------------------------------------------------------


{#IFNDEF DY_SP
if (covariant) myM=>pgeo(igrid)%m
}

if (typegridfill=="linear") then
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

{#IFNDEF DY_SP
! code test erase first
   if (amrentropy) then
      ix^L=ixM^LL^LADD1;
      call e_to_rhos(ixG^LL,ix^L,pw(igrid)%w,px(igrid)%x)
   else if (prolongprimitive) then
      ix^L=ixM^LL^LADD1;
      call primitive(ixG^LL,ix^L,pw(igrid)%w,px(igrid)%x)
   end if
}


   xComin^D=rnode(rpxmin^D_,igrid)\
   dxCo^D=rnode(rpdx^D_,igrid)\
end if


{#IFDEF STAGGERED
call old_neighbors(child_igrid,child_ipe,igrid,ipe)
}
{do ic^DB=1,2\}
   ichild=child_igrid(ic^D)


   {#IFNDEF DY_SP
   if (covariant) myM=>pgeo(ichild)%m
   }
   
   ixComin^D=ixMlo^D+(ic^D-1)*(ixMhi^D-ixMlo^D+1)/2\
   ixComax^D=ixMhi^D+(ic^D-2)*(ixMhi^D-ixMlo^D+1)/2\

   if (typegridfill=="linear") then
      xFimin^D=rnode(rpxmin^D_,ichild)\
      dxFi^D=rnode(rpdx^D_,ichild)\

      call prolong_2nd(ps(igrid),px(igrid)%x,ixCo^L,ps(ichild),px(ichild)%x, &
                   dxCo^D,xComin^D,dxFi^D,xFimin^D,ichild)
   else
      call prolong_1st(ps(igrid),ixCo^L,ps(ichild),px(ichild)%x)
   end if
{end do\}
   
{#IFNDEF DY_SP
if (covariant) myM=>pgeo(igrid)%m
}
! code test erase first
{#IFNDEF DY_SP
if (typegridfill=="linear") then
   if (amrentropy) then
      call rhos_to_e(ixG^LL,ix^L,pw(igrid)%w,px(igrid)%x)
   else if (prolongprimitive) then
      call conserve(ixG^LL,ix^L,pw(igrid)%w,px(igrid)%x,patchfalse)
   end if
end if
}

end subroutine prolong_grid
!=============================================================================
subroutine prolong_2nd(sCo,xCo,ixCo^L,sFi,xFi,dxCo^D,xComin^D,dxFi^D,xFimin^D,igridFi)
  {#IFDEF STAGGERED
  use mod_amr_fct
  }
use mod_imhd_intermediate, only: imhd_handle_small_values
include 'amrvacdef.f'

integer, intent(in) :: ixCo^L, igridFi
double precision, intent(in) :: dxCo^D, xComin^D, dxFi^D, xFimin^D
double precision, intent(in) :: xCo(ixG^T,1:ndim), xFi(ixG^T,1:ndim)
type(state), intent(in)      :: sCo
type(state), intent(inout)   :: sFi

integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, idim, iw
integer :: ixFi^L
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nw,ndim)
double precision :: xCo^D, xFi^D, eta^D, invdxCo^D
{#IFDEF STAGGERED
logical :: fine_^L
}
!-----------------------------------------------------------------------------
associate(wCo=>sCo%w%w, wFi=>sFi%w%w)
invdxCo^D=1.d0/dxCo^D;



{do ixCo^DB = ixCo^LIM^DB
   ! cell-centered coordinates of coarse grid point
   xCo^DB=xComin^DB+(dble(ixCo^DB-dixB)-half)*dxCo^DB

   ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixMlo^DB\}

   do idim=1,ndim
      hxCo^D=ixCo^D-kr(^D,idim)\
      jxCo^D=ixCo^D+kr(^D,idim)\

      do iw=1,nw
         slopeL=wCo(ixCo^D,iw)-wCo(hxCo^D,iw)
         slopeR=wCo(jxCo^D,iw)-wCo(ixCo^D,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idim)=signR*max(zero,min(dabs(slopeR), &
                                             signR*slopeL))
         case('woodward')
           slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), &
                              signR*slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idim)=signR*max(zero,min(mcbeta*dabs(slopeR), &
                              mcbeta*signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idim)=signR*max(zero,min(two*signR*slopeL, &
            (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idim)=signC*max(zero,min(dabs(slopeC), &
                             signC*slopeL,signC*slopeR))
         end select
      end do
   end do
   {do ix^DB=ixFi^DB,ixFi^DB+1
      ! cell-centered coordinates of fine grid point
      xFi^DB=xFimin^DB+(dble(ix^DB-dixB)-half)*dxFi^DB\}

      ! normalized distance between fine/coarse cell center
      ! in coarse cell: ranges from -0.5 to 0.5 in each direction
      ! (origin is coarse cell center)
      if (slab) then
         eta^D=(xFi^D-xCo^D)*invdxCo^D;
      else
         {eta^D=(xFi^D-xCo^D)*invdxCo^D &
               *two*(one-pgeo(igridFi)%dvolume(ix^DD) &
               /sum(pgeo(igridFi)%dvolume(ixFi^D:ixFi^D+1^D%ix^DD))) \}
      end if

      wFi(ix^D,1:nw) = wCo(ixCo^D,1:nw) &
                            + {(slope(1:nw,^D)*eta^D)+}
   {end do\}
{end do\}
{#IFDEF STAGGERED
call already_fine(sFi,igridFi,fine_^L)
call prolong_2nd_stg(sCo,sFi,ixCo^L,ixM^LL,dxCo^D,xComin^D,dxFi^D,xFimin^D,.false.,fine_^L)
}
! code test erase first
{#IFNDEF DY_SP
if (amrentropy) then
   call rhos_to_e(ixG^LL,ixM^LL,wFi,xFi)
else if (prolongprimitive) then
   call conserve(ixG^LL,ixM^LL,wFi,xFi,patchfalse)
end if
}

!call imhd_handle_small_values(wFi(ixG^T,1:nw),xFi(ixG^T,1:ndim),ixG^LL,ixM^LL, .True.)

end associate
end subroutine prolong_2nd
!=============================================================================
subroutine prolong_1st(sCo,ixCo^L,sFi,xFi)

include 'amrvacdef.f'

integer, intent(in)          :: ixCo^L
double precision, intent(in) :: xFi(ixG^T,1:ndim)
type(state), intent(in)      :: sCo
type(state), intent(out)     :: sFi
! .. local ..
integer                      :: ixCo^D, ixFi^D, iw
integer                      :: ixFi^L
!-----------------------------------------------------------------------------
associate(wCo=>sCo%w%w, wFi=>sFi%w%w)
  
{do ixCo^DB = ixCo^LIM^DB
   ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixMlo^DB\}
   forall(iw=1:nw) wFi(ixFi^D:ixFi^D+1,iw)=wCo(ixCo^D,iw)
{end do\}

end associate
end subroutine prolong_1st
!=============================================================================
