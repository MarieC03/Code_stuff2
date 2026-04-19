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
subroutine fill_gridvars(mygridvars,mypw)

use mod_physaux
!use mod_gridvars
use mod_metric
use constants
include 'amrvacdef.f'

type(walloc), dimension(ngridshi), intent(in)   :: mypw
type(walloc), dimension(ngridshi), intent(inout)  :: mygridvars
integer                                   :: igrid, iigrid, idir
double precision, dimension(ixG^T,1:ndir) :: beta, vu, vd, bd
double precision, dimension(ixG^T,1:nw)   :: w  ! will be primitive now
double precision, dimension(ixG^T,1:ndim)   :: x
!-----------------------------------------------------------------------------
call mpistop("stop this subroutinee  fillgridvar!")

do iigrid=1,igridstail; igrid=igrids(iigrid);
   call set_tmpGlobals(igrid)

   x(ixG^T,1:ndim) = px(igrid)%x(ixG^T,1:ndim)

   w(ixG^T,1:nw)   = mypw(igrid)%w(ixG^T,1:nw)
   call primitive(ixG^LL,ixG^LL,w,x)

!   mygridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0


   
   
 {#IFDEF PARTICLES_ADVECT
   ! Code units fields
   ! fill with Transport-velocity:
   do idir = 1, ndir
      call getv(w,x,ixG^LL,ixG^LL,idir,mygridvars(igrid)%w(ixG^T,vp1_-1+idir))
      mygridvars(igrid)%w(ixG^T,vp1_-1+idir) = mygridvars(igrid)%w(ixG^T,vp1_-1+idir)
   end do
   ! Fill density:
   mygridvars(igrid)%w(ixG^T,rhop_) = w(ixG^T,rho_)
   ! Fill pressure:
   call Pressuren(w,ixG^LL,ixG^LL,.false.,mygridvars(igrid)%w(ixG^T,pressp_),patchfalse)
   mygridvars(igrid)%w(ixG^T,pressp_) = mygridvars(igrid)%w(ixG^T,pressp_)
   ! Fill b:
   call get_b2(ixG^LL,ixG^LL,w,x,mygridvars(igrid)%w(ixG^T,bb2p_),w_is_primitive=.true.)
   mygridvars(igrid)%w(ixG^T,bb2p_) = mygridvars(igrid)%w(ixG^T,bb2p_)
}




{#IFDEF PARTICLES_LORENTZ
! Only Cartesian case.

   ! Scaled fields.
   mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = w(ixG^T,b1_:b^NC_) &
        * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY)

   do idir = 1, ndir
      beta(ixG^T,idir) =  w(ixG^T,u0_+idir)/w(ixG^T,lfac_) * UNIT_VELOCITY/const_c
   end do

   mygridvars(igrid)%w(ixG^T,ep1_) = mygridvars(igrid)%w(ixG^T,bp2_) * beta(ixG^T,3) &
        - mygridvars(igrid)%w(ixG^T,bp3_) * beta(ixG^T,2)
   
   mygridvars(igrid)%w(ixG^T,ep2_) = mygridvars(igrid)%w(ixG^T,bp3_) * beta(ixG^T,1) &
        - mygridvars(igrid)%w(ixG^T,bp1_) * beta(ixG^T,3)
   
   mygridvars(igrid)%w(ixG^T,ep3_) = mygridvars(igrid)%w(ixG^T,bp1_) * beta(ixG^T,2) &
        - mygridvars(igrid)%w(ixG^T,bp2_) * beta(ixG^T,1)
}

   


{#IFDEF PARTICLES_GEODESIC
mygridvars(igrid)%w(ixG^T,emiss_) = w(ixG^T,rho_)**2
}




{#IFDEF PARTICLES_GRLORENTZ
! fill with magnetic field:
   mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = w(ixG^T,b1_:b^NC_)

   do idir = 1, ndir
      vu(ixG^T,idir) = w(ixG^T,u0_+idir)/w(ixG^T,lfac_)
   end do

   call lower3(ixG^LL,ixG^LL,myM,vu,vd)
   call lower3(ixG^LL,ixG^LL,myM,w(ixG^T,b1_:b^NC_),bd)

   do idir = 1, ndir
      call acrossbU(ixG^LL,ixG^LL,myM,idir,bd,vd,mygridvars(igrid)%w(ixG^T,ep1_+idir-1))
   end do
}

end do

end subroutine fill_gridvars
!=============================================================================
