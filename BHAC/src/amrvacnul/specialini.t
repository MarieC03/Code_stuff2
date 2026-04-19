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
subroutine initglobaldata_usr

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

! initialize one grid within ixO^L

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical patchw(ixG^T)
patchw(ixO^S) = .false.
!-----------------------------------------------------------------------------

w(ixO^S,1:nw)=zero

w(ixO^S,rho_)=1.0d0

{#IFDEF ENERGY
   w(ixO^S,pp_) = 1.0d0
}

call conserve(ixI^L,ixO^L,w,x,patchw)
end subroutine initonegrid_usr
!=============================================================================
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()


include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixC^L,idir
double precision, intent(in)       :: xC(ixI^S,1:ndim)
double precision, intent(out)      :: A(ixI^S)

!-----------------------------------------------------------------------------

A(ixC^S) = zero


end subroutine initvecpot_usr
!=============================================================================
