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
subroutine printlog_special

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

include 'amrvacdef.f'

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,nwmax,w,s,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L,nwmax
double precision                   :: w(ixI^S,1:nwmax)
type(state)                        :: s
double precision                   :: normconv(0:nwmax)
!-----------------------------------------------------------------------------
associate(x=>s%x%x{#IFDEF STAGGERED ,ws=>s%ws%w})


end associate
end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

call mpistop("special wnames and primnames undefined")

! Example : as above in specialvar_output, assuming relativistic HD here...
! primnames= TRIM(primnames)//' '//'-rho'
! wnames=TRIM(wnames)//' '//'-d'

end subroutine specialvarnames_output
!=============================================================================
{#IFDEF UCONVERT
subroutine userspecialconvert(qunitconvert)
! Allow user to use their own data-postprocessing procedures

include 'amrvacdef.f'
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------

end subroutine userspecialconvert
!=============================================================================
}
{#IFDEF TRANSFORMW
subroutine transformw_usr(w,wtf,eqpar_tf,ixI^L,ixO^L)
! regenerate w and eqpar arrays to output into *tf.dat, e.g., add/remove e_
! variable
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: wtf(ixI^S,1:nwtf)
double precision, intent(out):: eqpar_tf(neqpartf)

!-----------------------------------------------------------------------------

end subroutine transformw_usr
!=============================================================================
}
{#IFDEF SPECIALTOLERANCE
subroutine special_tolerance(xlocal,tolerance)
!PURPOSE: use different tolerance in special regions for AMR to
!reduce/increase resolution there where nothing/something interesting happens.
include 'amrvacdef.f'

double precision, intent(in) :: xlocal(1:ndim)
double precision, intent(inout) :: tolerance

double precision :: bczone^D,addtol,tol_add
!-----------------------------------------------------------------------------
!amplitude of additional tolerance
addtol=0.3d0
! thickness of near-boundary region
bczone1=0.2d0*(xprobmax1-xprobmin1)
! linear changing of additional tolerance
if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
  tol_add=(1.d0-min(xlocal(1)-xprobmin1,xprobmax1-xlocal(1))/bczone1)*addtol
endif
bczone2=0.2d0*(xprobmax2-xprobmin2)
if(xlocal(2)-xprobmin2 < bczone2 .or. xprobmax2-xlocal(2) < bczone2) then
  tol_add=(1.d0-min(xlocal(2)-xprobmin2,xprobmax2-xlocal(2))/bczone2)*addtol
endif
bczone3=0.2d0*(xprobmax3-xprobmin3)
if(xprobmax3-xlocal(3) < bczone3) then
  tol_add=(1.d0-(xprobmax3-xlocal(3))/bczone3)*addtol
endif
tolerance=tolerance+tol_add

end subroutine special_tolerance
!=============================================================================
}
