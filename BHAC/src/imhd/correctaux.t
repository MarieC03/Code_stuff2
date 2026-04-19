!fixme   : use con2prim to find lfac, xi are weird!
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
subroutine getaux(clipping,w,x,ixI^L,ixO^L,subname)

! Calculate the auxiliary variables lfac and xi within ixO^L
! clipping is not used (yet) 
!use mod_con2prim, only: con2prim
use mod_metric, only: raise3, square3u
use mod_imhd_con2prim, only: imhd_con2prim
include 'amrvacdef.f'

logical, intent(in)            :: clipping
integer, intent(in)            :: ixI^L, ixO^L
double precision, intent(inout):: w(ixI^S,1:nw)
double precision, intent(in)   :: x(ixI^S,1:ndim)
character(len=*), intent(in)   :: subname
! .. local ..
integer::          err,ix^D, i^D
double precision :: {^C&sold^C_},{^C&bold^C_}, lfacold, xiold
integer          :: patchierror(ixG^T)
double precision, dimension(ixI^S) :: ssqr, bsqr, sdotb2
double precision, dimension(ixI^S,1:^NC) :: sU
logical        :: patchw(ixG^T)

!-----------------------------------------------------------------------------
    patchw = .false.    
!{do ix^DB= ixO^LIM^DB\}
!    if (tlow>zero) call fixp_usr(ixI^L,ix^D,ix^D,w,x)
!{enddo^D&\}


    ! 1. the error secuity done inside con2prim
    ! 2. only need new xi, lfac now, true the need_aux_only
!    call imhd_con2prim(ixI^L, ixO^L, x, w, patchw, .true.)
!    call imhd_con2prim(ixI^L, ixO^L, x, w, patchw, .false.)
end subroutine getaux
