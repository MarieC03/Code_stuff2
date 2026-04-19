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

module mod_connectivity
   use mod_indices, only: ngridshi
   implicit none

   integer, dimension(2,-1:1^D&,ngridshi), save :: neighbor
   integer, dimension(2,0:3^D&,ngridshi), save :: neighbor_child
   integer, dimension(-1:1^D&,ngridshi), save :: neighbor_type
   logical, dimension(-1:1^D&,ngridshi), save :: neighbor_active
   {^IFPHI integer, dimension(-1:1^D&,ngridshi), save :: neighbor_pole}

   integer, dimension(ngridshi), save :: igrids
   integer, save :: igridstail

   integer, dimension(ngridshi), save :: igrids_active
   integer, save :: igridstail_active
   integer, dimension(ngridshi), save :: igrids_passive
   integer, save :: igridstail_passive

   integer, dimension(^ND), save :: nrecv_fc, nsend_fc
{#IFDEF STAGGERED
   integer, dimension(^ND), save :: nrecv_cc, nsend_cc
}
   integer, save :: nrecv_bc_srl_13, nsend_bc_srl_13, &
                    nrecv_bc_srl_2, nsend_bc_srl_2, &
                    nrecv_bc_r, nsend_bc_r, &
                    nrecv_bc_r_13, nsend_bc_r_13, &
                    nrecv_bc_r_2, nsend_bc_r_2, &
                    nrecv_bc_p, nsend_bc_p

end module mod_connectivity
