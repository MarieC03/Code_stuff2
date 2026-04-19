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

module mod_fix_conserve
   use mod_indices, only: ngridshi
   implicit none

   type fluxalloc
      double precision, dimension(:^D&,:), pointer:: flux => null()  
      double precision, dimension(:^D&,:), pointer:: edge => null()
      !double precision, dimension(:^D&,:), pointer:: flux => null  
      !double precision, dimension(:^D&,:), pointer:: edge => null
   end type fluxalloc
   type(fluxalloc), dimension(2,^ND,ngridshi), save :: pflux

   integer, save :: nrecv, nsend
   double precision, dimension(:), allocatable, asynchronous, save :: recvbuffer, sendbuffer
   integer, dimension(^ND), save :: isize
   integer, save                 :: ibuf,ibuf_send ! todo: make buffer handling threadsafe
{#IFDEF STAGGERED
   integer, save :: nrecv_ct, nsend_ct
   double precision, dimension(:), allocatable, asynchronous, save :: recvbuffer_cc, sendbuffer_cc
   integer, dimension(^ND), save :: isize_stg
   integer, save                 :: ibuf_cc,ibuf_cc_send
}
end module mod_fix_conserve
