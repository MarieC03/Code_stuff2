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

module mod_comm_gc
   use mod_indices, only: ngridshi
   use mod_physicaldata

   ! === Pointers to the blocks and auxiliary buffers ===

   type(state), pointer, dimension(:), save :: pstate, pstateCo
   type(walloc), save :: pole_buf, pole_buf_stg !! Buffers to reverse phi at poles

   ! === Extents of blocks and other useful integers ===
   integer, save :: ixCoG^L, ixCoM^L
   integer, save :: nx^D, nxCo^D
   integer, save :: dixBCo, interpolation_order

   ! === Time for boundary conditions ===
   double precision, save :: time

   ! === Ranges and sizes of communications ===

   ! same resolution level
   integer, save :: nrecv_gc_srl, nsend_gc_srl
   integer, save :: nbuff_gc_recv_srl, nbuff_gc_send_srl
   integer, save :: ibuf_send_srl, ibuf_recv_srl
   integer, save :: isend_srl, irecv_srl

   integer, dimension(-1:1), save :: ixS_srl_^L, ixR_srl_^L
   integer, dimension(^ND,-1:1), save :: ixS_srl_stg_^L, ixR_srl_stg_^L

   integer, dimension(-1:1^D&), save     :: sizes_srl_send, sizes_srl_recv
   integer, dimension(^ND,-1:1^D&), save :: sizes_srl_send_stg, sizes_srl_recv_stg
   integer, dimension(-1:1^D&), save     :: sizes_srl_send_total, sizes_srl_recv_total

   integer, dimension(:), allocatable, save :: recvrequest_srl, sendrequest_srl
   integer, dimension(:,:), allocatable, save :: recvstatus_srl, sendstatus_srl
 
   double precision, dimension(:), allocatable, asynchronous, save :: recvbuffer_srl, sendbuffer_srl

   ! restriction
   integer, save :: nrecv_gc_r, nsend_gc_r
   integer, save :: nbuff_gc_recv_r, nbuff_gc_send_r
   integer, save :: ibuf_send_r, ibuf_recv_r
   integer, save :: isend_r, irecv_r

   integer, dimension(-1:1), save :: ixS_r_^L
   integer, dimension(^ND,-1:1), save :: ixS_r_stg_^L
   integer, dimension(0:3), save  :: ixR_r_^L
   integer, dimension(^ND,0:3), save  :: ixR_r_stg_^L

   integer, dimension(-1:1^D&), save     :: sizes_r_send
   integer, dimension(0:3^D&), save      :: sizes_r_recv
   integer, dimension(^ND,-1:1^D&), save :: sizes_r_send_stg
   integer, dimension(^ND,0:3^D&), save  :: sizes_r_recv_stg
   integer, dimension(-1:1^D&), save     :: sizes_r_send_total
   integer, dimension(0:3^D&), save      :: sizes_r_recv_total

   integer, dimension(:), allocatable, save :: recvrequest_r, sendrequest_r
   integer, dimension(:,:), allocatable, save :: recvstatus_r, sendstatus_r

   double precision, dimension(:), allocatable, asynchronous, save :: recvbuffer_r, sendbuffer_r

   ! prolongation
   integer, save :: nrecv_gc_p, nsend_gc_p
   integer, save :: nbuff_gc_recv_p, nbuff_gc_send_p
   integer, save :: ibuf_send_p, ibuf_recv_p
   integer, save :: isend_p, irecv_p

   integer, dimension(0:3), save  :: ixS_p_^L, ixR_p_^L
   integer, dimension(^ND,0:3), save  :: ixS_p_stg_^L, ixR_p_stg_^L

   integer, dimension(0:3^D&), save     :: sizes_p_send, sizes_p_recv
   integer, dimension(^ND,0:3^D&), save :: sizes_p_send_stg, sizes_p_recv_stg
   integer, dimension(0:3^D&), save     :: sizes_p_send_total, sizes_p_recv_total

   integer, dimension(:), allocatable, save :: recvrequest_p, sendrequest_p
   integer, dimension(:,:), allocatable, save :: recvstatus_p, sendstatus_p

   double precision, dimension(:), allocatable, asynchronous, save :: recvbuffer_p, sendbuffer_p

end module mod_comm_gc

