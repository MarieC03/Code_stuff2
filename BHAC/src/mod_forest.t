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

module mod_forest
   implicit none

   type tree_node_ptr
      type(tree_node), pointer :: node
   end type tree_node_ptr

   type tree_node
      integer :: ig^D, level, igrid, ipe, id
      logical :: leaf, active
      type(tree_node_ptr) :: parent, child(2^D&), neighbor(2,^ND), next, prev
   end type tree_node

   type(tree_node_ptr), dimension(:^D&), allocatable, save :: tree_root
   type(tree_node_ptr), dimension(:,:), allocatable, save :: igrid_to_node
   type(tree_node_ptr), dimension(:), allocatable, save :: level_head, &
                                                           level_tail
   integer, dimension(:,:), allocatable, save :: sfc

      !> Space filling curve for level 1 grid. sfc_iglevel1(^D, MN) gives ig^D (the
   !> spatial index of the grid) 
   integer, dimension(:,:), allocatable, save :: sfc_iglevel1

   !> iglevel1_sfc(ig^D) gives the Morton number for grid ig^D
   integer, dimension(:^D&), allocatable, save :: iglevel1_sfc
   
   integer, dimension(:), allocatable, save :: sfc_to_igrid
   integer, dimension(:), allocatable, save :: igrid_to_sfc
   integer, dimension(:), allocatable, save :: Morton_start, Morton_stop

   integer, dimension(:), allocatable, save :: Morton_sub_start, Morton_sub_stop

   logical, dimension(:,:), allocatable, save :: coarsen, refine, buffer, igrid_inuse

   !> Number of parent blocks
   integer, save :: nparents
   
   integer, save :: nleafs, nleafs_active
   integer       :: nglev1
   integer, dimension(:), allocatable, save :: nleafs_level

end module mod_forest
