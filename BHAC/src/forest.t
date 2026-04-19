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
!{#IFNDEF SFC_ON_BASE
!subroutine init_forest_root
!use mod_forest
!include 'amrvacdef.f'
!
!integer :: ig^D, level, igrid, ipe
!integer :: iside, i^D, Morton_no
!
!integer, external :: getnode
!!-----------------------------------------------------------------------------
!level=1
!Morton_no=0
!ipe=0
!nleafs={ng^D(1)*}
!nleafs_active=nleafs
!nleafs_level(1)={ng^D(1)*}
!nleafs_level(2:nlevelshi)=0
!call get_Morton_range
!{do ig^DB=1,ng^DB(1)\}
!   Morton_no=Morton_no+1
!   if (Morton_no>Morton_stop(ipe)) ipe=ipe+1
!   igrid=getnode(ipe)
!   if (ipe==mype) sfc_to_igrid(Morton_no)=igrid
!   call init_tree_leaf(tree_root(ig^D),ig^D,level,igrid,ipe,.true.)
!{end do\}
!
!! update root neighbor
!{do ig^DB=1,ng^DB(1)\}
!   {do iside=1,2
!      i^DD=kr(^DD,^D)*(2*iside-3);
!      call find_root_neighbor(tree_root(ig^DD)%node%neighbor(iside,^D), &
!                              tree_root(ig^DD),i^DD)
!   end do\}
!{end do\}
!
!end subroutine init_forest_root
!!=============================================================================
!}{#IFDEF SFC_ON_BASE
subroutine init_forest_root
use mod_forest
include 'amrvacdef.f'

integer :: ig^D, level, igrid, ipe
integer :: iside, i^D, Morton_no, isfc

integer, external :: getnode
!-----------------------------------------------------------------------------
level=1
Morton_no=0
ipe=0
nparents=0
nleafs={ng^D(1)*}
nleafs_active=nleafs
nleafs_level(1)={ng^D(1)*}
nleafs_level(2:nlevelshi)=0
call get_Morton_range
! Generate Morton-order space-filling curve to connect level 1 blocks
call level1_Morton_order
do isfc=1,nglev1
   ig^D=sfc_iglevel1(^D,isfc)\
   Morton_no=Morton_no+1
   if (Morton_no>Morton_stop(ipe)) ipe=ipe+1
   igrid=getnode(ipe)
   if (ipe==mype) then
      sfc_to_igrid(Morton_no)=igrid
      igrid_to_sfc(igrid)=Morton_no
   end if
   call init_tree_leaf(tree_root(ig^D),ig^D,level,igrid,ipe,.true.)
end do

! update root neighbor
{do ig^DB=1,ng^DB(1)\}
   {do iside=1,2
      i^DD=kr(^DD,^D)*(2*iside-3);
      call find_root_neighbor(tree_root(ig^DD)%node%neighbor(iside,^D), &
                              tree_root(ig^DD),i^DD)
   end do\}
{end do\}

! This call is here to ensure the sfc array is initialized
call amr_Morton_order()

end subroutine init_forest_root
!=============================================================================
subroutine level1_Morton_order
! use Morton curve to connect level 1 grid blocks
use mod_forest
include 'amrvacdef.f'

integer, allocatable :: gsq_sfc(:^D&)
integer :: ig^D, ngsq^D, Morton_no
integer(kind=8), external :: mortonEncode
!-----------------------------------------------------------------------------
! use the smallest square/cube to cover the full domain 
ngsq^D=2**ceiling(log(real(ng^D(1)))/log(2.0));
{^NOONED
{ngsq^D=max(ngsq^DD) \}
}
allocate(gsq_sfc(ngsq^D))
! get Morton-order numbers in the square/cube
{do ig^DB=1,ngsq^DB\}
   gsq_sfc(ig^D)=int(mortonEncode(ig^D-1,ndim))+1
{end do\}
! delete Morton blocks that are out of the domain
{do ig^DB=1,ngsq^DB\}
   if (ig^D>ng^D(1)|.or.) then
      where(gsq_sfc>=gsq_sfc(ig^D))
         gsq_sfc=gsq_sfc-1
      end where
   end if
{end do\}
! copy the modified Morton numbers to the blocks in the domain
allocate(iglevel1_sfc(ng^D(1)))
allocate(sfc_iglevel1(ndim,nglev1))
{do ig^DB=1,ng^DB(1)\}
   iglevel1_sfc(ig^D)=gsq_sfc(ig^D)
   {sfc_iglevel1(^D,iglevel1_sfc(ig^DD))=ig^D \}
{end do\}

end subroutine level1_Morton_order
!=============================================================================
integer(kind=8) function mortonEncode(ig^D,ndim)
use iso_fortran_env, only : int64
implicit none
integer(kind=4) :: i,ig^D,ndim
integer(kind=8) :: answer
answer = 0
do i=0,64/ndim
  {^IFONED answer=ig1}
  {^IFTWOD
   answer=ior(answer,ior(ishft(iand(ig1,ishft(1_int64,i)),i),&
          ishft(iand(ig2,ishft(1_int64,i)),i+1)))
  \}
  {^IFTHREED
   answer=ior(answer,ior(ishft(iand(ig1,ishft(1_int64,i)),2*i),ior(ishft(&
   iand(ig2,ishft(1_int64,i)),2*i+1),ishft(iand(ig3,ishft(1_int64,i)),2*i+2))))
  \}
end do
mortonEncode=answer
return
end function mortonEncode
!=============================================================================
!}
subroutine init_tree_leaf(tree,ig^D,level,igrid,ipe,active)
use mod_forest
implicit none

type(tree_node_ptr) :: tree
integer, intent(in) :: ig^D, level, igrid, ipe
logical, intent(in) :: active
integer :: ic^D
!-----------------------------------------------------------------------------
allocate(tree%node)

tree%node%ig^D=ig^D;
tree%node%level=level
tree%node%igrid=igrid
tree%node%ipe=ipe

tree%node%leaf=.true.
tree%node%active=active

nullify(tree%node%parent%node)
{do ic^DB=1,2\}
   nullify(tree%node%child(ic^D)%node)
{end do\}

call add_to_linked_list(level,tree)

! initialize neighbor pointers
nullify({tree%node%neighbor(1,^D)%node},{tree%node%neighbor(2,^D)%node})

igrid_to_node(igrid,ipe)%node => tree%node

end subroutine init_tree_leaf
!=============================================================================
subroutine coarsen_tree_leaf(igrid,ipe,child_igrid,child_ipe,active)
use mod_forest
implicit none

integer, intent(in) :: igrid, ipe
integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
logical, intent(out) :: active


integer :: level, ic^D, child_level, iside, iotherside, vote
type(tree_node_ptr) :: tree, child, child_neighbor
!-----------------------------------------------------------------------------
tree%node => igrid_to_node(child_igrid(1^D&),child_ipe(1^D&))%node%parent%node
level=tree%node%level

call add_to_linked_list(level,tree)

child_level=level+1
vote=0

{do ic^DB=1,2\}
   child%node => tree%node%child(ic^D)%node

!  vote for active:
   if(child%node%active) vote=vote+1

   call delete_from_linked_list(child_level,child)

   ! update neighbor pointers
   {iside=ic^D
   child_neighbor%node => child%node%neighbor(iside,^D)%node
   if (associated(child_neighbor%node)) then
      if (child%node%ig^D==child_neighbor%node%ig^D) then ! pole
         nullify(child_neighbor%node%neighbor(iside,^D)%node)
      else
         iotherside=3-iside
         nullify(child_neighbor%node%neighbor(iotherside,^D)%node)
      end if
   end if\}

   nullify(tree%node%child(ic^D)%node)
   deallocate(igrid_to_node(child_igrid(ic^D),child_ipe(ic^D))%node)
{end do\}

tree%node%leaf=.true.
tree%node%igrid=igrid
tree%node%ipe=ipe
igrid_to_node(igrid,ipe)%node => tree%node

!  Count the vote and set active/passive state:

if (vote /= 2**^ND) then 
!if (vote == 0) then 
   tree%node%active = .false.
   nleafs_active = nleafs_active - vote
else
   tree%node%active = .true.
   nleafs_active = nleafs_active - vote + 1
end if
active = tree%node%active

nleafs=nleafs-2**^ND+1
nparents=nparents-1
nleafs_level(child_level)=nleafs_level(child_level)-2**^ND
nleafs_level(level)=nleafs_level(level)+1

end subroutine coarsen_tree_leaf
!=============================================================================
subroutine refine_tree_leaf(child_igrid,child_ipe,igrid,ipe,active)
use mod_forest
include 'amrvacdef.f'

integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe
logical, intent(out):: active

integer :: ig^D, level, i^D, ic^D, child_ig^D, child_level, iside
integer :: my_neighbor_type
logical, dimension(ndim) :: pole
type(tree_node_ptr) :: tree, child, my_neighbor
!-----------------------------------------------------------------------------
tree%node => igrid_to_node(igrid,ipe)%node
ig^D=tree%node%ig^D;
level=tree%node%level
active=tree%node%active

tree%node%ipe=-1
tree%node%igrid=0
tree%node%leaf=.false.
tree%node%active=.true.

call delete_from_linked_list(level,tree)

child_level=level+1

{do ic^DB=1,2\}
   child_ig^D=2*(ig^D-1)+ic^D;
   call init_tree_leaf(child,child_ig^D,child_level, &
                             child_igrid(ic^D),child_ipe(ic^D),active)

   igrid_to_node(child_igrid(ic^D),child_ipe(ic^D))%node => child%node

   tree%node%child(ic^D)%node => child%node
   child%node%parent%node => tree%node
{end do\}

! update neighbor pointers
{do ic^DB=1,2\}
   child%node => tree%node%child(ic^D)%node
   {iside=ic^D
   i^DD=kr(^DD,^D)*(2*iside-3);
   call find_neighbor(my_neighbor,my_neighbor_type,child,i^DD,pole)
   select case (my_neighbor_type)
   case (3,4)
      child%node%neighbor(iside,^D)%node => my_neighbor%node
      if (pole(^D)) then
         my_neighbor%node%neighbor(iside,^D)%node => child%node
      else
         my_neighbor%node%neighbor(3-iside,^D)%node => child%node
      end if
   case default
      nullify(child%node%neighbor(iside,^D)%node)
   end select
   child%node%neighbor(3-ic^D,^D)%node=>tree%node%child(3-ic^D^D%ic^DD)%node\}
{end do\}

nleafs=nleafs+2**^ND-1
nparents=nparents+1
nleafs_level(child_level)=nleafs_level(child_level)+2**^ND
nleafs_level(level)=nleafs_level(level)-1

if (active) nleafs_active = nleafs_active + 2**^ND-1

end subroutine refine_tree_leaf
!=============================================================================
subroutine change_ipe_tree_leaf(recv_igrid,recv_ipe,send_igrid,send_ipe)
use mod_forest
implicit none

integer, intent(in) :: recv_igrid, recv_ipe, send_igrid, send_ipe

type(tree_node_ptr) :: tree
!-----------------------------------------------------------------------------
tree%node => igrid_to_node(send_igrid,send_ipe)%node

tree%node%igrid=recv_igrid
tree%node%ipe=recv_ipe

nullify(igrid_to_node(send_igrid,send_ipe)%node)
igrid_to_node(recv_igrid,recv_ipe)%node => tree%node

end subroutine change_ipe_tree_leaf
!=============================================================================
subroutine add_to_linked_list(level,tree)
use mod_forest
implicit none

integer, intent(in) :: level
type(tree_node_ptr) :: tree
!-----------------------------------------------------------------------------
nullify(tree%node%next%node)
if (associated(level_head(level)%node)) then
   tree%node%prev%node => level_tail(level)%node
   level_tail(level)%node%next%node => tree%node
   level_tail(level)%node => tree%node
else
   level_head(level)%node => tree%node
   level_tail(level)%node => tree%node
   nullify(tree%node%prev%node)
end if

end subroutine add_to_linked_list
!=============================================================================
subroutine delete_from_linked_list(level,tree)
use mod_forest
implicit none

integer, intent(in) :: level
type(tree_node_ptr) :: tree

type(tree_node_ptr) :: next, prev
!-----------------------------------------------------------------------------
prev%node => tree%node%prev%node
next%node => tree%node%next%node
if (associated(next%node).and.associated(prev%node)) then
   prev%node%next%node => next%node
   next%node%prev%node => prev%node
else if (associated(prev%node)) then
   level_tail(level)%node => prev%node
   nullify(prev%node%next%node)
else if (associated(next%node)) then
   level_head(level)%node => next%node
   nullify(next%node%prev%node)
else
   nullify(level_head(level)%node)
   nullify(level_tail(level)%node)
end if

end subroutine delete_from_linked_list
!=============================================================================
subroutine write_forest(file_handle, count_forest)
use mod_forest
include 'amrvacdef.f'

integer, intent(in) :: file_handle
integer, intent(out) :: count_forest

integer, dimension(MPI_STATUS_SIZE) :: status
integer :: ig^D, isfc
!-----------------------------------------------------------------------------
count_forest = 0
do isfc=1,nglev1
   ig^D=sfc_iglevel1(^D,isfc)\
   call write_node(tree_root(ig^D))
end do

!{#IFNDEF SFC_ON_BASE
!{do ig^DB=1,ng^DB(1)\}
!   call write_node(tree_root(ig^D))
!{end do\}
!}{#IFDEF SFC_ON_BASE
!do isfc=1,nglev1
!   ig^D=sfc_iglevel1(^D,isfc)\
!   call write_node(tree_root(ig^D))
!end do
!}
contains
!=============================================================================
! internal procedures
!=============================================================================
recursive subroutine write_node(tree)
implicit none

type(tree_node_ptr) :: tree

integer :: ic^D
!-----------------------------------------------------------------------------
count_forest = count_forest+1
if(typeparIO /= -1) then
  call MPI_FILE_WRITE(file_handle,tree%node%leaf,1,MPI_LOGICAL, &
                    status,ierrmpi)
else
  write(file_handle) tree%node%leaf
end if

if (.not.tree%node%leaf) then
   {do ic^DB=1,2\}
      call write_node(tree%node%child(ic^D))
   {end do\}
end if

end subroutine write_node
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine write_forest
!=============================================================================
subroutine read_forest(file_handle, N_flat_forest)
use mod_forest
include 'amrvacdef.f'

integer, intent(in) :: file_handle, N_flat_forest

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(MPI_STATUS_SIZE) :: status
!integer :: ig^D, level, size_logical, Morton_no, igrid, ipe
integer :: ig^D, level, Morton_no, igrid, ipe, isfc, array_offset
integer(kind=MPI_ADDRESS_KIND) :: size_logical, lb
logical, allocatable :: flatten_forest(:)

integer, external :: getnode
!-----------------------------------------------------------------------------
!call MPI_TYPE_EXTENT(MPI_LOGICAL,size_logical,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,size_logical,ierrmpi)
offset=int(size_block_io{#IFDEF STAGGERED + size_block_stg_io},kind=MPI_OFFSET_KIND)*int(nleafs,kind=MPI_OFFSET_KIND)

allocate(flatten_forest(1:N_flat_forest))
array_offset = 0
! read the flattened forest all at once
if (typeparIO==1) then
   call MPI_FILE_READ_AT_ALL(file_handle,offset+1,flatten_forest,N_flat_forest,MPI_LOGICAL, &
                              status,ierrmpi)
else
 if (mype==0) then
   call MPI_FILE_READ_AT(file_handle,offset+1,flatten_forest,N_flat_forest,MPI_LOGICAL, &
                             status,ierrmpi)
 end if
 if (npe>1)  call MPI_BCAST(flatten_forest,N_flat_forest,MPI_LOGICAL,0,icomm,ierrmpi)
end if

Morton_no=0
ipe=0
level=1
nleafs_level(1:nlevelshi) = 0
nparents = 0

call get_Morton_range
call level1_Morton_order
do isfc=1,nglev1
   ig^D=sfc_iglevel1(^D,isfc)\
   allocate(tree_root(ig^D)%node)
   nullify(tree_root(ig^D)%node%parent%node)
   call read_node(tree_root(ig^D),ig^D,level)
end do


!{#IFNDEF SFC_ON_BASE
!{do ig^DB=1,ng^DB(1)\}
!   allocate(tree_root(ig^D)%node)
!   nullify(tree_root(ig^D)%node%parent%node)
!   call read_node(tree_root(ig^D),ig^D,level)
!{end do\}
!}{#IFDEF SFC_ON_BASE
!call level1_Morton_order
!do isfc=1,nglev1
!   ig^D=sfc_iglevel1(^D,isfc)\
!   allocate(tree_root(ig^D)%node)
!   nullify(tree_root(ig^D)%node%parent%node)
!   call read_node(tree_root(ig^D),ig^D,level)
!end do
!}

call get_level_range

! Rebuild tree connectivity
call getigrids
call build_connectivity

! This call is here to ensure the sfc array is initialized
call amr_Morton_order()

deallocate(flatten_forest)

contains
!=============================================================================
! internal procedures
!=============================================================================
recursive subroutine read_node(tree,ig^D,level)
implicit none

type(tree_node_ptr) :: tree
integer, intent(in) :: ig^D, level

logical :: leaf
integer :: ic^D, child_ig^D, child_level
!-----------------------------------------------------------------------------

array_offset=array_offset+1
leaf = flatten_forest(array_offset)

tree%node%leaf=leaf
tree%node%ig^D=ig^D;
tree%node%level=level
tree%node%active=.true. .and. leaf

{do ic^DB=1,2\}
   nullify(tree%node%child(ic^D)%node)
{end do\}
nullify({tree%node%neighbor(1,^D)%node},{tree%node%neighbor(2,^D)%node})
nullify(tree%node%next%node,tree%node%prev%node)

call asign_tree_neighbor(tree)

if (leaf) then
   call add_to_linked_list(level,tree)
   nleafs_level(level) = nleafs_level(level) + 1

   Morton_no=Morton_no+1
   if (Morton_no>Morton_stop(ipe)) ipe=ipe+1
   igrid=getnode(ipe)
   tree%node%igrid=igrid
   tree%node%ipe=ipe
   igrid_to_node(igrid,ipe)%node => tree%node
   if (ipe==mype) sfc_to_igrid(Morton_no)=igrid
else
   nparents = nparents + 1
   tree%node%igrid=0
   tree%node%ipe=-1
   child_level=level+1
   {do ic^DB=1,2\}
      child_ig^D=2*(ig^D-1)+ic^D;
      allocate(tree%node%child(ic^D)%node)
      tree%node%child(ic^D)%node%parent%node => tree%node
      call read_node(tree%node%child(ic^D),child_ig^D,child_level)
   {end do\}
end if

end subroutine read_node
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine read_forest
!=============================================================================
{#IFDEF HDF5
subroutine read_forest_hdf5(file_id)
use mod_forest
use HDF5
include 'amrvacdef.f'

integer(hid_t), intent(in) :: file_id
integer(hid_t) :: dset_handle,dspace_handle

integer :: error,alloc_status
integer, dimension(:), allocatable :: all_levels
integer(hsize_t), dimension(1) :: dims,maxdims
integer :: offset=1
integer :: ig^D, level, Morton_no, igrid, ipe

integer, external :: getnode
!-----------------------------------------------------------------------------
call H5Dopen_f(file_id,"Levels",dset_handle,error)
call H5Dget_space_f(dset_handle,dspace_handle,error)
call H5Sget_simple_extent_dims_f(dspace_handle, dims, maxdims, error)
allocate(all_levels(1:dims(1)), stat=alloc_status)
call H5Dread_f(dset_handle,H5T_NATIVE_INTEGER,all_levels,dims,error)
call H5Sclose_f(dspace_handle,error)
call H5Dclose_f(dset_handle,error)

Morton_no=0
ipe=0
level=1
nleafs_level(1:nlevelshi) = 0

call get_Morton_range
{do ig^DB=1,ng^DB(1)\}
   allocate(tree_root(ig^D)%node)
   nullify(tree_root(ig^D)%node%parent%node)
   call read_node(tree_root(ig^D),ig^D,level)
{end do\}

deallocate(all_levels,stat=alloc_status)

call get_level_range

! Rebuild tree connectivity
call getigrids
call build_connectivity

contains
!=============================================================================
! internal procedures
!=============================================================================
recursive subroutine read_node(tree,ig^D,level)
implicit none

type(tree_node_ptr) :: tree
integer, intent(in) :: ig^D, level
integer :: level_in

logical :: leaf
integer :: ic^D, child_ig^D, child_level
!-----------------------------------------------------------------------------
level_in=all_levels(offset)
if(level_in==level) then
  leaf=.true.
  offset=offset+1
else
  leaf=.false.
end if

tree%node%leaf=leaf
tree%node%ig^D=ig^D;
tree%node%level=level
tree%node%active=.true. .and. leaf

{do ic^DB=1,2\}
   nullify(tree%node%child(ic^D)%node)
{end do\}
nullify({tree%node%neighbor(1,^D)%node},{tree%node%neighbor(2,^D)%node})
nullify(tree%node%next%node,tree%node%prev%node)

call asign_tree_neighbor(tree)

if (leaf) then
   call add_to_linked_list(level,tree)
   nleafs_level(level) = nleafs_level(level) + 1

   Morton_no=Morton_no+1
   if (Morton_no>Morton_stop(ipe)) ipe=ipe+1
   igrid=getnode(ipe)
   tree%node%igrid=igrid
   tree%node%ipe=ipe
   igrid_to_node(igrid,ipe)%node => tree%node
   if (ipe==mype) sfc_to_igrid(Morton_no)=igrid
else
   tree%node%igrid=0
   tree%node%ipe=-1
   child_level=level+1
   {do ic^DB=1,2\}
      child_ig^D=2*(ig^D-1)+ic^D;
      allocate(tree%node%child(ic^D)%node)
      tree%node%child(ic^D)%node%parent%node => tree%node
      call read_node(tree%node%child(ic^D),child_ig^D,child_level)
   {end do\}
end if

end subroutine read_node
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine read_forest_hdf5
!
\}
