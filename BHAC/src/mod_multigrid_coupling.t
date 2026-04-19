!> Module to couple the octree-mg library to GMUNU. This file uses the VACPP
!> preprocessor, but its use is kept to a minimum.
module mod_multigrid_coupling
  {^IFONED
  use m_octree_mg_1d
  }
  {^IFTWOD
  use m_octree_mg_2d
  }
  {^IFTHREED
  use m_octree_mg_3d
  }


  implicit none
  public

  !> Data structure containing the multigrid tree.
  type(mg_t) :: mg

  !> If defined, this routine is called after a new multigrid tree is
  !> constructed.
  procedure(after_new_tree), pointer :: mg_after_new_tree => null()

  interface
     subroutine after_new_tree()
     end subroutine after_new_tree
  end interface

contains

  !> Setup multigrid for usage
  subroutine mg_setup_multigrid()
    use mod_metric
    use mod_cfc_parameters
    include 'amrvacdef.f'

    if (ndim /= mg_ndim) &
         error stop "Multigrid module was compiled for different ndim"
{#IFDEF DY_SP
    select case (coordinate)
    case (Cartesian)
       mg%geometry_type = mg_cartesian
       continue
    case (cylindrical)
       !if (ndim == 3) error stop "Multigrid does not support cylindrical 3D"
       mg%geometry_type = mg_cylindrical
    case default
       mg%geometry_type = mg_spherical
    end select
}
{#IFNDEF DY_SP
    select case (typeaxial)
    case ("slab")
       if (ndim == 1) error stop "Multigrid only support 2D, 3D"
    case ("cylindrical")
       if (ndim == 3) error stop "Multigrid does not support cylindrical 3D"
       mg%geometry_type = mg_cylindrical
    case default
       error stop "Multigrid does not support your geometry"
    end select
}

    if (any([ block_nx^D ] /= block_nx1)) &
         error stop "Multigrid requires all block_nx to be equal"

    call mg_comm_init(mg)
    !call mg_set_methods(mg)
    call mg_tree_from_bhac(mg)
  end subroutine mg_setup_multigrid

  !> Set multigrid boundary conditions for the solution according to variable iw
  subroutine mg_copy_boundary_conditions(mg, iw)
    include 'amrvacdef.f'
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iw
    integer                   :: n

    do n = 1, mg_num_neighbors
       select case (typeB(iw, n))
       case ('symm')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('cont')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_continuous
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp ! Not needed
       case ('periodic')
          ! Nothing to do here
       case default
          print *, "Not a standard: ", trim(typeB(iw, n))
          error stop "You have to set a user-defined boundary method"
       end select
    end do
  end subroutine mg_copy_boundary_conditions

  !> If the grid has changed, rebuild the full multigrid tree
  subroutine mg_update_refinement(n_coarsen, n_refine)
    include 'amrvacdef.f'
    integer, intent(in) :: n_coarsen
    integer, intent(in) :: n_refine

    ! Don't build multigrid tree while doing initial refinement
    if (.not. time_advance) return

    if (.not. mg%is_allocated) then
       call mg_tree_from_bhac(mg)
    else if (n_coarsen + n_refine > 0) then
       call mg_deallocate_storage(mg)
       call mg_tree_from_bhac(mg)
    end if
  end subroutine mg_update_refinement

  !> Copy a variable to the multigrid tree, including a layer of ghost cells
  subroutine mg_copy_to_tree(iw_from, iw_to, restrict, restrict_gc, factor, state_from)
    use mod_forest
    include 'amrvacdef.f'
    integer, intent(in)            :: iw_from     !< Variable to use as right-hand side
    integer, intent(in)            :: iw_to       !< Copy to this variable
    logical, intent(in), optional  :: restrict    !< Restrict variable on multigrid tree
    logical, intent(in), optional  :: restrict_gc !< Fill ghost cells after restrict
    real(dp), intent(in), optional :: factor      !< out = factor * in
    !> If present, use this as the input state
    type(state), intent(in), optional, target :: state_from(ngridshi)
    integer                        :: iigrid, igrid, id
    integer                        :: nc, lvl, ilo^D, ihi^D
    type(tree_node), pointer       :: pnode
    real(dp)                       :: fac
    logical                        :: do_restrict, do_gc

    if (.not. mg%is_allocated) &
         error stop "mg_copy_to_tree: tree not allocated yet"

    fac = 1.0_dp; if (present(factor)) fac = factor
    do_restrict = .false.; if (present(restrict)) do_restrict = restrict
    do_gc = .false.; if (present(restrict_gc)) do_gc = restrict_gc

    ilo^D=ixMlo^D-1;
    ihi^D=ixMhi^D+1;

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Include one layer of ghost cells on grid leaves
       if (present(state_from)) then
          mg%boxes(id)%cc({0:nc+1}, iw_to) = &
               fac * state_from(igrid)%w%w({ilo^D:ihi^D}, iw_from)
       else
          mg%boxes(id)%cc({0:nc+1}, iw_to) = fac * pw(igrid)%w({ilo^D:ihi^D}, iw_from)
       end if
    end do

    if (do_restrict) then
       call mg_restrict(mg, iw_to)
       if (do_gc) call mg_fill_ghost_cells(mg, iw_to)
    end if

  end subroutine mg_copy_to_tree

  !> Copy a variable from the multigrid tree
  subroutine mg_copy_from_tree(iw_from, iw_to, state_to)
    use mod_forest
    include 'amrvacdef.f'
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    !> If present, use this as the output state
    type(state), intent(inout), optional, target :: state_to(ngridshi)
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "mg_copy_from_tree: tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       if (present(state_to)) then
          state_to(igrid)%w%w(ixM^T, iw_to) = mg%boxes(id)%cc({1:nc}, iw_from)
       else
          pw(igrid)%w(ixM^T, iw_to) = mg%boxes(id)%cc({1:nc}, iw_from)
       end if
    end do
  end subroutine mg_copy_from_tree

  !> Copy from multigrid tree with one layer of ghost cells. Corner ghost cells
  !> are not used/set.
  subroutine mg_copy_from_tree_gc(iw_from, iw_to, state_to)
    use mod_forest
    include 'amrvacdef.f'
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    !> If present, use this as the output state
    type(state), intent(inout), optional, target :: state_to(ngridshi)
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl, ilo^D, ihi^D
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "mg_copy_from_tree_gc: tree not allocated yet"

    ilo^D=ixMlo^D-1;
    ihi^D=ixMhi^D+1;

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       if (present(state_to)) then
          state_to(igrid)%w%w({ilo^D:ihi^D}, iw_to) = mg%boxes(id)%cc({0:nc+1}, iw_from)
       else
          pw(igrid)%w({ilo^D:ihi^D}, iw_to) = mg%boxes(id)%cc({0:nc+1}, iw_from)
       end if
    end do
  end subroutine mg_copy_from_tree_gc

  !> Generate a multigrid tree that includes the bhac tree, but also contains
  !> coarser grid levels. A number of checks has already been performed in
  !> mg_setup_multigrid, so we don't repeat these checks here.
  subroutine mg_tree_from_bhac(mg)
    use mod_forest
    include 'amrvacdef.f'
    type(mg_t), intent(inout)        :: mg
    integer                          :: i, n, id, ix(ndim)
    integer                          :: n_boxes_total, i_c, c_id, c_ix(ndim)
    integer                          :: min_lvl, lvl
    integer                          :: nb, nb_ix, nb_dim
    integer                          :: n_finer
    integer                          :: nxlone(ndim)

    type(tree_node), pointer         :: pnode, pnode_ch
    type(tree_node_ptr), allocatable :: id_to_node(:)
    real(dp)                         :: dr_coarse

    ! Estimate number of finer blocks
    n_finer = nparents+nleafs

    nxlone = [ ng^D(1) ] * block_nx1

    call mg_build_rectangle(mg, nxlone, block_nx1, &
         dx(:,1), [ xprobmin^D ], periodB, n_finer)

    mg%highest_lvl = levmax
    n_boxes_total = mg%n_boxes + n_finer

    ! To link the two trees
    allocate(id_to_node(n_boxes_total))

    ! Link base level
    do i = 1, size(mg%lvls(1)%ids)
       id = mg%lvls(1)%ids(i)
       ix = mg%boxes(id)%ix

       pnode               => tree_root({ix(^D)})%node
       pnode%id            =  id
       id_to_node(id)%node => pnode
       mg%boxes(id)%rank   =  pnode%ipe
    end do

    ! Add refinement
    do lvl = 1, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          pnode => id_to_node(id)%node

          if (.not. pnode%leaf) then
             call mg_add_children(mg, id)

             do i_c = 1, mg_num_children
                c_id = mg%boxes(id)%children(i_c)
                c_ix = mg_child_dix(:, i_c) + 1
                pnode_ch => pnode%child({c_ix(^D)})%node
                id_to_node(c_id)%node => pnode_ch
                pnode_ch%id = c_id
                mg%boxes(c_id)%rank = pnode_ch%ipe
             end do
          end if
       end do

       call mg_set_leaves_parents(mg%boxes, mg%lvls(lvl))

       if (lvl < mg%highest_lvl) then
          call mg_set_next_level_ids(mg, lvl)
          call mg_set_neighbors_lvl(mg, lvl+1)
       end if
    end do

    ! Store boxes with refinement boundaries (from the coarse side)
    do lvl = 1, mg%highest_lvl
       call mg_set_refinement_boundaries(mg%boxes, mg%lvls(lvl))
    end do

    ! Assign boxes to MPI processes
    call mg_load_balance_parents(mg)

    ! Allocate storage for boxes owned by this process
    call mg_allocate_storage(mg)

    if (associated(mg_after_new_tree)) then
       call mg_after_new_tree()
    end if

  end subroutine mg_tree_from_bhac


  subroutine clean_divb_multigrid()
    use mod_forest
    include 'amrvacdef.f'
    integer                     :: iigrid, igrid, id
    integer                     :: n, nc, lvl, ix^L, idim
    type(tree_node), pointer    :: pnode
    double precision            :: tmp(ixG^T), grad(ixG^T, ndim)
    double precision, parameter :: residual_reduction = 1d-10
    integer, parameter          :: max_its            = 50
    double precision            :: residual_it(max_its)
    integer                     :: i, j, i0, j0
    real(dp)                    :: phi_grad, max_divb
    {^IFTHREED
    integer                     :: k, k0
    }

    mg%operator_type = mg_laplacian

    ! TODO (?) Set boundary conditions
    mg%bc(:, mg_iphi)%bc_type = mg_bc_dirichlet
    mg%bc(:, mg_iphi)%bc_value = 0.0_dp

    ix^L=ixM^LL^LADD1;
    max_divb = 0.0d0

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       {#IFDEF STAGGERED
       call div_staggered(ixG^LL,ixM^LL,ps(igrid),tmp(ixG^T))
       }{#IFNDEF STAGGERED
       ! Reduce output array size, +1 was added for eventual pointdata output
       call get_divb(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),tmp)
       }
       mg%boxes(id)%cc({1:nc}, mg_irhs) = tmp(ixM^T)
       max_divb = max(max_divb, maxval(abs(tmp(ixM^T))))
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, max_divb, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, icomm, ierrmpi)

    ! Solve laplacian(phi) = divB
    if (mype == 0) print *, "Performing multigrid divB cleaning"
    if (mype == 0) print *, "iteration vs residual"

    do n = 1, max_its
       call mg_fas_fmg(mg, n>1, max_res=residual_it(n))

       ! V-cycles are cheaper but converge less quickly
       ! call mg_fas_vcycle(mg, max_res=residual_it(n))

       if (mype == 0) write(*, "(I4,E11.3)") n, residual_it(n)
       if (residual_it(n) < residual_reduction * max_divb) exit
    end do

    if (mype == 0 .and. n > max_its) then
       print *, "divb_multigrid warning: not fully converged"
       print *, "current amplitude of divb: ", residual_it(max_its)
       print *, "multigrid smallest grid: ", &
            mg%domain_size_lvl(:, mg%lowest_lvl)
       print *, "note: smallest grid ideally has <= 8 cells"
       print *, "multigrid dx/dy/dz ratio: ", mg%dr(:, 1)/mg%dr(1, 1)
       print *, "note: dx/dy/dz should be similar"
    end if

    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       nc    =  mg%box_size

       ! Geometry subroutines expect this to be set
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       ! Compute the gradient of phi
       tmp(ix^S) = mg%boxes(id)%cc({:,}, mg_iphi)

{#IFDEF STAGGERED
       {^IFTWOD
       ! first dimension
       do j = 1, nc
          do i = 0, nc
             phi_grad = (mg%boxes(id)%cc(i+1, j, mg_iphi) - &
                  mg%boxes(id)%cc(i, j, mg_iphi)) / dxlevel(1)
             i0 = i + dixB
             j0 = j + dixB
             pws(igrid)%w(i0, j0, bs1_) = &
                  pws(igrid)%w(i0, j0, bs1_) - phi_grad
          end do
       end do

       ! second dimension
       do j = 0, nc
          do i = 1, nc
             phi_grad = (mg%boxes(id)%cc(i, j+1, mg_iphi) - &
                  mg%boxes(id)%cc(i, j, mg_iphi)) / dxlevel(2)
             i0 = i + dixB
             j0 = j + dixB
             pws(igrid)%w(i0, j0, bs2_) = &
                  pws(igrid)%w(i0, j0, bs2_) - phi_grad
          end do
       end do
       }

       {^IFTHREED
       ! first dimension
       do k = 1, nc
          do j = 1, nc
             do i = 0, nc
                phi_grad = (mg%boxes(id)%cc(i+1, j, k, mg_iphi) - &
                     mg%boxes(id)%cc(i, j, k, mg_iphi)) / dxlevel(1)
                i0 = i + dixB
                j0 = j + dixB
                k0 = k + dixB
                pws(igrid)%w(i0, j0, k0, bs1_) = &
                     pws(igrid)%w(i0, j0, k0, bs1_) - phi_grad
             end do
          end do
       end do

       ! second dimension
       do k = 1, nc
          do j = 0, nc
             do i = 1, nc
                phi_grad = (mg%boxes(id)%cc(i, j+1, k, mg_iphi) - &
                     mg%boxes(id)%cc(i, j, k, mg_iphi)) / dxlevel(2)
                i0 = i + dixB
                j0 = j + dixB
                k0 = k + dixB
                pws(igrid)%w(i0, j0, k0, bs2_) = &
                     pws(igrid)%w(i0, j0, k0, bs2_) - phi_grad
             end do
          end do
       end do

       ! third dimension
       do k = 0, nc
          do j = 1, nc
             do i = 1, nc
                phi_grad = (mg%boxes(id)%cc(i, j, k+1, mg_iphi) - &
                     mg%boxes(id)%cc(i, j, k, mg_iphi)) / dxlevel(3)
                i0 = i + dixB
                j0 = j + dixB
                k0 = k + dixB
                pws(igrid)%w(i0, j0, k0, bs3_) = &
                     pws(igrid)%w(i0, j0, k0, bs3_) - phi_grad
             end do
          end do
       end do

       call faces2centers(ixM^LL, ps(igrid))
    }
}
    end do

  end subroutine clean_divb_multigrid

end module mod_multigrid_coupling

