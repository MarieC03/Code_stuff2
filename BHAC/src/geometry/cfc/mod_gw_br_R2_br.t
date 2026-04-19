module mod_gw_br_R2_br
  use mod_cfc_parameters
  use mod_imhd_intermediate

  implicit none
  private

  ! public methods
  public :: gw_br_solve_R2_br

  contains

  subroutine gw_br_R2_br_solver_init(R_br, sigma, gamma, w_i)
    use mod_multigrid_coupling
    use mod_forest
    use mod_metric
    include 'amrvacdef.f'

    double precision, intent(in)   :: R_br(ixG^T,1:igridstail)
    double precision, intent(in)   :: sigma(ixG^T,1:igridstail)
    double precision, intent(in)   :: gamma(ixG^T,1:3,1:3,1:igridstail), w_i(ixG^T,1:3,1:igridstail)
    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: rhs(ixG^T)
    integer                        :: idir, i ,j , ix^D
    double precision               :: dU(ixG^T,1:3,1:igridstail), d_U(ixG^T,1:3,1:igridstail), &
                                      Q3_ij_x_dU(ixG^T,1:igridstail), Q3_ij_local(ixG^T,1:3,1:3,1:igridstail), &
                                      Q3_ij_wiwj(ixG^T,1:igridstail)
    double precision               :: U(ixG^T,1:igridstail)
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;
    mg%operator_type = mg_laplacian
    call mg_set_methods(mg)

    !select case (coordinate)
    !case (cartesian)
    !   do nb=1, 2*ndim
    !      select case (typeB(U5_br_, nb)) 
    !      case ('cont','noinflow')
    !         ! outer boundary, use Robin B. C.
    !         !mg%bc(nb, mg_iphi)%bc_type = mg_bc_dirichlet
    !         !mg%bc(nb, mg_iphi)%bc_value = 0.0d0
    !         mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
    !         mg%bc(nb, mg_iphi)%bc_value = 1.0d0
    !      case ('symm')
    !         ! inner boundary, use neumann !fixme: is it okay?
    !         mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
    !         mg%bc(nb, mg_iphi)%bc_value = 0.0d0
    !      case default
    !         call mpistop("Error: The boundary conditions of U5_br are not set correctly.")
    !      end select
    !   end do
    !case (cylindrical)
    !   call mpistop('GW_BR does not support cylindrical yet')
    !case (spherical)
    !   call mpistop('GW_BR does not support spherical yet')
    !end select
    do nb=1, 2*ndim
       if (nb .ne. 5) then
          mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
          mg%bc(nb, mg_iphi)%bc_value = 1.0d0
       else
          mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(nb, mg_iphi)%bc_value = 0.0d0
       endif
    enddo

    ! copy the data into MG solver, and set up the source terms
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! layer of ghost cells on grid leaves are not included here

       call set_tmpGlobals(igrid)

       ! initial guess of (U5_br-1)
       !mg%boxes(id)%cc({:,}, mg_iphi) = (pw(igrid)%w(ix^S, U5_br_) - 1.0d0)

       ! source term 1 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp1) = 0.0d0

       ! source term 2 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp2) = 0.0d0


       d_U(ixG^T,:,iigrid) = 0.0d0
       U(ixG^T,iigrid) = pw(igrid)%w(ixG^T, U_br_)
   
       do idir = 1, ndim
         call partial_d( U(ixG^T,iigrid),ixG^LL,ixM^LL,px(igrid)%x,idir,d_U(ixG^T,idir,iigrid) )
       enddo

       if (use_index_contract) then
         do idir = 1, ndir
            dU(ixM^T,idir,iigrid) = d_U(ixM^T,idir,iigrid) / gamma(ixM^T,idir,idir,iigrid)
         enddo
       else
         dU(ixG^T,:,iigrid) = d_U(ixG^T,:,iigrid)
       endif

       Q3_ij_local(ixM^T,1,1,iigrid) = I3ij_new_global(1)
       Q3_ij_local(ixM^T,1,2,iigrid) = I3ij_new_global(2)
       Q3_ij_local(ixM^T,1,3,iigrid) = I3ij_new_global(3)
       Q3_ij_local(ixM^T,2,2,iigrid) = I3ij_new_global(4)
       Q3_ij_local(ixM^T,2,3,iigrid) = I3ij_new_global(5)
       Q3_ij_local(ixM^T,3,3,iigrid) = I3ij_new_global(6)

       Q3_ij_local(ixM^T,2,1,iigrid) = I3ij_new_global(2)
       Q3_ij_local(ixM^T,3,1,iigrid) = I3ij_new_global(3)
       Q3_ij_local(ixM^T,3,2,iigrid) = I3ij_new_global(5)

       Q3_ij_x_dU(ixG^T,iigrid) = 0.0d0
       do i = 1, 3
          do j = 1, 3
             Q3_ij_x_dU(ixM^T,iigrid) = Q3_ij_x_dU(ixM^T,iigrid) + &
                                 Q3_ij_local(ixM^T,i,j,iigrid) * px(igrid)%x(ixM^T, i) * dU(ixM^T,j,iigrid)
          enddo
       enddo

       Q3_ij_wiwj(ixM^T,iigrid) = 0.0d0
       do i = 1, 3
          do j = 1, 3
             Q3_ij_wiwj(ixM^T,iigrid) = Q3_ij_wiwj(ixM^T,iigrid) + &
                                 3.0d0 * w_i(ixM^T, i, iigrid) * w_i(ixM^T, j, iigrid) * Q3_ij_local(ixM^T,i,j,iigrid)
          enddo
       enddo
   
       rhs(ixM^T) = -4.0d0 * dpi * sigma(ixM^T,iigrid) *&
       !rhs(ixM^T) = -4.0d0 * dpi * pw(igrid)%w(ixM^T,D_) * pw(igrid)%w(ixM^T,psi_metric_)**6 *&
                     ( Q3_ij_x_dU(ixM^T,iigrid) - R_br(ixM^T,iigrid) + Q3_ij_wiwj(ixM^T,iigrid) )
       mg%boxes(id)%cc({1:nc}, mg_irhs) = rhs(ixM^T)
    end do
  end subroutine gw_br_R2_br_solver_init

  subroutine gw_br_solve_R2_br(R2_br, R_br, sigma, gamma, w_i)
    use mod_multigrid_coupling
    use mod_forest
    include 'amrvacdef.f'

    double precision, intent(inout):: R2_br(ixG^T,1:igridstail)
    double precision, intent(in)   :: R_br(ixG^T,1:igridstail)
    double precision, intent(in)   :: sigma(ixG^T,1:igridstail), w_i(ixG^T,1:3,1:igridstail)
    double precision, intent(in)   :: gamma(ixG^T,1:3,1:3,1:igridstail)
    integer                        :: id, iigrid, igrid
    double precision               :: res = huge(1.0d0)
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    integer                        :: ix^L

    call gw_br_R2_br_solver_init(R_br(ixG^T,1:igridstail), &
                                 sigma(ixG^T,1:igridstail), gamma(ixG^T,1:3,1:3,1:igridstail), w_i(ixG^T,1:3,1:igridstail))

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, .True., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       if (res <= gw_br_tol(1)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving R2_br:', it, mg_it, res
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge R2_br."
          write(*,*) "! it=", it,"N = ", mg_it, " Res = ", res
          write(*,*) "! Maybe the gw_br_tol is too small or cfc_it_max is too small"
       end if
    end if

    ! copy the data from MG solver
    ix^L=ixM^LL^LADD1;
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       R2_br(ix^S,iigrid) = mg%boxes(id)%cc({:,}, mg_iphi)
       !pw(igrid)%w(ix^S, U5_br_) = mg%boxes(id)%cc({:,}, mg_iphi)
       !pw(igrid)%w(ix^S, U5_br_) = mg%boxes(id)%cc({:,}, mg_iphi) + 1.0d0
    end do

  end subroutine gw_br_solve_R2_br

end module mod_gw_br_R2_br
