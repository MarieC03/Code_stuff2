! fixme: remember the BCs for psi is reducing to 1, but BCs for U and U5 reducing to 0
module mod_gw_br_R_br
  use mod_cfc_parameters
  use mod_imhd_intermediate

  implicit none
  private

  ! public methods
  public :: gw_br_solve_R_br

  contains

  subroutine gw_br_R_br_solver_init(sigma, gamma)
    use mod_multigrid_coupling
    use mod_forest
    use mod_metric
    include 'amrvacdef.f'

    double precision, intent(in)   :: sigma(ixG^T,1:igridstail)
    double precision, intent(in)   :: gamma(ixG^T,1:3,1:3,1:igridstail)
    integer                        :: id, iigrid, igrid, i, j, idir
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: rhs(ixG^T), rho_st(ixG^T,1:igridstail), &
                                      d_sigma(ixG^T,1:ndir,1:igridstail), &
                                      dsigma(ixG^T,1:ndir,1:igridstail), Q3_ij_x_dsigma(ixG^T,1:igridstail)
    double precision               :: Q3_ij_local(ixG^T,1:3,1:3,1:igridstail)
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;
    mg%operator_type = mg_laplacian
    call mg_set_methods(mg)

    !select case (coordinate)
    !case (cartesian)
    !   do nb=1, 2*ndim
    !      select case (typeB(R_br_, nb))
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
    !         call mpistop("Error: The boundary conditions of R_br are not set correctly.")
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

       call set_tmpGlobals(igrid)

       ! layer of ghost cells on grid leaves are not included here

       ! initial guess of (R_br-1)
       !mg%boxes(id)%cc({:,}, mg_iphi) = (pw(igrid)%w(ix^S, R_br_) - 1.0d0)

       Q3_ij_local(ixM^T,1,1,iigrid) = I3ij_new_global(1)
       Q3_ij_local(ixM^T,1,2,iigrid) = I3ij_new_global(2)
       Q3_ij_local(ixM^T,1,3,iigrid) = I3ij_new_global(3)
       Q3_ij_local(ixM^T,2,2,iigrid) = I3ij_new_global(4)
       Q3_ij_local(ixM^T,2,3,iigrid) = I3ij_new_global(5)
       Q3_ij_local(ixM^T,3,3,iigrid) = I3ij_new_global(6)
 
       Q3_ij_local(ixM^T,2,1,iigrid) = I3ij_new_global(2)
       Q3_ij_local(ixM^T,3,1,iigrid) = I3ij_new_global(3)
       Q3_ij_local(ixM^T,3,2,iigrid) = I3ij_new_global(5)

       rho_st(ixG^T,iigrid) = pw(igrid)%w(ixG^T,D_) * pw(igrid)%w(ixG^T,psi_metric_)**6
     
       d_sigma(ixG^T,:,iigrid) = 0.0d0 
       do idir = 1, ndim
         ! actually it is d_rho_st not d_sigma
         !call partial_d( rho_st(ixG^T,iigrid),ixG^LL,ixM^LL,px(igrid)%x(ixG^T,1:ndim),idir,d_sigma(ixG^T,idir,iigrid) )
         call partial_d( sigma(ixG^T,iigrid),ixG^LL,ixM^LL,px(igrid)%x(ixG^T,1:ndim),idir,d_sigma(ixG^T,idir,iigrid) )
       enddo

       if (use_index_contract) then
         do idir = 1, ndir
            dsigma(ixM^T, idir,iigrid) = d_sigma(ixM^T,idir,iigrid) / gamma(ixM^T,idir,idir,iigrid)
         enddo
       else
         dsigma(ixM^T,:,iigrid) = d_sigma(ixM^T,:,iigrid)
       endif
       
       ! it's Q_ij x^i d^j sigma
       Q3_ij_x_dsigma(ixG^T,iigrid) = 0.0d0
       do i = 1, 3
         do j = 1, 3
            Q3_ij_x_dsigma(ixM^T,iigrid) = Q3_ij_x_dsigma(ixM^T,iigrid) +&
                                    Q3_ij_local(ixM^T,i,j,iigrid) * px(igrid)%x(ixM^T, i) * dsigma(ixM^T, j, iigrid)
         enddo
       enddo

       ! source term 1 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp1) = 0.0d0

       ! source term 2 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp2) = 0.0d0

       rhs(ixM^T) = -4.0d0 * dpi * Q3_ij_x_dsigma(ixM^T,iigrid)  
       mg%boxes(id)%cc({1:nc}, mg_irhs) = rhs(ixM^T)
    end do
  end subroutine gw_br_R_br_solver_init

  subroutine gw_br_solve_R_br(R_br, sigma, gamma)
    use mod_multigrid_coupling
    use mod_forest
    include 'amrvacdef.f'

    double precision, intent(inout):: R_br(ixG^T,1:igridstail)
    double precision, intent(in)   :: sigma(ixG^T,1:igridstail)
    double precision, intent(in)   :: gamma(ixG^T,1:3,1:3,1:igridstail)
    integer                        :: id, iigrid, igrid
    double precision               :: res = huge(1.0d0)
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    integer                        :: ix^L

    call gw_br_R_br_solver_init(sigma(ixG^T,1:igridstail), gamma(ixG^T,1:3,1:3,1:igridstail))

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, .True., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       if (res <= gw_br_tol(3)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving R_br:', it, mg_it, res
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge R_br."
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

       R_br(ix^S, iigrid) = mg%boxes(id)%cc({:,}, mg_iphi)
       !pw(igrid)%w(ix^S, R_br_) = mg%boxes(id)%cc({:,}, mg_iphi)
       !pw(igrid)%w(ix^S, R_br_) = mg%boxes(id)%cc({:,}, mg_iphi) + 1.0d0
    end do

  end subroutine gw_br_solve_R_br
  
  ! fixme: add a checker to check GW emission low and skip it.

end module mod_gw_br_R_br
