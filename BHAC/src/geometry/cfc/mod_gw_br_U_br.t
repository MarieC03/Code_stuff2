! fixme: remember the BCs for psi is reducing to 1, but BCs for U and U5 reducing to 0
module mod_gw_br_U_br
  use mod_cfc_parameters
  implicit none
  private

  ! public methods
  public :: gw_br_solve_U_br

  contains

  subroutine gw_br_U_br_solver_init(sigma)
    use mod_multigrid_coupling
    use mod_forest
    use mod_metric
    include 'amrvacdef.f'

    double precision, intent(in)   :: sigma(ixG^T,1:igridstail)
    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: rhs(ixG^T,1:igridstail)
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;
    mg%operator_type = mg_laplacian
    call mg_set_methods(mg)

    select case (coordinate)
    case (cartesian)
       do nb=1, 2*ndim
          select case (typeB(U_br_, nb)) 
          case ('cont','noinflow')
             ! outer boundary, use Robin B. C.
             !mg%bc(nb, mg_iphi)%bc_type = mg_bc_dirichlet
             !mg%bc(nb, mg_iphi)%bc_value = 0.0d0
             ! psi, alp BCs
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
             mg%bc(nb, mg_iphi)%bc_value = 1.0d0
          case ('symm')
             ! inner boundary, use neumann !fixme: is it okay?
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case default
             call mpistop("Error: The boundary conditions of U_br are not set correctly.")
          end select
       end do
    case (cylindrical)
       call mpistop('GW_BR does not support cylindrical yet')
    case (spherical)
       call mpistop('GW_BR does not support spherical yet')
    end select

    ! copy the data into MG solver, and set up the source terms
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! layer of ghost cells on grid leaves are not included here

       ! initial guess of (U_br-1)
       !mg%boxes(id)%cc({:,}, mg_iphi) = (pw(igrid)%w(ix^S, U_br_) - 1.0d0)

       ! source term 1 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp1) = 0.0d0

       ! source term 2 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp2) = 0.0d0

       !rhs(ixM^T,iigrid) = -4.0d0 * dpi * pw(igrid)%w(ixM^T,D_) * pw(igrid)%w(ixM^T,psi_metric_)**6
       rhs(ixM^T,iigrid) = -4.0d0 * dpi * sigma(ixM^T,iigrid) 
       mg%boxes(id)%cc({1:nc}, mg_irhs) = rhs(ixM^T,iigrid)
    end do
  end subroutine gw_br_U_br_solver_init

  subroutine gw_br_solve_U_br(sigma)
    use mod_multigrid_coupling
    use mod_forest
    include 'amrvacdef.f'

    double precision, intent(in)   :: sigma(ixG^T,1:igridstail)
    integer                        :: id, iigrid, igrid
    double precision               :: res = huge(1.0d0)
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    integer                        :: ix^L

    call gw_br_U_br_solver_init(sigma(ixG^T,1:igridstail))

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, .True., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       if (res <= gw_br_tol(1)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving U_br:', it, mg_it, res
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge U_br."
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

       pw(igrid)%w(ix^S, U_br_) = mg%boxes(id)%cc({:,}, mg_iphi)
       !pw(igrid)%w(ix^S, U_br_) = mg%boxes(id)%cc({:,}, mg_iphi) + 1.0d0
    end do

  end subroutine gw_br_solve_U_br
  
  ! fixme: add a checker to check GW emission low and skip it.

end module mod_gw_br_U_br
