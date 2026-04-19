! the loop of idim  is very suspicious
module mod_gw_br_U_bri
  use mod_cfc_parameters
  use mod_imhd_intermediate
  implicit none
  private

  ! public methods
  public :: gw_br_solve_U_br1
  public :: gw_br_solve_U_br2
  public :: gw_br_solve_U_br3

  contains

  subroutine gw_br_solve_U_br1(sigma,w_i)
    use mod_multigrid_coupling
    use mod_forest
    use mod_metric
    include 'amrvacdef.f'

    double precision, intent(in)   :: sigma(ixG^T,1:igridstail), w_i(ixG^T,1:3,1:igridstail)
    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: res = huge(1.0d0)
    double precision               :: rhs(ixG^T)
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;
    mg%operator_type = mg_laplacian
    call mg_set_methods(mg)

    do nb = 1, 2*ndim
          select case (typeB(U_br1_, nb))
          case ('symm')
             ! inner boundary, but not reflecting:
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case ('asymm')
             ! reflecting inner boundary:
             !mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             ! fixme: very weird behavior for extremely low res or tightly patched mesh
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case default
             ! outer boundary:
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          end select
          !write(*,*) nb, idir, typeboundary(beta(idir), nb)
    end do

    ! If separated into three normal poisson, so should I treat each of them a scalar? --> then use this
    !do nb=1, 2*ndim
    !   if (nb .ne. 5) then
    !      mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
    !      mg%bc(nb, mg_iphi)%bc_value = 1.0d0
    !   else
    !      mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
    !      mg%bc(nb, mg_iphi)%bc_value = 0.0d0
    !   endif
    !enddo


    ! copy the data into MG solver, and set up the source terms
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)


       ! layer of ghost cells on grid leaves are not included here
       ! initial guess of (U_br^i-1)
       !mg%boxes(id)%cc({:,}, mg_iphi) = (pw(igrid)%w(ix^S, U_br^C_) - 1.0d0)

       ! source term 1 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp1) = 0.0d0

       ! source term 2 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp2) = 0.0d0

       !rhs(ixM^T) = -4.0d0 * dpi * pw(igrid)%w(ixM^T,D_) * pw(igrid)%w(ixM^T,psi_metric_)**6 * w_i(ixM^T,1,iigrid)
       rhs(ixM^T) = -4.0d0 * dpi * sigma(ixM^T,iigrid) * w_i(ixM^T,1,iigrid)
       mg%boxes(id)%cc({1:nc}, mg_irhs) = rhs(ixM^T)
    end do

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, .True., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       if (res <= gw_br_tol(2)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving U_bri_:', it, mg_it, res, 1
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge U_bri_."
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
       pw(igrid)%w(ix^S, U_br1_) = mg%boxes(id)%cc({:,}, mg_iphi)
    end do
  end subroutine gw_br_solve_U_br1

  subroutine gw_br_solve_U_br2(sigma,w_i)
    use mod_multigrid_coupling
    use mod_forest
    use mod_metric
    include 'amrvacdef.f'

    double precision, intent(in)   :: sigma(ixG^T,1:igridstail), w_i(ixG^T,1:3,1:igridstail)
    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: res = huge(1.0d0)
    double precision               :: rhs(ixG^T)
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;
    mg%operator_type = mg_laplacian
    call mg_set_methods(mg)

    do nb = 1, 2*ndim
          select case (typeB(U_br2_, nb))
          case ('symm')
             ! inner boundary, but not reflecting:
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case ('asymm')
             ! reflecting inner boundary:
             !mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             ! fixme: very weird behavior for extremely low res or tightly patched mesh
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case default
             ! outer boundary:
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          end select
          !write(*,*) nb, idir, typeboundary(beta(idir), nb)
    end do

    ! If separated into three normal poisson, so should I treat each of them a scalar? --> then use this
    !do nb=1, 2*ndim
    !   if (nb .ne. 5) then
    !      mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
    !      mg%bc(nb, mg_iphi)%bc_value = 1.0d0
    !   else
    !      mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
    !      mg%bc(nb, mg_iphi)%bc_value = 0.0d0
    !   endif
    !enddo


    ! copy the data into MG solver, and set up the source terms
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)


       ! layer of ghost cells on grid leaves are not included here
       ! initial guess of (U_br^i-1)
       !mg%boxes(id)%cc({:,}, mg_iphi) = (pw(igrid)%w(ix^S, U_br^C_) - 1.0d0)

       ! source term 1 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp1) = 0.0d0

       ! source term 2 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp2) = 0.0d0

       rhs(ixM^T) = -4.0d0 * dpi * sigma(ixM^T,iigrid) * w_i(ixM^T,2,iigrid)
       mg%boxes(id)%cc({1:nc}, mg_irhs) = rhs(ixM^T)
    end do

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, .True., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       if (res <= gw_br_tol(2)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving U_bri_:', it, mg_it, res, 2
       end if
    end do

    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge U_bri_."
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
       pw(igrid)%w(ix^S, U_br2_) = mg%boxes(id)%cc({:,}, mg_iphi)
    end do
  end subroutine gw_br_solve_U_br2

  subroutine gw_br_solve_U_br3(sigma,w_i)
    use mod_multigrid_coupling
    use mod_forest
    use mod_metric
    include 'amrvacdef.f'

    double precision, intent(in)   :: sigma(ixG^T,1:igridstail), w_i(ixG^T,1:3,1:igridstail)
    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: res = huge(1.0d0)
    double precision               :: rhs(ixG^T)
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;
    mg%operator_type = mg_laplacian
    call mg_set_methods(mg)

    do nb = 1, 2*ndim
          select case (typeB(U_br3_, nb))
          case ('symm')
             ! inner boundary, but not reflecting:
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case ('asymm')
             ! reflecting inner boundary:
             !mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             ! fixme: very weird behavior for extremely low res or tightly patched mesh
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case default
             ! outer boundary:
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          end select
          !write(*,*) nb, idir, typeboundary(beta(idir), nb)
    end do

    ! If separated into three normal poisson, so should I treat each of them a scalar? --> then use this
    !do nb=1, 2*ndim
    !   if (nb .ne. 5) then
    !      mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
    !      mg%bc(nb, mg_iphi)%bc_value = 1.0d0
    !   else
    !      mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
    !      mg%bc(nb, mg_iphi)%bc_value = 0.0d0
    !   endif
    !enddo


    ! copy the data into MG solver, and set up the source terms
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)


       ! layer of ghost cells on grid leaves are not included here
       ! initial guess of (U_br^i-1)
       !mg%boxes(id)%cc({:,}, mg_iphi) = (pw(igrid)%w(ix^S, U_br^C_) - 1.0d0)

       ! source term 1 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp1) = 0.0d0

       ! source term 2 = 0
       mg%boxes(id)%cc({1:nc}, mg_itmp2) = 0.0d0

       rhs(ixM^T) = -4.0d0 * dpi * sigma(ixM^T,iigrid) * w_i(ixM^T,3,iigrid)
       mg%boxes(id)%cc({1:nc}, mg_irhs) = rhs(ixM^T)
    end do

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, .True., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       if (res <= gw_br_tol(2)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving U_bri_:', it, mg_it, res, 3
       end if
    end do

    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge U_bri_."
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
       pw(igrid)%w(ix^S, U_br3_) = mg%boxes(id)%cc({:,}, mg_iphi)
    end do
  end subroutine gw_br_solve_U_br3


end module mod_gw_br_U_bri
