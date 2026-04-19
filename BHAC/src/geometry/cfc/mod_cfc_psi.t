module mod_cfc_psi
  use mod_cfc_parameters
  implicit none
  private
!  public
  ! public methods
  public :: cfc_solve_psi
  public :: cfc_check_psi

  contains

  subroutine cfc_psi_solver_init(A2)
    use mod_multigrid_coupling
    use mod_forest
    use mod_metric
    use mod_imhd_intermediate
    include 'amrvacdef.f'

    double precision, intent(in)   :: A2(ixG^T,1:igridstail)
    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: src1(ixG^T)
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;
    mg%operator_type = mg_cfc_psi
    call mg_set_methods(mg)

    select case (coordinate)
    case (cartesian)
       do nb=1, 2*ndim
          select case (typeB(psi_metric_, nb)) 
          case ('cont','noinflow')
             ! outer boundary, use Robin B. C.
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
             mg%bc(nb, mg_iphi)%bc_value = 1.0d0
          case ('symm')
             ! inner boundary, use neumann
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case default
             call mpistop("Error: The boundary conditions of psi/alp are not set correctly.")
          end select
       end do
    case (cylindrical)
       ! Neumann boundary condition for r=0
       mg%bc(1, mg_iphi)%bc_type = mg_bc_neumann
       mg%bc(1, mg_iphi)%bc_value = 0.0d0
       ! Robin boundary condition for r=rmax
       mg%bc(2, mg_iphi)%bc_type = mg_bc_robin
       mg%bc(2, mg_iphi)%bc_value = 1.0d0
       {^NOONED
       ! B.C. for z=zmin
       if( abs(xprobmin2) < smalldouble ) then
          mg%bc(3, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(3, mg_iphi)%bc_value = 0.0d0
       else
          mg%bc(3, mg_iphi)%bc_type = mg_bc_robin
          mg%bc(3, mg_iphi)%bc_value = 1.0d0
       end if
       ! B.C. for z = zmax
       mg%bc(4, mg_iphi)%bc_type = mg_bc_robin
       mg%bc(4, mg_iphi)%bc_value = 1.0d0
       }
       {^IFTHREED
       ! B.C. for phi
       if ( poleB(1,1) ) then
          ! pi-periodic conditions is applied at r=0
          ! nothing to do here
       else
          ! Neumann boundary condition for phi = phi_min and phi_max
          mg%bc(5:6, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(5:6, mg_iphi)%bc_value = 0.0d0
       end if
       }
    case (spherical)
       ! Neumann boundary condition for r=0
       mg%bc(1, mg_iphi)%bc_type = mg_bc_neumann
       mg%bc(1, mg_iphi)%bc_value = 0.0d0
       ! Robin boundary condition for r=rmax
       mg%bc(2, mg_iphi)%bc_type = mg_bc_robin
       mg%bc(2, mg_iphi)%bc_value = 1.0d0
       {^NOONED
       ! Neumann boundary condition for every boundary
       do nb=3, mg_num_neighbors
          mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(nb, mg_iphi)%bc_value = 0.0d0
       end do
       }
       {^IFTHREED
       ! B.C. for phi
       if ( poleB(1,2) .or. poleB(2,2) ) then
          ! pi-periodic conditions is applied at northpole or southpole
          ! nothing to do here
       else
          ! Neumann boundary condition for phi = phi_min and phi_max
          mg%bc(5:6, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(5:6, mg_iphi)%bc_value = 0.0d0
       end if
       }
    end select

    ! copy the data into MG solver, and set up the source terms
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! layer of ghost cells on grid leaves are not included here

       ! initial guess of (psi-1)
       mg%boxes(id)%cc({:,}, mg_iphi) = (pw(igrid)%w(ix^S, psi_metric_) - 1.0d0)

       ! source term 1
       !src1(ixM^T) = 2.0d0 * dpi &
       !      * ( pw(igrid)%w(ixM^T, tau_) + pw(igrid)%w(ixM^T, D_) ) 
       !mg%boxes(id)%cc({1:nc}, mg_itmp1) = src1(ixM^T)
       call imhd_get_tilde_U(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, src1(ixG^T))
       mg%boxes(id)%cc({1:nc}, mg_itmp1) = 2.0d0 * dpi * src1(ixM^T)

       ! source term 2, A2
       mg%boxes(id)%cc({1:nc}, mg_itmp2) = 1.25d-1 * A2(ixM^T, iigrid)

       ! No RHS 
       mg%boxes(id)%cc({1:nc}, mg_irhs) = 0.0d0

       !!write(27,*) "igrid -------------------------",igrid
       !!write(27,*) "psi",(pw(igrid)%w(ix^S, psi_metric_))
       !!write(27,*) "27:A2",A2(ixM^T, iigrid)
       !!write(27,*) "mg_iphi",mg%boxes(id)%cc({:,}, mg_iphi)
       !!write(27,*) "mg_itmp1",mg%boxes(id)%cc({1:nc}, mg_itmp1)
       !!write(27,*) "mg_itmp2",mg%boxes(id)%cc({1:nc}, mg_itmp2)
!!
       !!write(*,*)"file27------------------------------------------------"
       !!write(*,*) "27:igrid -------------------------",igrid
       !!write(*,*) "27:psi",(pw(igrid)%w(ix^S, psi_metric_))
       !!write(*,*) "27:A2",A2(ixM^T, iigrid)
       !!write(*,*) "27:mg_iphi",mg%boxes(id)%cc({:,}, mg_iphi)
       !!write(*,*) "27:mg_itmp1",mg%boxes(id)%cc({1:nc}, mg_itmp1)
       !!write(*,*) "27:mg_itmp2",mg%boxes(id)%cc({1:nc}, mg_itmp2)
    end do
      !  write(27,*) "Cfc_psi_solver_init"
  end subroutine cfc_psi_solver_init

  subroutine cfc_solve_psi(A2)
    use mod_multigrid_coupling
    use mod_forest
    include 'amrvacdef.f'

    double precision, intent(in)   :: A2(ixG^T,1:igridstail)
    integer                        :: id, iigrid, igrid
    double precision               :: res = huge(1.0d0)
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    integer                        :: ix^L

   ! write(*,*)"before cfc psi solver init"

    call cfc_psi_solver_init(A2(ixG^T,1:igridstail))

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, .True., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       !write(*,*) 'solving psi:', mg_it, res
       !!me
       !!write(26,*) 'mg_it:solving psi:', it, mg_it, res
       !!write(26,*) "mg_it:boxes",mg%boxes(id)%cc({:,}, mg_iphi) + 1.0d0
       !!
       if (res <= cfc_tol(2)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving psi:', it, mg_it, res
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge psi."
          write(*,*) "! it=", it,"N = ", mg_it, " Res = ", res
          write(*,*) "! Maybe the cfc_tol is too small or cfc_it_max is too small"
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
       ! Include one layer of ghost cells on grid leaves
       ! the extra layer is needed, as we need to workout d/dr (alp/psi**6) !
       pw(igrid)%w(ix^S, psi_metric_) = mg%boxes(id)%cc({:,}, mg_iphi) + 1.0d0
       !!me
      !! write(26,*) "igrid---------",igrid
      !! write(26,*) "mg:boxes",mg%boxes(id)%cc({:,}, mg_iphi) + 1.0d0
      !! write(26,*) "psi",pw(igrid)%w(ix^S, psi_metric_)
      !! 
      !! write(*,*)"file26-------------------------------"
      !! write(*,*) "26:igrid---------",igrid
      !! write(*,*) "26:mg:boxes",mg%boxes(id)%cc({:,}, mg_iphi) + 1.0d0
      !! write(*,*) "26:psi",pw(igrid)%w(ix^S, psi_metric_)
       !!
    end do
    
   ! write(26,*) "cfc_solve_psi"
  end subroutine cfc_solve_psi

  subroutine cfc_check_psi(res, A2)
    use mod_multigrid_coupling
    use mod_forest
    include 'amrvacdef.f'
    double precision, intent(in)   :: A2(ixG^T,1:igridstail)
    double precision, intent(out)  :: res
    call cfc_psi_solver_init(A2(ixG^T,1:igridstail))
    call mg_max_residual(mg, res)
  end subroutine cfc_check_psi

end module mod_cfc_psi
