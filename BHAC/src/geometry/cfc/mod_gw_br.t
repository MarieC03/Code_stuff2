!> Module for cfc
!< Note that this metric solver support with mod_grhd only!
! fixme: if correct psi as well, meaning we need to do another c2p again?
module mod_gw_br
  use mod_cfc_parameters
  use mod_gw_br_U_br
  use mod_gw_br_R2_br
  use mod_gw_br_R_br
  use mod_gw_br_U_bri
  use mod_gw_br_Iij_Qij
  use mod_full_gr_source
  use mod_eos
  implicit none
  public

  ! public methods
  public :: solve_gw_br

  contains

  subroutine gw_br_activate()
    use mod_multigrid_coupling
    include 'amrvacdef.f'
    integer                      :: n

    if (mype==0) then
       print*,'-----------------------------------------------------------------------------'
       print*,'------------------------GW back reaction activated---------------------------'
       print*,'-----------------------------------------------------------------------------'
    end if

    if (mype == 0) then
       print*,'gw_br_dt_update:', gw_br_dt_update
    endif

    if (coordinate /= Cartesian) call mpistop('Current GW backreaction scheme only for 3D Cart')

    use_multigrid = .True.
  end subroutine gw_br_activate

  logical function gw_br_update_flag()
    include 'amrvacdef.f'
    gw_br_update_flag = .False.

    if (it_start < 5) then
       gw_br_update_flag=.True.
    endif

    if ( t >= gw_br_t_last_update + gw_br_dt_update ) then
       gw_br_update_flag=.True.
    end if

    if (gw_br_use_I3_ij) then
       if (it_start < 3) then
           dt1 = dt
           dt2 = dt_old
           dt3 = dt_old2
           gw_br_t_last_time = t+dt
       else
   
         if (gw_br_update_flag) then
           ! dt between last last update and last update
           dt3 = dt2
           dt2 = dt1
           ! dt between last update and now, but the current value of t is not updated.
           ! the actual t we should use is t+dt, since t is updated after all procedure, but now the time is t + dt
           dt1 = t + dt - gw_br_t_last_time
           gw_br_t_last_time = t + dt
           !if (mype==0) then
           !   write(*,*) 't, dt, dt1, dt2, dt3, gw_br_t_last_time, gw_br_t_last_update'
           !   write(*,*) t, dt, dt1, dt2, dt3,  gw_br_t_last_time, gw_br_t_last_update
           !endif 
           gw_br_t_last_update = t
         endif
       endif
    else
      if (gw_br_update_flag) gw_br_t_last_update = t
      !if (mype==0) then
      !  write(*,*) gw_br_t_last_update, t
      !endif
    endif

  end function gw_br_update_flag

  subroutine solve_gw_br()
    use mod_multigrid_coupling
    use mod_forest
    use mod_imhd_intermediate, only: conformal_transformation, inverse_conformal_trans
    include 'amrvacdef.f'
    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl, ix^D, idim, i, j
    type(tree_node), pointer       :: pnode

    double precision               :: R_br(ixG^T,1:igridstail), R2_br(ixG^T,1:igridstail)
    double precision               :: sigma(ixG^T,1:igridstail), h_00(ixG^T,1:igridstail)
    double precision, allocatable  :: Iij_dot_new(:^D&,:,:), Iij_new(:^D&,:,:), Q3_ij(:^D&,:,:) ! Q_ij^[3]
    double precision               :: array_max_local(1:4), array_max(1:4), alp_min_local, alp_min
    double precision               :: eps, rho, ye, temp, dummy, v2, U, vi(1:3), htot(ixG^T,1:igridstail)
    double precision               :: lfac(ixG^T,1:igridstail), gamma(ixG^T,1:3,1:3,1:igridstail), alp_prime(ixG^T,1:igridstail)
    double precision               :: wi(ixG^T,1:3), psi6_prime(ixG^T), b2(ixG^T)
    double precision               :: w_i(ixG^T,1:3,1:igridstail), g(ixG^T,0:3,0:3), beta(ixG^T,1:3)
    double precision               :: Iij_dot_new_local(1:6), Iij_dot_new_global(1:6)
    double precision               :: Iij_new_local(1:6), Iij_new_global(1:6), Iij_dotdot(1:6), Iij_dotdot_old(1:6)
   
    double precision               :: Q3_ij_test_pt, Q3_ij_test_pt_old, U_bri_local(1:3), U_bri_global(1:3)
    double precision               :: eps_tmp, rho_tmp, prs_tmp, temp_tmp

    double precision, parameter    :: rel_alp = 2.0d-1

    I3ij_new_global(:) = 0.0d0

    ! We need d_i sigma
    call getbc(t,ps,psCoarse)

    sigma = 0.0d0
    h_00  = 0.0d0

    if (gw_br_use_I3_ij) then
      if (use_3rd_d_Iij) then
        allocate(Iij_new(ixG^T,1:6,1:igridstail))
        Iij_new = 0.0d0
      else
        allocate(Iij_dot_new(ixG^T,1:6,1:igridstail))
        Iij_dot_new = 0.0d0
      endif
    endif
    allocate(Q3_ij(ixG^T,1:6,1:igridstail))
    Q3_ij = 0.0d0

    ! initialize the solver following the set eqts of cfc
    if (cfc_redblack) then
       mg%smoother_type = mg_smoother_gsrb
    else
       mg%smoother_type = mg_smoother_gs
    end if

    mg%n_cycle_up   = cfc_n_cycle(1)
    mg%n_cycle_down = cfc_n_cycle(2)

    !-----------------------------------------------------------------------
    ! Step 1: Find sigma := Tmu_mu for the RHS of modified Newtonian potential
    !         U_br_ and later PN corrections
    !-----------------------------------------------------------------------
    ! obtain the sigma: rho_st > sigma > rho
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid = 1, igridstail;  igrid =  igrids(iigrid);
       call set_tmpGlobals(igrid)
       ! will have d_sigma, we need ghostcells as well
       ! case I
       call gw_br_get_sigma_grid(pw(igrid)%w, px(igrid)%x, ixG^LL, ixG^LL, &
                 sigma(ixG^T,iigrid))

       ! case II: using rho_st
       !sigma(ixG^T,iigrid) = pw(igrid)%w(ixG^T,D_) * pw(igrid)%w(ixG^T,psi_metric_)**6 
    end do
    !$OMP END PARALLEL DO
    call gw_br_solve_U_br(sigma(ixG^T,1:igridstail))

    !-----------------------------------------------------------------------
    ! Step 2: Solve for U_br^C (Modified Newtonian momentum???)
    !         and then passing ghostcells before calculation for 
    !         later first/second order derivatives 
    !-----------------------------------------------------------------------
    do iigrid = 1, igridstail;  igrid = igrids(iigrid);
      call set_tmpGlobals(igrid)
      ! including ghostcells

      ! including ghostcells
      if (gw_br_couple_weakly) then
        call imhd_get_intermediate_variables(ixG^LL, ixG^LL, pw(igrid)%w, px(igrid)%x, &
                                             gamma=gamma(ixG^T,1:3,1:3,iigrid), lfac=lfac(ixG^T,iigrid), &
                                              htot=htot(ixG^T,iigrid))
      
      else
        call imhd_get_intermediate_variables(ixG^LL, ixG^LL, pw(igrid)%w, px(igrid)%x, &
                                             gamma=gamma(ixG^T,1:3,1:3,iigrid), lfac=lfac(ixG^T,iigrid), &
                                             b2=b2(ixG^T))
        {do ix^D = ixG^LLIM^D \}
           rho_tmp  = sigma(ix^D,iigrid)
           !rho_tmp  = pw(igrid)%w(ix^D, D_) * pw(igrid)%w(ix^D, psi_metric_)**6
           if (eos_type == tabulated) then
             temp_tmp = pw(igrid)%w(ix^D, T_eps_)

             ! fixme: HOTFIX bound the rho_max for khi proj.
             rho_tmp = max( min( eos_rhomax, rho_tmp), eos_rhomin) 

             call eos_temp_get_all_one_grid(rho_tmp, temp_tmp, pw(igrid)%w(ix^D,ye_),&
                                            eps_tmp, prs=prs_tmp)
             !call eos_get_pressure_one_grid(prs_tmp, rho_tmp, dummy, temp=temp_tmp, ye=pw(igrid)%w(ix^D,ye_))
             !call eos_get_pressure_one_grid(prs_tmp, pw(igrid)%w(ix^D, rho_), dummy, temp=temp_tmp, ye=pw(igrid)%w(ix^D,ye_))

             ! use normal rho to find eps
             !call eos_get_eps_one_grid(prs_tmp, rho_tmp, eps_tmp, temp=temp_tmp, ye=pw(igrid)%w(ix^D,ye_))
             !call eos_get_eps_one_grid(prs_tmp, pw(igrid)%w(ix^D, rho_), eps_tmp, temp=temp_tmp, ye=pw(igrid)%w(ix^D,ye_))
           else
             eps_tmp = pw(igrid)%w(ix^D, T_eps_)
             call eos_get_pressure_one_grid(prs_tmp, rho_tmp, eps_tmp)
           endif
           !htot(ix^D,iigrid) = 1.0d0 + eps_tmp + (prs_tmp + b2(ix^D))/pw(igrid)%w(ix^D, rho_)
           htot(ix^D,iigrid) = 1.0d0 + eps_tmp + (prs_tmp + b2(ix^D))/rho_tmp
        {enddo^D&\}
      endif

      !if (use_h00_coupling) then
      !   alp_prime(ixG^T,iigrid) = dsqrt(pw(igrid)%w(ixG^T, alp_metric_)**2 - pw(igrid)%w(ixG^T,h_00_) )
      !   !alp_prime(ixG^T,iigrid) = pw(igrid)%w(ixG^T,alp_metric_) 
      !else
      !   alp_prime(ixG^T,iigrid) = pw(igrid)%w(ixG^T,alp_metric_) 
      !endif
 
      if (use_index_contract) then
        !{^C&  w_i(ixG^T,^C,iigrid) = pw(igrid)%w(ixG^T, u^C_) / lfac(ixG^T,iigrid) * gamma(ixG^T,^C,^C,iigrid) \}
        ! case II: vi = u^i/u_0 following Lioutas2024
        !{^C& w_i(ixG^T,^C,iigrid) = ( alp_prime(ixG^T,iigrid)*pw(igrid)%w(ixG^T, u^C_)/lfac(ixG^T,iigrid) -&
        !                              pw(igrid)%w(ixG^T, beta_metric^C_) ) * gamma(ixG^T,^C,^C,iigrid) \}
        !{^C& w_i(ixG^T,^C,iigrid) = ( alp_prime(ixG^T,iigrid)*pw(igrid)%w(ixG^T, u^C_)/lfac(ixG^T,iigrid) -&
        !                              pw(igrid)%w(ixG^T, beta_metric^C_) ) \}

        !{^C& w_i(ixG^T,^C,iigrid) = (1.0d0 + htot(ixG^T,iigrid) ) * pw(igrid)%w(ixG^T,u^C_) * gamma(ixG^T,^C,^C,iigrid) \}
        ! case IV: the defintion of htot is without rest mass!!!!!!  (Blanchet 1990)
        !{^C& w_i(ixG^T,^C,iigrid) = (htot(ixG^T,iigrid) ) * pw(igrid)%w(ixG^T,u^C_)  \}
        {^C& w_i(ixG^T,^C,iigrid) = (htot(ixG^T,iigrid) ) * pw(igrid)%w(ixG^T,u^C_) * gamma(ixG^T,^C,^C,iigrid) \}
      else
        !{^C&  w_i(ixG^T,^C,iigrid) = pw(igrid)%w(ixG^T, u^C_) / lfac(ixG^T,iigrid)  \}
        ! case II: vi = u^i/u_0 following Lioutas2024
        !{^C& w_i(ixG^T,^C,iigrid) = alp_prime(ixG^T,iigrid)*pw(igrid)%w(ixG^T, u^C_)/lfac(ixG^T,iigrid) - pw(igrid)%w(ixG^T, beta_metric^C_)\}
        ! case III: u4_i = Wv^i * gamma
        !{^C& w_i(ixG^T,^C,iigrid) = (1.0d0 + htot(ixG^T,iigrid) ) * pw(igrid)%w(ixG^T,u^C_) * gamma(ixG^T,^C,^C,iigrid) \}
        !{^C& w_i(ixG^T,^C,iigrid) = (1.0d0 + htot(ixG^T,iigrid) ) * pw(igrid)%w(ixG^T,u^C_) \}
        ! case IV: the defintion of htot is without rest mass!!!!!!  (Blanchet 1990)
        !{^C& w_i(ixG^T,^C,iigrid) = (htot(ixG^T,iigrid) ) * pw(igrid)%w(ixG^T,u^C_)  \}
        {^C& w_i(ixG^T,^C,iigrid) = (htot(ixG^T,iigrid) ) * pw(igrid)%w(ixG^T,u^C_) * gamma(ixG^T,^C,^C,iigrid) \}
      endif
    enddo

    if (.not. gw_br_use_I3_ij .and. gw_br_include_dU_i) then
      call mpistop('If you want gw_br_include_dU_i, pls add U_brC inside mod_variables')
      ! looping U_br^C
      !do idim = 1, ^NC
        call gw_br_solve_U_br1(sigma(ixG^T,1:igridstail),w_i(ixG^T,1:3,1:igridstail))
        call gw_br_solve_U_br2(sigma(ixG^T,1:igridstail),w_i(ixG^T,1:3,1:igridstail))
        call gw_br_solve_U_br3(sigma(ixG^T,1:igridstail),w_i(ixG^T,1:3,1:igridstail))
      !end do !end do of idim of U_br^C
    endif
    call getbc(t,ps,psCoarse)

    !-----------------------------------------------------------------------
    ! Step 3: Compute the functional of the instantaneous state of matter 
    !         Q_ij^[3] or mass quadropule I_ij^[3]
    !-----------------------------------------------------------------------
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid = 1, igridstail;  igrid = igrids(iigrid);
      call set_tmpGlobals(igrid)

      if (gw_br_use_I3_ij) then
        if (use_3rd_d_Iij) then
            call compute_Iij_grid(ixG^LL, ixM^LL, ps(igrid), Iij_new(ixG^T,1:6,iigrid))
        else
          ! No need for t=0 to get I3_dot and Iij_dot, as we do not have pw_t_dot_old
          if (it_start == 0) then
            Q3_ij(ixM^T,1:6,iigrid) = 0.0d0

            ! store D_st and eps for next time slice, when it does not pass 'compute_Iij_dot'
            pw_t_dot_old(igrid)%w(ixM^T,7)    = pw(igrid)%w(ixM^T,D_) * pw(igrid)%w(ixM^T,psi_metric_)**6

            call dysp_get_lfac(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, lfac(ixG^T,iigrid))

            {do ix^D = ixM^LLIM^D \}
              if (eos_type == tabulated) then
                rho  = pw(igrid)%w(ix^D,rho_)
                temp = pw(igrid)%w(ix^D,T_eps_)
                ye   = pw(igrid)%w(ix^D,ye_)
                call eos_get_eps_one_grid(dummy,rho,eps,temp=temp,ye=ye)
              else
                eps  = pw(igrid)%w(ix^D,T_eps_)
              endif

              {^C& vi(^C) = pw(igrid)%w(ix^D,alp_metric_)*pw(igrid)%w(ix^D, u^C_)/lfac(ix^D,iigrid) - pw(igrid)%w(ix^D, beta_metric^C_)\}
     
              v2 = 0.0d0
              v2 = {^C& vi(^C) * vi(^C) +}
              ! Newtonian potential, which is > U_br
              U = 2.0d0 * (pw(igrid)%w(ix^D,psi_metric_) - 1.0d0)

              ! bracket = (1 + v2/2 - U + eps)
              pw_t_dot_old(igrid)%w(ix^D,8) = 1.0d0 + v2/2.0d0 - U + eps

            {enddo^D&\}

          else
            if (dt1 == 0.0d0) then
               write(*,*) dt1, dt2
               call mpistop('dt1 when computing Iij_dot')
            endif
            ! find the new Iij and store the new w for t_dot_old
            !call compute_Iij_dot_grid(ixG^LL, ixM^LL, ps(igrid), Iij_dot_new(ixG^T,1:6,iigrid), &
            call compute_Iij_dot_grid_Oechslin2007(ixG^LL, ixM^LL, ps(igrid), Iij_dot_new(ixG^T,1:6,iigrid), &
                                                            pw_t_dot_old(igrid)%w(ixG^T,7:12), &
                                                            !pw_t_dot_old(igrid)%w(ixG^T,7:8), &
                                                            dt1, dt2) 
          endif ! endif of it_start = 0
        endif ! endif of use_3rd_d_Iij
      else
        call compute_Q3_ij(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim), &
                           Q3_ij(ixG^T,1:6,iigrid), sigma(ixG^T,iigrid), &
                           gamma(ixG^T,1:3,1:3,iigrid), lfac(ixG^T,iigrid), w_i(ixG^T,1:3,iigrid))
      endif ! endif of gw_br_use_I3_ij
    end do
    !$OMP END PARALLEL DO

    if (gw_br_use_I3_ij) then
      if (use_3rd_d_Iij) then
        Iij_new_local   = 0.0d0
        Iij_new_global  = 0.0d0
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          do i = 1, 6
            Iij_new_local(i)  = Iij_new_local(i) + sum(Iij_dot_new(ixM^T,i,iigrid))
          enddo
        end do
        !$OMP END PARALLEL DO
        call MPI_ALLREDUCE(Iij_new_local,  Iij_new_global, 6, mpi_double_precision, &
              MPI_SUM, icomm, ierrmpi)
        if (dt3 == 0.0d0 .or. it_start < 10) then
          Q3_ij(ixG^T,1:6,1:igridstail) = 0.0d0
        else
          ! We are getting I3dot(t) by backward differencing (need 4 time step)
          if ( (abs(dt2-dt3)/dt2 .le. tol_rel_dt) .and. (abs(dt1-dt2)/dt1 .le. tol_rel_dt) ) then
            !$OMP PARALLEL DO PRIVATE(igrid)
            do iigrid=1,igridstail; igrid=igrids(iigrid);
               do i = 1, 6
                 Iij_dotdot(i)     = 2.0d0*(Iij_new_global(i) - (1.0d0 + dt1/dt2)*Iij_old(i) + &
                              Iij_old2(i)*(dt1/dt2)) / (dt1*(dt1+dt2))

                 Iij_dotdot_old(i) = 2.0d0*(Iij_old(i) - (1.0d0 + dt2/dt3)*Iij_old2(i) + &
                              Iij_old3(i)*(dt2/dt3)) / (dt2*(dt2+dt3))

                 !Q3_ij(ixM^T,i,iigrid) = (Iij_dotdot(i) - Iij_dotdot_old(i))/ ( (dt1 + dt2 + dt3)/3.0d0 )
                 Q3_ij(ixM^T,i,iigrid) = (Iij_dotdot(i) - Iij_dotdot_old(i))/ dt1
               enddo
            end do
            !$OMP END PARALLEL DO
          else
            Q3_ij = 0.0d0
          endif
        endif ! endif of dt3 = 0 or it_start

        Iij_old3(1:6) = Iij_old2(1:6)
        Iij_old2(1:6) = Iij_old(1:6)
        Iij_old(1:6)  = Iij_new_global(1:6)
      else
        Iij_dot_new_local   = 0.0d0
        Iij_dot_new_global  = 0.0d0
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          do i = 1, 6
            Iij_dot_new_local(i)  = Iij_dot_new_local(i) + sum(Iij_dot_new(ixM^T,i,iigrid))
          enddo
        end do
        !$OMP END PARALLEL DO
        call MPI_ALLREDUCE(Iij_dot_new_local,  Iij_dot_new_global, 6, mpi_double_precision, &
              MPI_SUM, icomm, ierrmpi)

        if (dt2 == 0.0d0 .or. it_start < 10) then
          Q3_ij(ixG^T,1:6,1:igridstail) = 0.0d0
        else
          ! We are getting I3dot(t-1) instead of I3dot(t)
          !> I3dot(t-1) = [ Idot(t) - (1+dt1/dt2)* Idot(t-1) + (dt1/dt2)*Idot(t-2) ] / delta_t^2
          if (  abs(dt1-dt2)/dt1 .le. tol_rel_dt ) then
            !$OMP PARALLEL DO PRIVATE(igrid)
            do iigrid=1,igridstail; igrid=igrids(iigrid);
               do i = 1, 6
                 Q3_ij(ixM^T,i,iigrid) = &
                       2.0d0*(Iij_dot_new_global(i) - (1.0d0 + dt1/dt2)*Iij_dot_old(i) + &
                              Iij_dot_old2(i)*(dt1/dt2)) / (dt1*(dt1+dt2))
               enddo
            end do
            !$OMP END PARALLEL DO
          else
            Q3_ij = 0.0d0
          endif
        endif !endif of dt2 = 0 and it_start < 10

           Iij_dot_old2(1:6) = Iij_dot_old(1:6)
           Iij_dot_old(1:6)  = Iij_dot_new_global(1:6)
      endif !endif of use_3rd_d_Iij

    else
        Iij_new_local(:)   = 0.0d0
        Iij_new_global(:)  = 0.0d0
        do iigrid=1,igridstail; igrid=igrids(iigrid);
            Iij_new_local(1)  = Iij_new_local(1) + sum(Q3_ij(ixM^T,1,iigrid))
            Iij_new_local(2)  = Iij_new_local(2) + sum(Q3_ij(ixM^T,2,iigrid))
            Iij_new_local(3)  = Iij_new_local(3) + sum(Q3_ij(ixM^T,3,iigrid))
            Iij_new_local(4)  = Iij_new_local(4) + sum(Q3_ij(ixM^T,4,iigrid))
            Iij_new_local(5)  = Iij_new_local(5) + sum(Q3_ij(ixM^T,5,iigrid))
            Iij_new_local(6)  = Iij_new_local(6) + sum(Q3_ij(ixM^T,6,iigrid))
            ! YOU SHOULD NEVER FORGET HOW TO DO MPISUM --> YOU MUST NEED Iij = Iij + sim(xxx)!!!!!!FUCK
        end do
        call MPI_ALLREDUCE(Iij_new_local,  Iij_new_global, 6, mpi_double_precision, &
              MPI_SUM, icomm, ierrmpi)

        do i = 1, 6
          I3ij_new_global(i) = Iij_new_global(i)
        enddo
    endif !endif of gw_br_use_I3_ij

    !-----------------------------------------------------------------------
    ! Step 4: Solve R_br
    !-----------------------------------------------------------------------

    ! Be careful, becoz at it_start = 0, 1.  We do not have correct t_dot_old, t_dot_old2 --> wrong Qij
    if (it_start > 10 .or. .not. gw_br_use_I3_ij) then

    !if (it_start > 2 .or. .not. gw_br_use_I3_ij) then
      call gw_br_solve_R_br(R_br(ixG^T,1:igridstail), sigma(ixG^T,1:igridstail), gamma(ixG^T,1:3,1:3,1:igridstail))

      !-----------------------------------------------------------------------
      ! Step 4: Using all solutions to find R2_br
      !-----------------------------------------------------------------------
      call gw_br_solve_R2_br(R2_br(ixG^T,1:igridstail), R_br(ixG^T,1:igridstail), &
                             sigma(ixG^T,1:igridstail), gamma(ixG^T,1:3,1:3,1:igridstail), w_i(ixG^T,1:3,1:igridstail))

      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         call set_tmpGlobals(igrid)

         call gw_br_get_h_00_grid(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim), &
                                  h_00(ixG^T,iigrid), R_br(ixG^T,iigrid), R2_br(ixG^T,iigrid), gamma(ixG^T,1:3,1:3,iigrid))
         pw(igrid)%w(ixM^T,h_00_) = h_00(ixM^T,iigrid)

      end do
      !$OMP END PARALLEL DO
         ! GW BR coupling to metric ! FIXME: if we directly change alp_metric inside w, 
         ! it will change the next guess of alp in MG solver --> find weird alp next time
         ! Solution: couple alp in hydro when we need alp
    endif 

    if (gw_br_use_I3_ij) then
      if (use_3rd_d_Iij) then
        deallocate(Iij_new)
      else
        deallocate(Iij_dot_new)
      endif
    endif
    deallocate(Q3_ij)

  end subroutine solve_gw_br

  !> calculate sigma on a single grid
  subroutine gw_br_get_sigma_grid( w, x, ixI^L, ixO^L, sigma)
    use mod_metric
    use mod_imhd_intermediate
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: sigma(ixI^S)
    double precision                :: htot(ixI^S)

    integer                         :: idir, jdir, kdir
    double precision                :: Tmunu(ixI^S,0:3,0:3), Tmu_nu(ixI^S,0:3,0:3), u(ixI^S,0:3), Ptot(ixI^S), bmu(ixI^S, 0:3), alp(ixI^S)
    double precision                :: gamma(ixI^S,1:3,1:3), g(ixI^S,0:3, 0:3), g_munu(ixI^S,0:3,0:3), beta(ixI^S,1:3)

    call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, gamma=gamma(ixI^S,1:3,1:3), &
                                         bmu=bmu(ixI^S,0:3), Ptot=Ptot(ixI^S), htot=htot(ixI^S), u4=u(ixI^S,0:3))

    if (use_h00_coupling) then
      !alp(ixI^S) = dsqrt( w(ixI^S,alp_metric_)**2 - w(ixI^S,h_00_) )
      alp(ixI^S) = w(ixI^S, alp_metric_)
    else
      alp(ixI^S) = w(ixI^S, alp_metric_)
    endif

{^C& beta(ixI^S,^C) = w(ixI^S, beta_metric^C_) \}

    call get_g_up_munu(g(ixI^S,0:3,0:3), alp(ixI^S), beta(ixI^S,1:3), gamma(ixI^S,1:3,1:3), ixI^L, ixO^L)

    ! energy-momentum tensor T^{munu}
    do idir = 0,3
       do jdir = 0,3
          Tmunu(ixI^S,idir,jdir) = w(ixI^S, rho_) * htot(ixI^S) * u(ixI^S, idir) * u(ixI^S, jdir) &
                               + Ptot(ixI^S) * g(ixI^S,idir,jdir) &
                               - bmu(ixI^S, idir) * bmu(ixI^S, jdir)
       enddo
    enddo

    ! case I: T\mu_\mu
    !--> lead to unphysical high delta
    !call get_g_down_munu(g_munu(ixI^S,0:3,0:3), alp(ixI^S), beta(ixI^S,1:3), gamma(ixI^S,1:3,1:3), ixI^L, ixO^L)

    !Tmu_nu = 0.0d0
    !do idir = 0,3
    !   do jdir = 0,3
    !     do kdir = 0,3
    !       Tmu_nu(ixI^S,idir,jdir) = Tmu_nu(ixI^S,idir,jdir) + Tmunu(ixI^S,idir,kdir) * g_munu(ixI^S, kdir, jdir)
    !     enddo
    !   enddo
    !enddo

    !sigma = 0.0d0
    !do idir = 0, 3
    !  sigma(ixI^S) = sigma(ixI^S) + Tmu_nu(ixI^S,idir,idir)
    !enddo

    ! case II: T\mu\mu
    sigma = 0.0d0
    do idir = 0, 3
      sigma(ixI^S) = sigma(ixI^S) + Tmunu(ixI^S,idir,idir)
    enddo 
    

  end subroutine gw_br_get_sigma_grid

  subroutine gw_br_get_h_00_grid(ixI^L, ixO^L, w, x, h_00, R_br, R2_br, gamma)
    use mod_imhd_intermediate
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(in)    :: R_br(ixI^S), R2_br(ixI^S), gamma(ixI^S,1:3,1:3)
    double precision, intent(inout) :: h_00(ixI^S)
    double precision                :: U(ixI^S), d_U(ixI^S,1:3), dU(ixI^S,1:3)
    double precision                :: Q3_ij_x_dU(ixI^S), Q3_ij_local(ixI^S,1:3,1:3), alp_prime_local
    !code test
    double precision                :: h_00_real(ixI^S), beta_sqr(ixI^S), alp_prime_real(ixI^S)
    integer                         :: idir, i ,j, ix^D

      d_U = 0.0d0
  
      U(ixI^S) = w(ixI^S, U_br_)

      do idir = 1, 3
        call partial_d( U(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_U(ixI^S,idir) )
      enddo
    
      if (use_index_contract) then
        do idir = 1, ndir
           dU(ixO^S,idir) = d_U(ixO^S,idir) / gamma(ixO^S,idir,idir)
        enddo
      else
        dU = d_U
      endif
  
      Q3_ij_local(ixO^S,1,1) = I3ij_new_global(1)
      Q3_ij_local(ixO^S,1,2) = I3ij_new_global(2)
      Q3_ij_local(ixO^S,1,3) = I3ij_new_global(3)
      Q3_ij_local(ixO^S,2,2) = I3ij_new_global(4)
      Q3_ij_local(ixO^S,2,3) = I3ij_new_global(5)
      Q3_ij_local(ixO^S,3,3) = I3ij_new_global(6)
  
      Q3_ij_local(ixO^S,2,1) = I3ij_new_global(2)
      Q3_ij_local(ixO^S,3,1) = I3ij_new_global(3)
      Q3_ij_local(ixO^S,3,2) = I3ij_new_global(5)
  
      Q3_ij_x_dU = 0.0d0
      do i = 1, 3
         do j = 1, 3
            Q3_ij_x_dU(ixO^S) = Q3_ij_x_dU(ixO^S) + &
                                Q3_ij_local(ixO^S,i,j) * x(ixO^S, i) * dU(ixO^S,j)
         enddo
      enddo

      h_00(ixO^S) = -4.0d0/5.0d0 * (1.0d0 - 2.0d0 * U(ixO^S)) * (Q3_ij_x_dU(ixO^S) - R_br(ixO^S)) +&
                    4.0d0/5.0d0 * R2_br(ixO^S)


     {do ix^D = ixO^LIM^D \}
       !write(*,*) 4.0d0/5.0d0 * (1.0d0 - 2.0d0 * U(ix^D)) * (Q3_ij_x_dU(ix^D) - w(ix^D,R_br_)), 8.0d0/5.0d0 * w(ix^D,U5_br_)
       !! +ve h_00 code test
       if (h_00(ix^D) < 0.0d0) then
          h_00(ix^D) = -h_00(ix^D)
       endif

       ! get h_00_real after +ve treatment of h_00
       h_00_real(ix^D) = h_00(ix^D)

       if (t < 400.0d0) then
         ! direct bound of h00
         if (h_00(ix^D) > 0.4d-2) then
            h_00(ix^D) = 0.4d-2
         endif
       else
         if (h_00(ix^D) > 1.5d-2) then
            h_00(ix^D) = 1.5d-2
         endif
       endif !only for bounding first peak no crashing


       alp_prime_local = dsqrt (w(ix^D, alp_metric_)**2 - h_00(ix^D) )
       !alp_prime_local = w(ix^D, alp_metric_) - h_00(ix^D)/2.0d0/w(ix^D, alp_metric_)
       if ( alp_prime_local > w(ix^D,alp_metric_)) then
       ! bound for < 0.2 will avoid formation of BH
       !if ( alp_prime_local < 0.2d0 .or. alp_prime_local > w(ix^D,alp_metric_)) then
         !write(*,*) h_00(ix^D), (Q3_ij_x_dU(ix^D) - w(ix^D,R_br_)), Q3_ij_x_dU(ix^D), w(ix^D,R_br_)
         !call mpistop('h_00 < 0')
          h_00(ix^D) = 0.0d0
       endif
    
       beta_sqr(ix^D) = 0.0d0 
       beta_sqr(ix^D) = {^C& gamma(ix^D, ^C, ^C) * w(ix^D, beta_metric^C_) * w(ix^D, beta_metric^C_) +}

       {#IFDEF DELTA
       if ( (w(ix^D, rho_) .gt. 1.0d11 * rho_gf) .and. (beta_sqr(ix^D) .gt. 0.01d0) ) then
         w(ix^D, delta_br_) = (h_00_real(ix^D) - h_00(ix^D))/beta_sqr(ix^D) + 1.0d0
       else
         w(ix^D, delta_br_) = 1.0d0
       endif
       }
     {enddo^D&\}
  end subroutine gw_br_get_h_00_grid

end module mod_gw_br


