
!========================================================================================
!==========                       module mod_post_process                    ==========
!========================================================================================

! Utils of post processing
module mod_post_process
    implicit none

    contains

    ! find global extreme value in pw%w
    double precision function find_global_extreme(take_index, extreme_type)
        include 'amrvacdef.f'
        integer                           :: take_index   ! index where the variable locate
        character(len=3)                  :: extreme_type ! "MIN" or "MAX"
        ! local variables
        double precision                  :: local_value
        integer                           :: iigrid, igrid

        select case (extreme_type)
        case("MIN")
            local_value = 1.0d+99
            do iigrid = 1, igridstail
                igrid = igrids(iigrid)
                local_value = min(local_value, minval(pw(igrid)%w(ixM^T,take_index)) )
            end do
            call MPI_ALLREDUCE(local_value, find_global_extreme, 1, mpi_double_precision, &
                MPI_MIN, icomm, ierrmpi)
        case("MAX")
            local_value = -1.0d+99
            do iigrid = 1, igridstail
                igrid = igrids(iigrid)
                local_value = max(local_value, maxval(pw(igrid)%w(ixM^T,take_index)) )
            end do
            call MPI_ALLREDUCE(local_value, find_global_extreme, 1, mpi_double_precision, &
                MPI_MAX, icomm, ierrmpi)
        case default
            find_global_extreme = 0.d0
            write(*, *) "Error in function find_global_extreme, the extreme_type not supported!"
        end select
    end function find_global_extreme

    subroutine cal_parker_criterion(ixI^L, ixO^L, w, x, parker)
        ! calculate the parker criterion for the magnetic break out using formula
        ! in paper of Parker1966
        use mod_eos
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        include 'amrvacdef.f'
        integer, intent(in)                :: ixI^L, ixO^L
        double precision, intent(in)       :: w(ixI^S, 1:nw)
        double precision, intent(in)       :: x(ixI^S, 1:ndim)
        double precision, intent(out)      :: parker(ixI^S)
        ! local variables
        double precision                   :: b2_tmp(ixI^S)
        double precision                   :: cs2_tmp, ed_tmp, prs_tmp, temp_tmp
        double precision                   :: eps_tmp, gamma_tmp, inv_beta
        integer                            :: ix^D

        parker = 0
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, b2=b2_tmp(ixI^S))
        {do ix^D=ixOmin^D, ixOmax^D \}
            if (eos_type == tabulated) then
                temp_tmp = w(ix^D, T_eps_)
                call eos_get_eps_one_grid(eos_prsmin, w(ix^D, rho_), eps_tmp, temp_tmp, w(ix^D, ye_))
                call eos_get_pressure_one_grid(prs_tmp, w(ix^D, rho_), eps_tmp, temp_tmp, ye=w(ix^D, ye_))
                call eos_get_cs2_one_grid(cs2_tmp, w(ix^D, rho_), eps_tmp, temp_tmp, ye=w(ix^D, ye_))
            else
                eps_tmp = w(ix^D, T_eps_)
                call eos_get_pressure_one_grid(prs_tmp, w(ix^D, rho_), eps_tmp)
                call eos_get_cs2_one_grid(cs2_tmp, w(ix^D, rho_), eps_tmp)
            endif
            ed_tmp = (1+eps_tmp)*w(ix^D, rho_)
            inv_beta = b2_tmp(ix^D)/(2*prs_tmp)
            gamma_tmp = cs2_tmp*(ed_tmp+prs_tmp)/prs_tmp
            parker(ix^D) = gamma_tmp-1-(inv_beta*(1+2*inv_beta))/(2+3*inv_beta)
        {end do \}

    end subroutine cal_parker_criterion


{^IFTHREED
    subroutine MRI_diagnostics(ixI^L, ixO^L, w, x, res, massc_x, massc_y)
        ! See Kenta et al, arXiv: 1710.01311, Eq.(3.2-3.5)
        ! to avoid sigularity in R factor, it is a ratio of volume average on each block
        ! all the phi direction is defined relative to mass center
        use mod_cfc_parameters
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        use mod_metric
        use mod_eos
        include 'amrvacdef.f'
    
        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: w(ixI^S, 1:nw)
        double precision, intent(in)     :: x(ixI^S, 1:ndim)
        double precision, optional, intent(in)  :: massc_x, massc_y
        double precision, intent(out)    :: res(ixI^S, 1:4) ! Qz, Qphi, R, alpha
        ! local
        logical, parameter               :: use_vec = .false.
        double precision                 :: htot(ixI^S), b4_vec(ixI^S, 0:ndir), omega(ixI^S)
        double precision                 :: b4_form(ixI^S, 0:ndir), bsq(ixI^S)
        double precision                 :: delta_xphi, rad, bR, bphi, delta_z
        double precision                 :: sum_brsq, sum_bphisq, val_eps=0.0d0, var_temp, p_lc
        double precision                 :: nc_x, nc_y, cmx=0.0d0, cmy=0.0d0
        integer                          :: ix^D
    
        ! fixme: add other coordinate dependency
        res = 0.0d0
        sum_brsq = 0.0d0
        sum_bphisq = 0.0d0
        if (present(massc_x)) cmx = massc_x
        if (present(massc_y)) cmy = massc_y
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, bmu=b4_vec, b2=bsq, htot=htot)
        select case(coordinate)
        case (cartesian)
            call cartesian_3d_get_angular_omega(ixI^L, ixO^L, w, x, omega, cmx, cmy)
            htot(ixI^S) = htot(ixI^S)*w(ixI^S, rho_)
            if (use_vec) then
                b4_form = b4_vec
            else
                call lower4_dysp(ixI^L, ixO^L, w, x, b4_vec, b4_form)
            endif
            {do ix^D=ixOmin^D, ixOmax^D \}
                if (eos_type == tabulated) then
                    var_temp = w(ix^D, T_eps_)
                    call eos_get_eps_one_grid(0.0d0, w(ix^D, rho_), val_eps, var_temp, w(ix^D, ye_))
                    call eos_get_pressure_one_grid(p_lc, w(ix^D, rho_), val_eps, var_temp, ye=w(ix^D, ye_))
                else
                    call eos_get_pressure_one_grid(p_lc, w(ix^D, rho_), w(ix^D, T_eps_), ye=w(ix^D, ye_))
                endif
                nc_x = x(ix^D, 1)-cmx
                nc_y = x(ix^D, 2)-cmy
                rad = DSQRT(nc_x**2+nc_y**2)
                delta_xphi = (x(ix1+1, ix2, ix3, 1)-x(ix1-1, ix2, ix3, 1))/2. ! r*dphi
                delta_z = (x(ix1, ix2, ix3+1, 3)-x(ix1, ix2, ix3-1, 3))/2.
                if (rad>10*smalldouble) then
                    ! fixme: ambiguity in bphi^2? include the metric term?
                    bphi = (-nc_y*b4_form(ix^D, 1)+nc_x*b4_form(ix^D, 2))/rad**2
                    bR = (nc_x*b4_form(ix^D, 1)+nc_y*b4_form(ix^D, 2))/rad
                    sum_brsq = sum_brsq+bR**2
                    sum_bphisq = sum_bphisq+bphi**2
                    res(ix^D, 1) = 2*dpi*b4_form(ix^D, 3)/(delta_z*omega(ix^D)*DSQRT(htot(ix^D)))
                    res(ix^D, 2) = 2*dpi*bphi/(delta_xphi*omega(ix^D)*DSQRT(htot(ix^D)))
                    !res(ix^D, 4) = -bR*bphi/(bsq(ix^D)/2.) ! maxwell stress over mag pressure
                    res(ix^D, 4) = -bR*bphi/p_lc ! maxwell stress over fluid pressure
                else
                    res(ix^D, 1:4) = 0.0d0
                endif
            {end do \}
            res(ixO^S, 3) = sum_brsq/sum_bphisq
        case default
            call mpistop("MRI_diagnostics Not implemented yet!")
        end select
    
    end subroutine MRI_diagnostics

    
    subroutine imhd_get_em_energy_pt_decomp(ixI^L, ixO^L, w, x, E_em_pol, E_em_tor, massc_x, massc_y)
        ! get the decomposition of the EM energy into poloidal and toroidal
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        use mod_imhd_intermediate, only: imhd_get_Efield_Eulerian
        use mod_metric, only: square3u_dysp
        include 'amrvacdef.f'
        integer, intent(in)                          :: ixI^L, ixO^L
        double precision, intent(in)                 :: w(ixI^S, 1:nw)
        double precision, intent(in)                 :: x(ixI^S, 1:ndim)
        double precision, intent(in)                 :: massc_x, massc_y
        double precision, intent(out)                :: E_em_pol(ixI^S), E_em_tor(ixI^S)
        ! local
        double precision                             :: Br, Bth, Bz, Er, Eth, Ez, rad ! cylindrycal
        double precision                             :: Ef_eul(ixI^S, 1:ndir), x_ref_mc, y_ref_mc
        integer                                      :: ix^D

        call imhd_get_Efield_Eulerian(ixI^L, ixO^L, w, x, Ef_eul)
        {do ix^D=ixOmin^D, ixOmax^D \}
            x_ref_mc = x(ix^D, 1)-massc_x
            y_ref_mc = x(ix^D, 2)-massc_y
            rad = dsqrt(x_ref_mc**2+y_ref_mc**2)
            Bth = -y_ref_mc*w(ix^D, b1_)+x_ref_mc*w(ix^D, b2_) ! phi-component
            Br = x_ref_mc*w(ix^D, b1_)+y_ref_mc*w(ix^D, b2_)   ! r-component
            Bz = w(ix^D, b3_)
            Eth = -y_ref_mc*Ef_eul(ix^D, 1)+x_ref_mc*Ef_eul(ix^D, 2) ! phi-component
            Er = x_ref_mc*Ef_eul(ix^D, 1)+y_ref_mc*Ef_eul(ix^D, 2)   ! r-component
            Ez = Ef_eul(ix^D, 3)
            if (rad>10*smalldouble) then
                E_em_pol(ix^D) = w(ix^D, psi_metric_)**4*((Br**2+Er**2)/rad**2+Bz**2+Ez**2)/2.
                E_em_tor(ix^D) = w(ix^D, psi_metric_)**4*(Bth**2+Eth**2)/(2.*rad**2)
            else
                E_em_pol(ix^D) = 0.0d0
                E_em_tor(ix^D) = 0.0d0
            endif
        {end do \}
    end subroutine imhd_get_em_energy_pt_decomp


    subroutine imhd_get_em_energy_cylind_decomp(ixI^L, ixO^L, w, x, E_em_pol, E_em_tor, E_em_z, massc_x, massc_y)
        ! get the decomposition of the EM energy into poloidal, toroidal, and z component
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        use mod_imhd_intermediate, only: imhd_get_Efield_Eulerian
        use mod_metric, only: square3u_dysp
        include 'amrvacdef.f'
        integer, intent(in)                          :: ixI^L, ixO^L
        double precision, intent(in)                 :: w(ixI^S, 1:nw)
        double precision, intent(in)                 :: x(ixI^S, 1:ndim)
        double precision, intent(in)                 :: massc_x, massc_y
        double precision, intent(out)                :: E_em_pol(ixI^S), E_em_tor(ixI^S), E_em_z(ixI^S)
        ! local
        double precision                             :: Br, Bth, Bz, Er, Eth, Ez, rad ! cylindrycal
        double precision                             :: Ef_eul(ixI^S, 1:ndir), x_ref_mc, y_ref_mc
        integer                                      :: ix^D

        call imhd_get_Efield_Eulerian(ixI^L, ixO^L, w, x, Ef_eul)
        {do ix^D=ixOmin^D, ixOmax^D \}
            x_ref_mc = x(ix^D, 1)-massc_x
            y_ref_mc = x(ix^D, 2)-massc_y
            rad = dsqrt(x_ref_mc**2+y_ref_mc**2)
            Bth = -y_ref_mc*w(ix^D, b1_)+x_ref_mc*w(ix^D, b2_) ! phi-component
            Br = x_ref_mc*w(ix^D, b1_)+y_ref_mc*w(ix^D, b2_)   ! r-component
            Bz = w(ix^D, b3_)
            Eth = -y_ref_mc*Ef_eul(ix^D, 1)+x_ref_mc*Ef_eul(ix^D, 2) ! phi-component
            Er = x_ref_mc*Ef_eul(ix^D, 1)+y_ref_mc*Ef_eul(ix^D, 2)   ! r-component
            Ez = Ef_eul(ix^D, 3)
            if (rad>10*smalldouble) then
                E_em_pol(ix^D) = w(ix^D, psi_metric_)**4*((Br**2+Er**2)/rad**2+Bz**2+Ez**2)/2.
                E_em_tor(ix^D) = w(ix^D, psi_metric_)**4*(Bth**2+Eth**2)/(2.*rad**2)
            else
                E_em_pol(ix^D) = 0.0d0
                E_em_tor(ix^D) = 0.0d0
            endif
            E_em_z(ix^D) = w(ix^D, psi_metric_)**4*(Bz**2+Ez**2)/2.
        {end do \}
    end subroutine imhd_get_em_energy_cylind_decomp

}

    subroutine cartesian_3d_get_toroidal_contraction(ixI^L, ixO^L, w, x, vec3d, phi_contr)
        ! fixme: myabe consider shift of mass center later on
        ! vec^phi*vec_phi using cylindrycal coordinate, assuming CFC metric form
        include 'amrvacdef.f'
        integer, intent(in)                :: ixI^L, ixO^L
        double precision, intent(in)       :: w(ixI^S, 1:nw)
        double precision, intent(in)       :: x(ixI^S, 1:ndim)
        double precision, intent(in)       :: vec3d(ixI^S, 1:3)          ! contravariant 1 vector in 3 space
        double precision, intent(out)      :: phi_contr(ixI^S)
        double precision                   :: vec3d_cylind_psi, rad
        integer                            :: ix^D

        {do ix^D=ixOmin^D, ixOmax^D \}
            rad = dsqrt(x(ix^D, 1)**2+x(ix^D, 2)**2)
            vec3d_cylind_psi = -x(ix^D, 2)*vec3d(ix^D, 1)+x(ix^D, 1)*vec3d(ix^D, 2)! phi-component
            if (rad>10*smalldouble) then
                phi_contr(ix^D) = w(ix^D, psi_metric_)**4*vec3d_cylind_psi**2/rad**2
            else
                phi_contr(ix^D) = 0.
            endif
        {end do \}

    end subroutine cartesian_3d_get_toroidal_contraction


    subroutine cartesian_3d_get_angular_omega(ixI^L, ixO^L, w, x, omega, massc_x, massc_y)
        ! Calculate Omega = u^phi over u^t
        use mod_imhd_intermediate, only: dysp_get_lfac
        include 'amrvacdef.f'
        integer, intent(in)                     :: ixI^L, ixO^L
        double precision, intent(in)            :: w(ixI^S, 1:nw)
        double precision, intent(in)            :: x(ixI^S, 1:ndim)
        double precision, optional, intent(in)  :: massc_x, massc_y
        double precision, intent(out)           :: omega(ixI^S)
        ! local variables
        double precision                   :: u4_vec(ixI^S, 0:3), lfac_temp(ixI^S)
        double precision                   :: vec3d_cylind_psi, rad_sq, cmx=0.0d0, cmy=0.0d0
        integer                            :: ix^D

        call dysp_get_lfac(ixI^L, ixO^L, w, x, lfac_temp(ixI^S))
        u4_vec(ixO^S, 0) = lfac_temp(ixO^S)/w(ixO^S, alp_metric_)
        {^C& u4_vec(ixO^S, ^C) = (w(ixO^S, u^C_)/lfac_temp(ixO^S)-w(ixO^S, beta_metric^C_)/w(ixO^S, alp_metric_))*lfac_temp(ixO^S) \}

        if (present(massc_x)) cmx = massc_x
        if (present(massc_y)) cmy = massc_y
        {do ix^D=ixOmin^D, ixOmax^D \}
            rad_sq = (x(ix^D, 1)-cmx)**2+(x(ix^D, 2)-cmy)**2
            vec3d_cylind_psi = -(x(ix^D, 2)-cmy)*u4_vec(ix^D, 1)+(x(ix^D, 1)-cmx)*u4_vec(ix^D, 2)
            if (rad_sq>10*smalldouble) then
                ! assuming orthonormal coordinate base of vec{e}_phi = partial_phi/r, as if in cylindrical coordinate
                omega(ix^D) = vec3d_cylind_psi/(u4_vec(ix^D, 0)*rad_sq)
            else
                omega(ix^D) = 0.
            endif
        {end do \}
    end subroutine cartesian_3d_get_angular_omega


    subroutine imhd_get_total_em_energy_alternative(ixI^L, ixO^L, w, x, E_em_tot)
        ! Total EM energy: (E^2+B^2)/(2)
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        use mod_imhd_intermediate, only: imhd_get_Efield_Eulerian
        use mod_metric, only: square3u_dysp
        include 'amrvacdef.f'
        integer, intent(in)                          :: ixI^L, ixO^L
        double precision, intent(in)                 :: w(ixI^S, 1:nw)
        double precision, intent(in)                 :: x(ixI^S, 1:ndim)
        double precision, intent(out)                :: E_em_tot(ixI^S)
        double precision, dimension(ixI^S, 1:ndir)   :: Ef_eul           ! contravariant 1 vector in 3 space
        double precision, dimension(ixI^S)           :: Bsq_eul, Esq_eul ! Eulerian frame: B^2, E^2

        call square3u_dysp(ixI^L, ixO^L, w, x, w(ixI^S, b1_:b^NC_), Bsq_eul)
        call imhd_get_Efield_Eulerian(ixI^L, ixO^L, w, x, Ef_eul)
        call square3u_dysp(ixI^L, ixO^L, w, x, Ef_eul, Esq_eul)
        E_em_tot(ixO^S) = (Bsq_eul(ixO^S)+Esq_eul(ixO^S))/2.

    end subroutine imhd_get_total_em_energy_alternative


    subroutine imhd_get_total_em_energy(ixI^L, ixO^L, w, x, E_em_tot)
        ! Total EM energy: (E^2+B^2)/(2)
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        use mod_metric, only: square3u_dysp
        include 'amrvacdef.f'
        integer, intent(in)                          :: ixI^L, ixO^L
        double precision, intent(in)                 :: w(ixI^S, 1:nw)
        double precision, intent(in)                 :: x(ixI^S, 1:ndim)
        double precision, intent(out)                :: E_em_tot(ixI^S)
        double precision, dimension(ixI^S)           :: Bsq_eul(ixI^S), lfac_temp(ixI^S), BdotU(ixI^S) ! Eulerian frame: B^2

        call square3u_dysp(ixI^L, ixO^L, w, x, w(ixI^S, b1_:b^NC_), Bsq_eul)
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, B_dot_v=BdotU(ixI^S), lfac=lfac_temp(ixI^S))
        E_em_tot(ixO^S) = (Bsq_eul(ixO^S)*(2.-1./lfac_temp(ixO^S)**2)-BdotU(ixO^S)**2)/2.

    end subroutine imhd_get_total_em_energy

    subroutine imhd_get_toroidal_em_energy(ixI^L, ixO^L, w, x, E_em_tor)
        ! Use imhd_get_em_energy_pt_decomp instead, which already take into account the shift of the mass center
        ! The toroidal component of EM energy: (E^phi*E_phi+B^phi*B_phi)/(2)
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        use mod_imhd_intermediate, only: imhd_get_Efield_Eulerian
        include 'amrvacdef.f'
        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(out)               :: E_em_tor(ixI^S)
        double precision, dimension(ixI^S, 1:ndir)  :: Ef_eul          ! contravariant
        double precision, dimension(ixI^S)          :: Ephi_contr, Bphi_contr

        call imhd_get_Efield_Eulerian(ixI^L, ixO^L, w, x, Ef_eul)
        call cartesian_3d_get_toroidal_contraction(ixI^L, ixO^L, w, x, Ef_eul, Ephi_contr)
        call cartesian_3d_get_toroidal_contraction(ixI^L, ixO^L, w, x, w(ixI^S, b1_:b^NC_), Bphi_contr)
        E_em_tor(ixO^S) = (Ephi_contr(ixO^S)+Bphi_contr(ixO^S))/2.
    end subroutine imhd_get_toroidal_em_energy


    subroutine imhd_get_poloidal_em_energy(ixI^L, ixO^L, w, x, E_em_pol)
        ! Use imhd_get_em_energy_pt_decomp instead
        ! E_total-E_toroidal
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        include 'amrvacdef.f'
        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(out)               :: E_em_pol(ixI^S)
        double precision, dimension(ixI^S)          :: E_em_total, E_em_toroidal

        call imhd_get_total_em_energy(ixI^L, ixO^L, w, x, E_em_total)
        call imhd_get_toroidal_em_energy(ixI^L, ixO^L, w, x, E_em_toroidal)
        E_em_pol(ixO^S) = E_em_total(ixO^S)-E_em_toroidal(ixO^S)
    end subroutine imhd_get_poloidal_em_energy


    subroutine four_accelaration_euler_observer(ixI^L, ixO^L, w, x, ai)
        ! Calculate a^{i}, See Eq. 4.19 of Eric G.'s book
        use mod_metric, only: raise3_dysp
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(inout)             :: ai(ixI^S, 1:ndir)
        ! local
        double precision                            :: a3_form(ixI^S, 1:ndir)

        ai = 0.0d0
        {^D& call partial_d(dlog(w(ixI^S, alp_metric_)), ixI^L, ixO^L, x, ^D, a3_form(ixI^S, ^D)) \}
        call raise3_dysp(ixI^L, ixO^L, w, x, a3_form, ai)
    end subroutine four_accelaration_euler_observer


    subroutine expansion_scalar(ixI^L, ixO^L, w, w_old, x, acce_euler, expansion, aux_Lambda_scalar, spacial_expansion)
        ! Expansion of four velocity: nabla_{nu} U^{nu}, See Chabanov 21, Eq. (62)-(64)
        use mod_imhd_intermediate, only: dysp_get_lfac
        use mod_metric, only: lower3_dysp
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw), w_old(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: acce_euler(ixI^S, 1:ndir)
        double precision, intent(inout)             :: expansion(ixI^S)
        double precision, intent(inout), optional   :: aux_Lambda_scalar(ixI^S)
        double precision, intent(inout), optional   :: spacial_expansion(ixI^S)
        ! internal variables
        double precision                            :: temp_vec(ixI^S, 1:ndir), lfac_temp(ixI^S), lfac_old(ixI^S)
        double precision                            :: sqrt_gamma(ixI^S), dx_lfac(ixI^S)
        integer                                     :: i

        expansion(ixI^S) = 0.

        call get_sqrt_gamma_hat(x, ixI^L, ixO^L^LADD1, sqrt_gamma(ixI^S))
        sqrt_gamma(ixI^S) = sqrt_gamma(ixI^S) * w(ixI^S, psi_metric_)**6
        {^C& temp_vec(ixI^S, ^C) = sqrt_gamma(ixI^S)*w(ixI^S, u0_+^C) \}
        call divvector(temp_vec(ixI^S, 1:ndir), ixI^L, ixO^L, expansion(ixI^S))
        expansion(ixO^S) = expansion(ixO^S)/sqrt_gamma(ixO^S)
        if (present(spacial_expansion)) spacial_expansion(ixO^S) = expansion(ixO^S) ! assume not the same array after assignment
        call dysp_get_lfac(ixI^L, ixO^L^LADD1, w, x, lfac_temp(ixI^S))
        call dysp_get_lfac(ixI^L, ixO^L^LADD1, w_old, x, lfac_old(ixI^S))
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), temp_vec(ixI^S, 1:ndir))
        if (present(aux_Lambda_scalar)) then
            aux_Lambda_scalar = 0.0d0
            aux_Lambda_scalar(ixO^S) = aux_Lambda_scalar(ixO^S)+(lfac_temp(ixO^S)-lfac_old(ixO^S))/(dt*w(ixO^S, alp_metric_))
            do i=1, ndim
                call partial_d(lfac_temp(ixI^S), ixI^L, ixO^L, x, i, dx_lfac(ixI^S))
                aux_Lambda_scalar(ixO^S) = aux_Lambda_scalar(ixO^S)-dx_lfac(ixO^S)*w(ixO^S, beta_metric1_+i-1)/ &
                    w(ixO^S, alp_metric_)+acce_euler(ixO^S, i)*temp_vec(ixO^S, i)
            enddo
            expansion(ixO^S) = expansion(ixO^S)+aux_Lambda_scalar(ixO^S)
        else
            expansion(ixO^S) = expansion(ixO^S)+(lfac_temp(ixO^S)-lfac_old(ixO^S))/(dt*w(ixO^S, alp_metric_))
            do i=1, ndim
                call partial_d(lfac_temp(ixI^S), ixI^L, ixO^L, x, i, dx_lfac(ixI^S))
                expansion(ixO^S) = expansion(ixO^S)-dx_lfac(ixO^S)*w(ixO^S, beta_metric1_+i-1)/w(ixO^S, alp_metric_)
                expansion(ixO^S) = expansion(ixO^S)+acce_euler(ixO^S, i)*temp_vec(ixO^S, i)
            enddo
        endif
    end subroutine expansion_scalar


    subroutine auxiliary_A_and_lambda(ixI^L, ixO^L, w, w_old, x, acce_euler, aux_A, aux_lambda)
        ! See M.C., 2022, Eq.(46-47): A^i=WU^j*D_j(WU^i), Lambda^i= 1/alpha * (∂t-Lie_beta) WU^i+W hat{a}^i
        use mod_imhd_intermediate, only: dysp_get_lfac
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw), w_old(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: acce_euler(ixI^S, 1:ndir)
        double precision, intent(inout)             :: aux_A(ixI^S, ndir)
        double precision, intent(inout)             :: aux_lambda(ixI^S, ndir)
        ! internal variables
        double precision                            :: temp_derivative_v(ixI^S), temp_derivative_beta(ixI^S)
        double precision                            :: temp_cov_derivative_v(ixI^S), lfac_temp(ixI^S)
        integer                                     :: k, l

        aux_A = 0.0d0
        aux_lambda = 0.0d0
        call dysp_get_lfac(ixI^L, ixO^L, w, x, lfac_temp(ixI^S))
        do k=1, ndim
            do l=1, ndim
                call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), l, k, .true., .true., temp_cov_derivative_v)
                call partial_d(w(ixI^S, u0_+k), ixI^L, ixO^L, x, l, temp_derivative_v(ixI^S))
                call partial_d(w(ixI^S, beta_metric1_+k-1), ixI^L, ixO^L, x, l, temp_derivative_beta(ixI^S))
                aux_A(ixO^S, k) = aux_A(ixO^S, k)+w(ixO^S, u0_+l)*temp_cov_derivative_v(ixO^S)
                aux_lambda(ixO^S, k) = aux_lambda(ixO^S, k)-(temp_derivative_v(ixO^S)*&
                    w(ixO^S, beta_metric1_+l-1)-temp_derivative_beta(ixO^S)*w(ixO^S, u0_+l))/w(ixO^S, alp_metric_)
            enddo
            aux_lambda(ixO^S, k) = aux_lambda(ixO^S, k)+(w(ixO^S, u0_+k)- &
                w_old(ixO^S, u0_+k))/(dt*w(ixO^S, alp_metric_))+lfac_temp(ixO^S)*acce_euler(ixO^S, k)
        enddo
    end subroutine auxiliary_A_and_lambda


    subroutine vorticity_spacial_projection_tensor(ixI^L, ixO^L, w, x, Kij, aux_A, aux_lambda, i, j, vort)
        ! Get pure spacial part of vorticity tensor vort^{ij}, see Eq. (73-74) in M.C.
        use mod_metric, only: lower3_dysp
        use mod_imhd_intermediate, only: dysp_get_lfac
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        integer, intent(in)                         :: i, j
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: Kij(ixI^S, 1:3, 1:3) ! K^ij = A^ij / psi^10
        double precision, intent(in)                :: aux_A(ixI^S, ndir)
        double precision, intent(in)                :: aux_lambda(ixI^S, ndir)
        double precision, intent(inout)             :: vort(ixI^S)
        ! internal variables
        double precision                            :: temp_derivative1(ixI^S), temp_derivative2(ixI^S)
        double precision                            :: temp_vec(ixI^S, 1:ndir), lfac_temp(ixI^S)
        integer                                     :: k

        vort = 0.0d0
        ! vorticity pure spacial part
        call dysp_get_lfac(ixI^L, ixO^L, w, x, lfac_temp(ixI^S))
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), i, j, .true., .false., temp_derivative1)
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), j, i, .true., .false., temp_derivative2)
        vort(ixO^S) = vort(ixO^S)+(temp_derivative1(ixO^S)-temp_derivative2(ixO^S))/2.
        vort(ixO^S) = vort(ixO^S)-(aux_A(ixO^S, i)*w(ixO^S, u0_+j)-aux_A(ixO^S, j)*w(ixO^S, u0_+i))/2.
        ! vorticity of mixed spacial and temporal part
        vort(ixO^S) = vort(ixO^S)-lfac_temp(ixO^S)*(aux_lambda(ixO^S, i)*w(ixO^S, u0_+j) &
            -aux_lambda(ixO^S, j)*w(ixO^S, u0_+i))/2.
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), temp_vec)
        ! vorticity of the curvature part
        do k=1, ndir
            vort(ixO^S) = vort(ixO^S)+lfac_temp(ixO^S)*temp_vec(ixO^S, k)*(Kij(ixO^S, i, k)&
                *w(ixO^S, u0_+j)-w(ixO^S, u0_+i)*Kij(ixO^S, j, k))
        enddo
    end subroutine vorticity_spacial_projection_tensor


    subroutine shear_tensor_spacial_projection(ixI^L, ixO^L, w, x, Kij, aux_A, aux_lambda, aux_Lambda_scalar, spacial_expansion, i, j, shear)
        ! Get pure spacial part shear^{ij}, notice K=0 are used, otherwise we miss some terms, see Eq. (66-69) in M.C.
        use mod_metric, only: lower3_dysp
        use mod_imhd_intermediate, only: dysp_get_lfac
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        integer, intent(in)                         :: i, j 
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: Kij(ixI^S, 1:3, 1:3)
        double precision, intent(in)                :: aux_A(ixI^S, ndir), aux_lambda(ixI^S, ndir)
        double precision, intent(in)                :: spacial_expansion(ixI^S), aux_Lambda_scalar(ixI^S)
        double precision, intent(inout)             :: shear(ixI^S)
        ! internal variables
        double precision                            :: temp_derivative1(ixI^S), temp_derivative2(ixI^S)
        double precision                            :: temp_vec(ixI^S, 1:ndir), lfac_temp(ixI^S)
        double precision                            :: inv_gammaij(ixI^S, 1:3, 1:3)
        integer                                     :: k, idir

        shear = 0.0d0
        ! shear contribution of pure spacial part
        call dysp_get_lfac(ixI^L, ixO^L, w, x, lfac_temp(ixI^S))
        call get_gammainvij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, inv_gammaij(ixI^S,1:3,1:3))                                
        do idir = 1, ndir
           inv_gammaij(ixO^S,idir,idir) = inv_gammaij(ixO^S,idir,idir) / w(ixO^S, psi_metric_)**4
        end do
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), i, j, .true., .false., temp_derivative1)
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), j, i, .true., .false., temp_derivative2)
        shear(ixO^S) = shear(ixO^S)+(temp_derivative1(ixO^S)+temp_derivative2(ixO^S))/2.
        shear(ixO^S) = shear(ixO^S)+(aux_A(ixO^S, i)*w(ixO^S, u0_+j)+aux_A(ixO^S, j)*w(ixO^S, u0_+i))/2. &
            -spacial_expansion(ixO^S)*(inv_gammaij(ixO^S, i, j)+w(ixO^S, u0_+i)*w(ixO^S, u0_+j))/3.
        ! shear contribution of mixed spacial and temporal part
        shear(ixO^S) = shear(ixO^S)+lfac_temp(ixO^S)*(aux_lambda(ixO^S, i)*w(ixO^S, u0_+j)+aux_lambda(ixO^S, j)*w(ixO^S, u0_+i))/2. &
            -aux_Lambda_scalar(ixO^S)*(inv_gammaij(ixO^S, i, j)+w(ixO^S, u0_+i)*w(ixO^S, u0_+j))/3.
        ! shear contribution of the curvature part
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), temp_vec)
        do k=1, ndir
            shear(ixO^S) = shear(ixO^S)-lfac_temp(ixO^S)*temp_vec(ixO^S, k)* &
                (Kij(ixO^S, i, k)*w(ixO^S, u0_+j)+w(ixO^S, u0_+i)*Kij(ixO^S, j, k))
        enddo
        shear(ixO^S) = shear(ixO^S)-lfac_temp(ixO^S)*Kij(ixO^S, i, j)
    end subroutine shear_tensor_spacial_projection


    subroutine source_mag(ixI^L, ixO^L, w, x, Kij, aux_A, aux_lambda, aux_Lambda_scalar, expansion, spacial_expansion, s_mag)
        ! calcualte source term of the magnetic field: s_mag = shear^mu^nu * b_mu * b_nu-expansion*b^2/6
        use mod_metric, only: lower3_dysp
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: Kij(ixI^S, 1:3, 1:3)
        double precision, intent(in)                :: aux_A(ixI^S, ndir), aux_lambda(ixI^S, ndir)
        double precision, intent(in)                :: expansion(ixI^S), spacial_expansion(ixI^S), aux_Lambda_scalar(ixI^S)
        double precision, intent(inout)             :: s_mag(ixI^S)
        ! internal variables
        double precision                            :: shear_spacial_tensor(ixI^S)
        double precision                            :: shear_scalar(ixI^S), shear_space_t_vector(ixI^S, 1:3)
        double precision                            :: temp_vec_U(ixI^S, 1:ndir), temp_vec_B(ixI^S, 1:ndir)
        double precision                            :: gamma_tmp(ixI^S, 1:3, 1:3), lfac_temp(ixI^S)
        double precision                            :: b_3p1_decompose(ixI^S, 0:3), BdotU(ixI^S), bsq_temp(ixI^S)
        integer                                     :: i, j

        s_mag = 0.0d0
        shear_space_t_vector = 0.0d0
        shear_scalar = 0.0d0
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, gamma=gamma_tmp(ixI^S, 1:3, 1:3), &
            B_dot_v=BdotU(ixI^S), lfac=lfac_temp(ixI^S), b2=bsq_temp(ixI^S))
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), temp_vec_U)
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, b1_:b3_), temp_vec_B)
        ! calculate the 3+1 decomposition of b, See 3+1 Book of E. G., Eq. (6.108)
        b_3p1_decompose = 0.0d0
        b_3p1_decompose(ixO^S, 0) = BdotU(ixO^S)*lfac_temp(ixO^S) ! -b^mu*n_mu
        ! spacial projection of b_nu
        do i=1, 3
            b_3p1_decompose(ixO^S, i) = temp_vec_U(ixO^S, i)*BdotU(ixO^S)+temp_vec_B(ixO^S, i)/lfac_temp(ixO^S)
        enddo
        ! get all projections of shear
        do i=1, 3
            do j=1, 3
                call shear_tensor_spacial_projection(ixI^L, ixO^L, w, x, Kij, aux_A, aux_lambda, &
                    aux_Lambda_scalar, spacial_expansion, i, j, shear_spacial_tensor(ixI^S))
                shear_space_t_vector(ixO^S, i) = shear_space_t_vector(ixO^S, i)+temp_vec_U(ixO^S, j)*&
                    shear_spacial_tensor(ixO^S)/lfac_temp(ixO^S)
                shear_scalar(ixO^S) = shear_scalar(ixO^S)+shear_spacial_tensor(ixO^S)*gamma_tmp(ixO^S, i, j)
                s_mag(ixO^S) = s_mag(ixO^S)+shear_spacial_tensor(ixO^S)*b_3p1_decompose(ixO^S, i)*&
                    b_3p1_decompose(ixO^S, j)
            enddo
            s_mag(ixO^S) = s_mag(ixO^S)-2*shear_space_t_vector(ixO^S, i)*b_3p1_decompose(ixO^S, i)*b_3p1_decompose(ixO^S, 0)
        enddo
        s_mag(ixO^S) = s_mag(ixO^S)+shear_scalar(ixO^S)*b_3p1_decompose(ixO^S, 0)**2-expansion(ixO^S)*bsq_temp(ixO^S)/6.
    end subroutine source_mag


    subroutine m1_outflow_flux(ixI^L, ixO^L, w, x, Outflux, speciesKSP)
        !> routine to call calculate_m1_outflow_flux from mod_m1_outflow.t
        !> calculate fluxes for the m1-outflow
        {#IFDEF M1
        use mod_m1_outflow
        }
        include 'amrvacdef.f'
        integer, intent(in)                         :: ixI^L, ixO^L, speciesKSP
        double precision, intent(inout)             :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(inout)             :: Outflux(ixI^S, 1:ndim)
        outflux = 0.0d0
        {#IFDEF M1
        call calculate_m1_outflow_flux(ixI^L,ixO^L, w, x, outflux, speciesKSP)
        }
    end subroutine m1_outflow_flux


    subroutine poynting_flux(ixI^L, ixO^L, w, x, poynt)
        ! See equation (11) of Bernard J. Kelly, et al 2017, S^i = alpha*(b^2*u^i*u_0-b^i*b_0)
        use mod_metric, only: lower3_dysp
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(inout)             :: poynt(ixI^S, 1:ndim)
        ! internal variables
        double precision                            :: b_0(ixI^S), u_0(ixI^S), bsq(ixI^S), lfac_temp(ixI^S)
        double precision                            :: b_4vec(ixI^S, 0:ndir), u4_vec(ixI^S, 0:ndir)
        double precision                            :: BdotU(ixI^S), B_3d_1form(ixI^S, 1:3), WU_3d_1form(ixI^S, 1:3)
        integer                                     :: i

        poynt = 0.0d0
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, b2=bsq, B_dot_v=BdotU, &
                bmu=b_4vec, u4=u4_vec, lfac=lfac_temp(ixI^S))
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form)
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, b1_:b3_), B_3d_1form)
        ! calculate b_0 and u_0
        b_0(ixO^S) = -BdotU(ixO^S)*w(ixO^S, alp_metric_)*lfac_temp(ixO^S)
        u_0(ixO^S) = -w(ixO^S, alp_metric_)*lfac_temp(ixO^S)
        do i=1, ndir
            b_0(ixO^S) = b_0(ixO^S)+w(ixO^S, beta_metric1_+i-1)*B_3d_1form(ixO^S, i)/lfac_temp(ixO^S)+&
                    BdotU(ixO^S)*WU_3d_1form(ixO^S, i)*w(ixO^S, beta_metric1_+i-1)
            u_0(ixO^S) = u_0(ixO^S)+WU_3d_1form(ixO^S, i)*w(ixO^S, beta_metric1_+i-1)
        enddo
        ! calculate poynting flux
        do i=1, ndim
            poynt(ixO^S, i) = w(ixO^S, alp_metric_)*bsq(ixO^S)*u4_vec(ixO^S, i)*u_0(ixO^S)-&
                    w(ixO^S, alp_metric_)*b_4vec(ixO^S, i)*b_0(ixO^S)
        enddo
    end subroutine poynting_flux

    subroutine poynting_r(ixI^L, ixO^L, w, x, poynt_r)
        include 'amrvacdef.f'
        ! project poynting flux in the r direction
        ! fixme: result looks strange
        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(inout)             :: poynt_r(ixI^S)
        ! local variables
        double precision                            :: poynt(ixI^S, 1:ndim)
        double precision                            :: aux_r(ixI^S)

        poynt_r = 0.0d0
        call poynting_flux(ixI^L, ixO^L, w, x, poynt)
        aux_r(ixO^S) = dsqrt({^D& x(ixO^S, ^D)**2+})
        poynt_r(ixO^S) = ({^D& poynt(ixO^S, ^D)*x(ixO^S, ^D)+})/aux_r(ixO^S)
    end subroutine poynting_r

    subroutine imhd_cal_current_mag(ixI^L, ixO^L, w, x, magJ)
        ! See M.C., 2022 [Crustal magnetic fields], Eq.(1) in Appendix, magnitude of the electric current
        use mod_metric, only: square3u_dysp, lower3_dysp
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(inout)             :: magJ(ixI^S)
        ! internal variables
        double precision                            :: currentJ(ixI^S, 1:ndir), temp_alphaB_down(ixI^S, 1:ndir)! Ji, alpha*B_i
        double precision                            :: temp_derivative1(ixI^S), temp_derivative2(ixI^S)

        currentJ = 0.0d0
        magJ = 0.0d0
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, b1_:b3_), temp_alphaB_down(ixI^S, 1:ndir))
        {^C& temp_alphaB_down(ixO^S, ^C) = temp_alphaB_down(ixO^S, ^C)*w(ixO^S, alp_metric_) \}
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, temp_alphaB_down(ixI^S, 1:ndir), 2, 3, .false., .true., temp_derivative1)
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, temp_alphaB_down(ixI^S, 1:ndir), 3, 2, .false., .true., temp_derivative2)
        currentJ(ixO^S, 1) = (temp_derivative1(ixO^S)-temp_derivative2(ixO^S))/w(ixO^S, alp_metric_)
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, temp_alphaB_down(ixI^S, 1:ndir), 3, 1, .false., .true., temp_derivative1)
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, temp_alphaB_down(ixI^S, 1:ndir), 1, 3, .false., .true., temp_derivative2)
        currentJ(ixO^S, 2) = (temp_derivative1(ixO^S)-temp_derivative2(ixO^S))/w(ixO^S, alp_metric_)
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, temp_alphaB_down(ixI^S, 1:ndir), 1, 2, .false., .true., temp_derivative1)
        call CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, temp_alphaB_down(ixI^S, 1:ndir), 2, 1, .false., .true., temp_derivative2)
        currentJ(ixO^S, 3) = (temp_derivative1(ixO^S)-temp_derivative2(ixO^S))/w(ixO^S, alp_metric_)
        call square3u_dysp(ixI^L, ixO^L, w, x, currentJ(ixI^S, 1:3), magJ(ixI^S))
        magJ(ixO^S) = dsqrt(magJ(ixO^S))
    end subroutine imhd_cal_current_mag


    subroutine kinematic_acceleration(ixI^L, ixO^L, w, x, Kij, aux_A, aux_lambda, kinacc)
        ! Calcualte 3+1 decomposition of kinematic accelaration: (a^0, a^i_{perp})
        use mod_metric, only: lower3_dysp
        use mod_imhd_intermediate, only: dysp_get_lfac
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: Kij(ixI^S, 1:3, 1:3)
        double precision, intent(in)                :: aux_A(ixI^S, 1:ndir), aux_lambda(ixI^S, 1:ndir)
        double precision, intent(inout)             :: kinacc(ixI^S, 0:3)
        ! local variables
        integer                                     :: i, j
        double precision                            :: lfac_temp(ixI^S), WU_3d_1form(ixI^S, 1:3) ! WU_i

        kinacc = 0.0d0
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form(ixI^S, 1:3))
        call dysp_get_lfac(ixI^L, ixO^L, w, x, lfac_temp(ixI^S))
        do i=1, ndir
            kinacc(ixO^S, i) = kinacc(ixO^S, i)+aux_A(ixO^S, i)+lfac_temp(ixO^S)*aux_lambda(ixO^S, i)
            do j=1, ndir
                kinacc(ixO^S, i) = kinacc(ixO^S, i)-2*lfac_temp(ixO^S)*WU_3d_1form(ixO^S, j)*Kij(ixO^S, i, j)
            enddo
            kinacc(ixO^S, 0) = kinacc(ixO^S, 0)+WU_3d_1form(ixO^S, i)*kinacc(ixO^S, i)/lfac_temp(ixO^S)
        enddo
    end subroutine kinematic_acceleration


    subroutine cal_adm_mass_volume_contribution(ixI^L, ixO^L, w, x, mm)
        ! Calculate the volume integration form of AMD mass in CFC using Eq. 8.52 of E.G.'s book
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        use mod_metric
        include 'amrvacdef.f'
        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: w(ixI^S, 1:nw)
        double precision, intent(in)     :: x(ixI^S, 1:ndim)
        double precision, intent(out)    :: mm(ixI^S)
        ! local
        double precision                 :: K_ij(ixI^S, 1:ndir, 1:ndir), gammainv(ixI^S,1:ndir,1:ndir)
        integer                          :: i, j, k, l, idir

        mm = 0.0d0
        call get_gammainvij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gammainv(ixI^S,1:ndir,1:ndir))
        do idir = 1, ndir
           gammainv(ixO^S,idir,idir) = gammainv(ixO^S,idir,idir) / w(ixO^S, psi_metric_)**4
        end do
        call get_K_ij_dysp(ixI^L, ixO^L, w, x, K_ij)
        do i=1, ndir
            do j=1, ndir
                do k=1, ndir
                    do l=1, ndir
                        mm(ixO^S) = mm(ixO^S)+K_ij(ixO^S, i, j)*K_ij(ixO^S, k, l)*&
                            gammainv(ixO^S, k, j)*gammainv(ixO^S, l, i)
                    enddo
                enddo
            enddo
        enddo
        mm(ixO^S) = (mm(ixO^S)/(16.*dpi)+w(ixO^S, tau_)+w(ixO^S, d_))/w(ixO^S, psi_metric_)
    end subroutine cal_adm_mass_volume_contribution


    subroutine total_kinetic_energy(ixI^L, ixO^L, w, x, kin)
        ! See Shibata21 (2109.08732) Eq. 2.25
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        use mod_metric, only: lower3_dysp
        include 'amrvacdef.f'

        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: w(ixI^S, 1:nw)
        double precision, intent(in)     :: x(ixI^S, 1:ndim)
        double precision, intent(out)    :: kin(ixI^S)
        ! local
        double precision                 :: WU_3d_1form(ixI^S, 1:ndir), h_thermal(ixI^S)
        double precision                 :: lfac_temp(ixI^S), wudotbeta(ixI^S)
        integer                          :: i

        kin = 0.0
        wudotbeta = 0.0
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form(ixI^S, 1:ndir))
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, lfac=lfac_temp, h_th=h_thermal)
        do i=1, 3
            wudotbeta(ixO^S) = wudotbeta(ixO^S)+WU_3d_1form(ixO^S, i)*w(ixO^S, beta_metric1_+i-1)
        enddo
        kin(ixO^S) = w(ixO^S, alp_metric_)*(lfac_temp(ixO^S)**2-1)-lfac_temp(ixO^S)*wudotbeta(ixO^S)
        kin(ixO^S) = w(ixO^S, rho_)*h_thermal(ixO^S)*kin(ixO^S)/2.
    end subroutine total_kinetic_energy


    subroutine total_z_kinetic_energy(ixI^L, ixO^L, w, x, kin_z)
        ! See Shibata21 (2109.08732) Eq. 2.25
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        use mod_metric, only: lower3_dysp
        include 'amrvacdef.f'

        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: w(ixI^S, 1:nw)
        double precision, intent(in)     :: x(ixI^S, 1:ndim)
        double precision, intent(out)    :: kin_z(ixI^S)
        ! local
        double precision                 :: WU_3d_1form(ixI^S, 1:ndir), h_thermal(ixI^S)
        double precision                 :: lfac_temp(ixI^S)

        kin_z = 0.0
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form(ixI^S, 1:ndir))
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, lfac=lfac_temp, h_th=h_thermal)
        kin_z(ixO^S) = w(ixO^S, alp_metric_)*w(ixO^S, u1_+ndim-1)-lfac_temp(ixO^S)*w(ixO^S, beta_metric1_+ndim-1)
        kin_z(ixO^S) = kin_z(ixO^S)*WU_3d_1form(ixO^S, ndim)*w(ixO^S, rho_)*h_thermal(ixO^S)/2.
    end subroutine total_z_kinetic_energy

    subroutine cal_turbulent_with_block_mean(ixI^L, ixO^L, w, x, turb_part, mean_part, var_idx)
        ! calculate the turbulent part of the vectors using average value from one block
        include 'amrvacdef.f'
        integer, intent(in)                :: ixI^L, ixO^L, var_idx
        double precision, intent(in)       :: w(ixI^S, 1:nw)
        double precision, intent(in)       :: x(ixI^S, 1:ndim)
        double precision, intent(inout)    :: turb_part(ixI^S)
        double precision, intent(inout)    :: mean_part(ixI^S)
        ! local
        integer                            :: num_value

        mean_part = 0.0d0
        turb_part = 0.0d0
        num_value = {(ixOmax^D-ixOmin^D+1)*}
        mean_part(ixO^S) = sum(w(ixO^S, var_idx))/num_value
        turb_part(ixO^S) = w(ixO^S, var_idx)-mean_part(ixO^S)

    end subroutine cal_turbulent_with_block_mean


    subroutine cartesian_3d_rotational_kinetic_energy(ixI^L, ixO^L, w, x, rot_kin, massc_x, massc_y)
        ! See Luciano's book Eq.12.48, including the EM angular
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        ! Note: not very clear this will include coordinate or not
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        use mod_metric, only: lower3_dysp
        include 'amrvacdef.f'

        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: w(ixI^S, 1:nw)
        double precision, intent(in)     :: x(ixI^S, 1:ndim)
        double precision, optional, intent(in)  :: massc_x, massc_y
        double precision, intent(out)    :: rot_kin(ixI^S)
        ! local
        double precision                 :: h_thermal(ixI^S), omega(ixI^S), rad
        double precision                 :: lfac_temp(ixI^S), WUdotB(ixI^S)
        double precision                 :: B_3d_1form(ixI^S, 1:ndir), WU_3d_1form(ixI^S, 1:ndir)
        double precision                 :: b_x(ixI^S), b_y(ixI^S), b_phi(ixI^S), u4_phi(ixI^S)
        double precision                 :: nc_x, nc_y, cmx=0.0d0, cmy=0.0d0
        integer                          :: i, ix^D

        rot_kin = 0.0
        WUdotB = 0.0
        if (present(massc_x)) cmx = massc_x
        if (present(massc_y)) cmy = massc_y
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form(ixI^S, 1:ndir))
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, b1_:b3_), B_3d_1form(ixI^S, 1:ndir))
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, lfac=lfac_temp, h_th=h_thermal)
        do i=1, ndir
            WUdotB(ixO^S) = WUdotB(ixO^S)+WU_3d_1form(ixO^S, i)*w(ixO^S, b1_+i-1)
        enddo
        b_x(ixO^S) = (B_3d_1form(ixO^S, 1)+WUdotB(ixO^S)*WU_3d_1form(ixO^S, 1))/lfac_temp(ixO^S)
        b_y(ixO^S) = (B_3d_1form(ixO^S, 2)+WUdotB(ixO^S)*WU_3d_1form(ixO^S, 2))/lfac_temp(ixO^S)
        {do ix^D=ixOmin^D, ixOmax^D \}
            nc_x = x(ix^D, 1)-cmx
            nc_y = x(ix^D, 2)-cmy
            rad = DSQRT(nc_x**2+nc_y**2)
            if (rad>10*smalldouble) then
                ! fixme: ambiguity in bphi^2? include the metric term?
                b_phi(ix^D) = (-nc_y*b_x(ix^D)+nc_x*b_y(ix^D))/rad**2
                u4_phi(ix^D) = (-nc_y*WU_3d_1form(ix^D, 1)+nc_x*WU_3d_1form(ix^D, 2))/rad**2
            else
                b_phi(ix^D) = 0.0d0
                u4_phi(ix^D) = 0.0d0
            endif
        {end do \}
        call cartesian_3d_get_angular_omega(ixI^L, ixO^L, w, x, omega, cmx, cmy)
        rot_kin(ixO^S) = (lfac_temp(ixO^S)*w(ixO^S, rho_)*h_thermal(ixO^S)*u4_phi(ixO^S)-&
            WUdotB(ixO^S)*b_phi(ixO^S))*omega(ixO^S)/2.
    end subroutine cartesian_3d_rotational_kinetic_energy


    subroutine cartesian_3d_total_angular_momentum(ixI^L, ixO^L, w, x, ang_mom, massc_x, massc_y)
        ! See Luciano's book Eq.12.47, including the EM angular
        ! Note: don't forget to multiply the sqrt(gamma) factor in the integration stage
        ! Note: the T(n, phi-down) projection differs with projection on phi direction by a factor of rsquare
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        use mod_metric, only: lower3_dysp
        include 'amrvacdef.f'

        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: w(ixI^S, 1:nw)
        double precision, intent(in)     :: x(ixI^S, 1:ndim)
        double precision, optional, intent(in)  :: massc_x, massc_y
        double precision, intent(out)    :: ang_mom(ixI^S)
        ! local
        double precision                 :: h_total(ixI^S), p_tot(ixI^S)
        double precision                 :: lfac_temp(ixI^S), WUdotB(ixI^S)
        double precision                 :: B_3d_1form(ixI^S, 1:ndir), WU_3d_1form(ixI^S, 1:ndir)
        double precision                 :: b_x(ixI^S), b_y(ixI^S), b_phi(ixI^S), u4_phi(ixI^S)
        double precision                 :: beta_3d_1form(ixI^S, 1:3), beta_phi(ixI^S)
        double precision                 :: nc_x, nc_y, cmx=0.0d0, cmy=0.0d0
        integer                          :: i, ix^D

        ang_mom = 0.0
        WUdotB = 0.0
        if (present(massc_x)) cmx = massc_x
        if (present(massc_y)) cmy = massc_y
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form(ixI^S, 1:ndir))
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, b1_:b3_), B_3d_1form(ixI^S, 1:ndir))
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, beta_metric1_:beta_metric3_), beta_3d_1form(ixI^S, 1:ndir))
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, lfac=lfac_temp, Ptot=p_tot, htot=h_total)
        do i=1, ndir
            WUdotB(ixO^S) = WUdotB(ixO^S)+WU_3d_1form(ixO^S, i)*w(ixO^S, b1_+i-1)
        enddo
        b_x(ixO^S) = (B_3d_1form(ixO^S, 1)+WUdotB(ixO^S)*WU_3d_1form(ixO^S, 1))/lfac_temp(ixO^S)
        b_y(ixO^S) = (B_3d_1form(ixO^S, 2)+WUdotB(ixO^S)*WU_3d_1form(ixO^S, 2))/lfac_temp(ixO^S)
        {do ix^D=ixOmin^D, ixOmax^D \}
            nc_x = x(ix^D, 1)-cmx
            nc_y = x(ix^D, 2)-cmy
            b_phi(ix^D) = (-nc_y*b_x(ix^D)+nc_x*b_y(ix^D))
            u4_phi(ix^D) = (-nc_y*WU_3d_1form(ix^D, 1)+nc_x*WU_3d_1form(ix^D, 2))
            beta_phi(ix^D) = (-nc_y*beta_3d_1form(ix^D, 1)+nc_x*beta_3d_1form(ix^D, 2))
        {end do \}
        ang_mom(ixO^S) = (lfac_temp(ixO^S)*w(ixO^S, rho_)*h_total(ixO^S)*u4_phi(ixO^S)+&
            p_tot(ixO^S)*beta_phi(ixO^S)-WUdotB(ixO^S)*b_phi(ixO^S))
    end subroutine cartesian_3d_total_angular_momentum


    subroutine cartesian_3d_specific_angular_momentum(ixI^L, ixO^L, w, x, specific_am, massc_x, massc_y)
        ! See E.G.'s book rotating NS E.Q.(3.92)
        use mod_imhd_intermediate, only: dysp_get_lfac
        use mod_metric, only: lower3_dysp
        include 'amrvacdef.f'

        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: w(ixI^S, 1:nw)
        double precision, intent(in)     :: x(ixI^S, 1:ndim)
        double precision, optional, intent(in)  :: massc_x, massc_y
        double precision, intent(out)    :: specific_am(ixI^S)
        ! local
        double precision                 :: rad_sq, u4_phi(ixI^S), u4_t(ixI^S), lfac_temp(ixI^S)
        double precision                 :: WU_3d_1form(ixI^S, 1:ndir)
        double precision                 :: nc_x, nc_y, cmx=0.0d0, cmy=0.0d0
        integer                          :: i, ix^D

        specific_am = 0.0
        call dysp_get_lfac(ixI^L, ixO^L, w, x, lfac_temp(ixI^S))
        if (present(massc_x)) cmx = massc_x
        if (present(massc_y)) cmy = massc_y
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form(ixI^S, 1:ndir))
        u4_t(ixO^S) = w(ixO^S, psi_metric_)**4*({^D& w(ixO^S, u^D_)*w(ixO^S, beta_metric^D_)+})&
            -w(ixO^S, alp_metric_)*lfac_temp(ixO^S)
        {do ix^D=ixOmin^D, ixOmax^D \}
            nc_x = x(ix^D, 1)-cmx
            nc_y = x(ix^D, 2)-cmy
            rad_sq = nc_x**2+nc_y**2
            if (rad_sq>10*smalldouble) then
                ! always consider the coordinate base vector, not the orthernormal one
                u4_phi(ix^D) = (-nc_y*WU_3d_1form(ix^D, 1)+nc_x*WU_3d_1form(ix^D, 2))!/rad_sq
            else
                u4_phi(ix^D) = 0.0d0
            endif
        {end do \}
        specific_am(ixO^S) = -u4_phi(ixO^S)/u4_t(ixO^S)
    end subroutine cartesian_3d_specific_angular_momentum


    subroutine CFC_partiald_3d1vecform(ixI^L, ixO^L, w, x, vecform, idir, jdir, is_vec, is_cov, derivative)
        ! Calculate Di(vj), vj can be vector or form (controlled by is_vec), Di can be covariant or
        ! contravariant (controlled by is_cov), see Eq. 7.32 of Eric, G.'s book, assuming diagonal gamma_hat
        include 'amrvacdef.f'

        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision, intent(in)    :: w(ixI^S, 1:nw)
        double precision, intent(in)    :: vecform(ixI^S, ndir) ! 3d 1vector or 1form, controlled by is_vec
        integer, intent(in)             :: idir, jdir
        logical, intent(in)             :: is_vec, is_cov
        double precision, intent(out)   :: derivative(ixI^S)
        ! local
        integer                         :: k
        double precision                :: dLln_psi(ixI^S, ndir) ! tilde{D}_i ln(psi)
        double precision                :: contract_vdpsi(ixI^S) ! v_i tilde{D}^i ln(psi)
        double precision                :: gamma_hat(ixI^S,1:ndir,1:ndir)

        derivative = 0.0d0
        dLln_psi = 0.0d0
        contract_vdpsi = 0.0d0
        call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_hat(ixI^S,1:ndir,1:ndir))
        call partial_d(vecform(ixI^S, jdir), ixI^L, ixO^L, x, idir, derivative)
        do k=1, ndir
            call partial_d(dlog(w(ixI^S, psi_metric_)), ixI^L, ixO^L, x, k, dLln_psi(ixI^S, k))
            if (is_vec) then
                contract_vdpsi(ixO^S) = contract_vdpsi(ixO^S)+vecform(ixO^S, k)*dLln_psi(ixO^S, k)
            else
                contract_vdpsi(ixO^S) = contract_vdpsi(ixO^S)+vecform(ixO^S, k)*dLln_psi(ixO^S, k)/gamma_hat(ixO^S, k, k)
            endif
        enddo
        if (is_vec) then
            derivative(ixO^S) = derivative(ixO^S)+2*(contract_vdpsi(ixO^S)*kr(idir, jdir)+&
                vecform(ixO^S, jdir)*dLln_psi(ixO^S, idir)-(dLln_psi(ixO^S, jdir)/gamma_hat(ixO^S, jdir, jdir))*&
                (vecform(ixO^S, idir)*gamma_hat(ixO^S, idir, idir)))
        else
            derivative(ixO^S) = derivative(ixO^S)-2*(vecform(ixO^S, jdir)*dLln_psi(ixO^S, idir)+&
                vecform(ixO^S, idir)*dLln_psi(ixO^S, jdir)-contract_vdpsi(ixO^S)*gamma_hat(ixO^S, idir, jdir))
        endif
        if (.not. is_cov) derivative(ixO^S) = derivative(ixO^S)/(gamma_hat(ixO^S, idir, idir)*w(ixO^S, psi_metric_)**4)
    end subroutine CFC_partiald_3d1vecform


    subroutine reconstruct_4d_2vector_shear_tensor(ixI^L, ixO^L, w, x, shearij, shearmunu)
        ! Calculate shear^munu from shear^ij
        use mod_metric, only: lower3_dysp
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: shearij(ixI^S, 1:3, 1:3)
        double precision, intent(inout)             :: shearmunu(ixI^S, 0:3, 0:3)
        ! local variables
        double precision                            :: euler_4d_1vecctor(ixI^S, 0:3), gamma_tmp(ixI^S, 1:3, 1:3)
        double precision                            :: shear_proj_3vec(ixI^S, 0:3), shear_spacial_trace(ixI^S)
        double precision                            :: lfac_temp(ixI^S), WU_3d_1form(ixI^S, 1:3) ! WU_i
        integer                                     :: i, j

        ! init and get auxiliary variables
        shearmunu = 0.0d0
        shear_proj_3vec = 0.0d0
        shear_spacial_trace = 0.0d0
        shearmunu(ixO^S, 1:3, 1:3) = shearij(ixO^S, 1:3, 1:3)
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form(ixI^S, 1:3))
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, lfac=lfac_temp, gamma=gamma_tmp)
        euler_4d_1vecctor(ixO^S, 0) = 1./w(ixO^S, alp_metric_)
        ! calculate vec projection and trace
        do i=1, ndir
            euler_4d_1vecctor(ixO^S, i) = -w(ixO^S, beta_metric1_+i-1)/w(ixO^S, alp_metric_)
            do j=1, ndir
                shear_proj_3vec(ixO^S, i) = shear_proj_3vec(ixO^S, i)+&
                    shearij(ixO^S, i, j)*WU_3d_1form(ixO^S, j)/lfac_temp(ixO^S)
                !shear_spacial_trace(ixO^S) = shear_spacial_trace(ixO^S)+shearij(ixO^S, i, j)*gamma_tmp(ixO^S, i, j)
            enddo
        enddo
        do i=1, ndir
            shear_spacial_trace(ixO^S) = shear_spacial_trace(ixO^S)+shear_proj_3vec(ixO^S, i)* &
                WU_3d_1form(ixO^S, i)/lfac_temp(ixO^S)
        enddo
        ! calculate shearmunu
        do i=0, ndir
            do j=0, ndir
                shearmunu(ixO^S, i, j) = shearmunu(ixO^S, i, j)+&
                    shear_spacial_trace(ixO^S)*euler_4d_1vecctor(ixO^S, i)*euler_4d_1vecctor(ixO^S, j)+ &
                    (shear_proj_3vec(ixO^S, i)*euler_4d_1vecctor(ixO^S, j)+shear_proj_3vec(ixO^S, j)*euler_4d_1vecctor(ixO^S, i))
            enddo
        enddo

    end subroutine reconstruct_4d_2vector_shear_tensor


    subroutine reconstruct_4d_2vector_vorticity_tensor(ixI^L, ixO^L, w, x, vorticityij, vorticitymunu)
        ! Calculate vorticity^munu from vorticity^ij
        use mod_metric, only: lower3_dysp
        use mod_imhd_intermediate, only: dysp_get_lfac
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: vorticityij(ixI^S, 1:3, 1:3)
        double precision, intent(inout)             :: vorticitymunu(ixI^S, 0:3, 0:3)
        ! local variables
        double precision                            :: euler_4d_1vecctor(ixI^S, 0:3)
        double precision                            :: vorticity_proj_3vec(ixI^S, 0:3)
        double precision                            :: lfac_temp(ixI^S), WU_3d_1form(ixI^S, 1:3) ! WU_i
        integer                                     :: i, j

        vorticitymunu = 0.0d0
        vorticity_proj_3vec = 0.0d0
        vorticitymunu(ixO^S, 1:3, 1:3) = vorticityij(ixO^S, 1:3, 1:3)
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, u1_:u3_), WU_3d_1form(ixI^S, 1:3))
        call dysp_get_lfac(ixI^L, ixO^L, w, x, lfac_temp(ixI^S))
        euler_4d_1vecctor(ixO^S, 0) = 1./w(ixO^S, alp_metric_)
        do i=1, ndir
            euler_4d_1vecctor(ixO^S, i) = -w(ixO^S, beta_metric1_+i-1)/w(ixO^S, alp_metric_)
            do j=1, ndir
                vorticity_proj_3vec(ixO^S, i) = vorticity_proj_3vec(ixO^S, i)+&
                    vorticityij(ixO^S, i, j)*WU_3d_1form(ixO^S, j)/lfac_temp(ixO^S)
            enddo
        enddo
        ! calculate vorticity munu
        do i=0, ndir
            do j=0, ndir
                vorticitymunu(ixO^S, i, j) = vorticitymunu(ixO^S, i, j)+(vorticity_proj_3vec(ixO^S, i)* &
                    euler_4d_1vecctor(ixO^S, j)-vorticity_proj_3vec(ixO^S, j)*euler_4d_1vecctor(ixO^S, i))
            enddo
        enddo
    end subroutine reconstruct_4d_2vector_vorticity_tensor


    subroutine check_vort_shear_correct(ixI^L, ixO^L, w, x, Kij, aux_A, aux_lambda, aux_Lambda_scalar, spacial_expansion, check_out)
        ! Check the constraint equations: vortmunu*u_nu, vort_trace, shear^ munu*u_nu, shear_trace, all should be zero!
        use mod_metric
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        include 'amrvacdef.f'

        integer, intent(in)                         :: ixI^L, ixO^L
        double precision, intent(in)                :: w(ixI^S, 1:nw)
        double precision, intent(in)                :: x(ixI^S, 1:ndim)
        double precision, intent(in)                :: Kij(ixI^S, 1:3, 1:3)
        double precision, intent(in)                :: aux_A(ixI^S, ndir), aux_lambda(ixI^S, ndir)
        double precision, intent(in)                :: spacial_expansion(ixI^S), aux_Lambda_scalar(ixI^S)
        double precision, intent(inout)             :: check_out(ixI^S, 1:10)
        ! local variables
        double precision                            :: vorticityij(ixI^S, 1:3, 1:3), vorticitymunu(ixI^S, 0:3, 0:3)
        double precision                            :: shearij(ixI^S, 1:3, 1:3), shearmunu(ixI^S, 0:3, 0:3)
        double precision                            :: u4_vec(ixI^S, 0:3), u4_form(ixI^S, 0:3)
        double precision                            :: gmunu(ixI^S, 0:3, 0:3), gmunu_cp(ixI^S, 0:3, 0:3)
        double precision                            :: beta_3d_1form(ixI^S, 1:3), beta_contract(ixI^S)
        integer                                     :: i, j

        check_out = 0.0d0
        gmunu = 0.0d0
        ! calculate pure spacial projection of shear and vorticity tensor
        do i=1, 3
            do j=1, 3
                call vorticity_spacial_projection_tensor(ixI^L, ixO^L, w, x, Kij, aux_A, aux_lambda, &
                    i, j, vorticityij(ixI^S, i, j))
                call shear_tensor_spacial_projection(ixI^L, ixO^L, w, x, Kij, aux_A, aux_lambda, &
                    aux_Lambda_scalar, spacial_expansion, i, j, shearij(ixI^S, i, j))
            enddo
        enddo
        ! get the 4dim 1vector of four velocity
        call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, u4=u4_vec, gamma=gmunu(ixI^S, 1:3, 1:3))
        call lower4_dysp(ixI^L, ixO^L, w, x, u4_vec, u4_form)
        ! calcualte the full 4dim shear and vorticity tensor
        call reconstruct_4d_2vector_vorticity_tensor(ixI^L, ixO^L, w, x, vorticityij, vorticitymunu)
        call reconstruct_4d_2vector_shear_tensor(ixI^L, ixO^L, w, x, shearij, shearmunu)
        ! get g_munu
        call lower3_dysp(ixI^L, ixO^L, w, x, w(ixI^S, beta_metric1_:beta_metric3_), beta_3d_1form(ixI^S, 1:3))
        beta_contract(ixO^S) = beta_3d_1form(ixO^S, 1)*w(ixO^S, beta_metric1_)+beta_3d_1form(ixO^S, 2)* &
            w(ixO^S, beta_metric2_)+beta_3d_1form(ixO^S, 3)*w(ixO^S, beta_metric3_)
        gmunu(ixO^S, 0, 0) = -w(ixO^S, alp_metric_)**2+beta_contract(ixO^S)
        do i=1, 3
            do j=1, 3
                gmunu(ixO^S, i, j) = gmunu(ixO^S, i, j)*w(ixO^S, psi_metric_)**4
            enddo
            gmunu(ixO^S, 0, i) = beta_3d_1form(ixO^S, i)
            gmunu(ixO^S, i, 0) = beta_3d_1form(ixO^S, i)
        end do
        do i=0, 3
            do j=0, 3
                ! shear should be perpendicular with 4 vector
                check_out(ixO^S, 1+i) = check_out(ixO^S, 1+i)+vorticitymunu(ixO^S, i, j)*u4_form(ixO^S, j)
                ! shear is traceless
                check_out(ixO^S, 5) = check_out(ixO^S, 5)+vorticitymunu(ixO^S, i, j)*gmunu(ixO^S, i, j)
                ! vort should be perpendicular with 4 vector
                check_out(ixO^S, 6+i) = check_out(ixO^S, 6+i)+shearmunu(ixO^S, i, j)*u4_form(ixO^S, j)
                ! shear is traceless
                check_out(ixO^S, 10) = check_out(ixO^S, 10)+shearmunu(ixO^S, i, j)*gmunu(ixO^S, i, j)
            enddo
        enddo
    end subroutine check_vort_shear_correct

end module mod_post_process

