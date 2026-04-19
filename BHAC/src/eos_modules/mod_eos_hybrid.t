!> T_eps for hybrid eos is eps
module mod_eos_hybrid


    implicit none
    public

    double precision, allocatable    :: eos_tables_cold_betaeq(:, :)
    integer, parameter               :: hybrid_tabulated = 1, hybrid_lorene=2, hybrid_poly = 3
    integer, parameter               :: eos_table_interp_order = 1
    integer, parameter               :: idx_ye=1, idx_eps=2, idx_prs=3, idx_cs2=4
    double precision                 :: dlog_rho = 0.0d0 ! assume the logrho_table is equidistance in log-scale
    logical                          :: debug_mod_eos_hybrid = .false.
    logical                          :: is_equidistant = .false.

contains


    subroutine eos_hybrid_activate()
        use mod_eos
        include 'amrvacdef.f'
        eos_type = hybrid
        ! check if the parameters are make sences
        if ((eos_gamma_th < 1.0d0) .or. (eos_gamma_th > 2.0d0+10*tiny(1.0d0))) call mpistop ("Error: eos_gamma_th < 1 or eos_gamma_th > 2")
        select case (table_type_input)
        case ("polytrope")
            if (eos_gamma <= 0.0d0) call mpistop ("Error: eos_gamma <= 0")
            if (eos_adiab < 0.0d0) call mpistop  ("Error: eos_adiab < 0")
            eos_get_pressure_one_grid         => hybrid_poly_get_pressure_one_grid
            eos_get_eps_one_grid              => hybrid_poly_get_eps_one_grid
            eos_get_cs2_one_grid              => hybrid_poly_get_cs2_one_grid
            eos_get_temp_one_grid             => hybrid_poly_get_temp_one_grid
            eos_get_eps_range                 => hybrid_poly_get_eps_range
            if (mype == 0 ) then
                write(*, *) 'eos_gamma_th: ', eos_gamma_th
                write(*, *) 'eos_gamma: ', eos_gamma 
                write(*, *) 'eos_adiab: ', eos_adiab
                print*, '****************************************'
                print*, '-------Hybrid_poly EOS activated--------'
                print*, '****************************************'
            endif
        case ("compose", "scollapse")
            eos_get_pressure_one_grid         => hybrid_tab_get_pressure_one_grid
            eos_get_eps_one_grid              => hybrid_tab_get_eps_one_grid
            eos_get_cs2_one_grid              => hybrid_tab_get_cs2_one_grid
            eos_get_temp_one_grid             => hybrid_tab_get_temp_one_grid
            eos_get_all_beta_eqm_one_grid     => hybrid_tab_get_all_beta_eqm_one_grid
            eos_get_eps_range                 => hybrid_tab_get_eps_range
            if (mype == 0 ) then
                write(*, *) 'eos_gamma_th: ', eos_gamma_th
                write(*, *) 'tabulated file: ', eos_table_name
                print*, '****************************************'
                print*, '-----Hybrid_tabulated EOS activated-----'
                print*, '****************************************'
            endif
            call get_zero_temp_betaeqm_sequence()
        case ("Lorene")
            eos_get_pressure_one_grid         => hybrid_tab_get_pressure_one_grid
            eos_get_eps_one_grid              => hybrid_tab_get_eps_one_grid
            eos_get_cs2_one_grid              => hybrid_tab_get_cs2_one_grid
            eos_get_temp_one_grid             => hybrid_tab_get_temp_one_grid
            eos_get_all_beta_eqm_one_grid     => hybrid_tab_get_all_beta_eqm_one_grid
            eos_get_eps_range                 => hybrid_tab_get_eps_range
            if (mype == 0 ) then
                write(*, *) 'eos_gamma_th: ', eos_gamma_th
                write(*, *) 'tabulated file: ', eos_table_name
                print*, '****************************************'
                print*, '-----Hybrid_lorene EOS activated-----'
                print*, '****************************************'
            endif
            call read_lorene_table()
        case default
            call mpistop("Unknown hybrid eos type!")
        end select

  
        call eos_check

    end subroutine eos_hybrid_activate


!==================================================================
!===================== auxiliary functions ========================
!==================================================================


    ! Subroutine to perform binary search and find the bounds
    subroutine binary_search_bounds(arr, n, svalue, lower_idx, upper_idx)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: arr(n), svalue
        integer, intent(out) :: lower_idx, upper_idx
        integer :: low, high, mid
    
        ! Initialize bounds to -1 (indicating no bound)
        low = 1
        high = n
        lower_idx = -1
        upper_idx = -1
        ! Binary search loop
        do while (low <= high)
            mid = (low + high) / 2
            if (arr(mid) == svalue) then
                lower_idx = mid
                upper_idx = mid
                return
            elseif (arr(mid) < svalue) then
                low = mid + 1
            else
                high = mid - 1
            end if
        end do
        ! Set the indices based on the binary search result
        if (high >= 1) lower_idx = high   ! Index of the largest element smaller than svalue
        if (low <= n) upper_idx = low     ! Index of the smallest element larger than svalue
    
    end subroutine binary_search_bounds

    subroutine check_equidistant(arr, eq_dist)
        implicit none
        double precision, intent(in) :: arr(:)
        logical, intent(out) :: eq_dist 
        integer :: i, n
        double precision :: diff, first_diff
        double precision, parameter :: epsilon = 1.0e-10 ! Tolerance for round-off error

        n = size(arr)
        eq_dist = .true.
        ! If array has less than 2 elements, we can't determine equidistance
        if (n < 2) then
            return
        end if
        first_diff = arr(2) - arr(1)
        do i = 2, n-1
            diff = arr(i+1) - arr(i)
            if (abs(diff - first_diff) > epsilon) then
                eq_dist = .false.
                return
            end if
        end do
    end subroutine check_equidistant

    subroutine read_lorene_table()
        use mod_eos
        include 'amrvacdef.f'
        integer :: i, ios, dummy
        character(len=100) :: line
        double precision, allocatable   :: temp_eosh(:)

        ! Open the file in SHARED mode for multithreading
        open(unit=10, file=eos_table_name, status="old", access="sequential", action="read", &
             position="rewind", shared, iostat=ios)
        ! skip the frist five rows in Lorene table
        do i = 1, 5
            read(10, '(A)', iostat=ios) line
        end do
        read(10, *, iostat=ios) nrho
        ! skip another 3 rows in Lorene table
        do i = 1, 3
            read(10, '(A)', iostat=ios) line
        end do
        if (allocated(logrho_table)) deallocate(logrho_table)
        if (allocated(eos_tables_cold_betaeq)) deallocate(eos_tables_cold_betaeq)
        allocate(logrho_table(nrho)) ! not log yet, take log at last step
        allocate(eos_tables_cold_betaeq(1:nrho, 1:4)) ! in the order of ye, eps, prs, cs^2, no log envolved
        allocate(temp_eosh(1:nrho))
        eos_tables_cold_betaeq = 0.0d0  ! all ye should be zero because they are not used
        do i = 1, nrho
            read(10, *, iostat=ios) dummy, logrho_table(i), eos_tables_cold_betaeq(i, idx_eps), eos_tables_cold_betaeq(i, idx_prs)
            if (ios /= 0) then
                call mpistop("Error reading Lorene eos table or end of file encountered")
            end if
        end do
        close(10)
        ! processing the table
        logrho_table = logrho_table*(massn_cgs*cm3_to_fm3*rho_gf) ! rho (not take log yet)
        eos_tables_cold_betaeq(:, idx_eps) = eos_tables_cold_betaeq(:, idx_eps)*rho_gf/logrho_table-1 ! eps
        eos_tables_cold_betaeq(:, idx_prs) = eos_tables_cold_betaeq(:, idx_prs)*press_gf ! pressure
        temp_eosh(1:nrho) = 1+eos_tables_cold_betaeq(1:nrho, idx_eps)+eos_tables_cold_betaeq(1:nrho, idx_prs)/logrho_table(1:nrho) ! specific enthalpy
        ! calculate square of sound speed
        eos_tables_cold_betaeq(1, idx_cs2) = (eos_tables_cold_betaeq(2, idx_prs)-eos_tables_cold_betaeq(1, idx_prs))/&
            (logrho_table(2)-logrho_table(1))/temp_eosh(1) ! first order at head 
        eos_tables_cold_betaeq(nrho, idx_cs2) = (eos_tables_cold_betaeq(nrho, idx_prs)-eos_tables_cold_betaeq(nrho-1, idx_prs))/&
            (logrho_table(nrho)-logrho_table(nrho-1))/temp_eosh(nrho) ! first order at tail 
        do i = 2, nrho-1
            eos_tables_cold_betaeq(i, idx_cs2) = (eos_tables_cold_betaeq(i+1, idx_prs)-eos_tables_cold_betaeq(i-1, idx_prs))/&
                (logrho_table(i+1)-logrho_table(i-1))/temp_eosh(i) ! second order at intermediate
        enddo
        eos_tables_cold_betaeq(:, idx_cs2) = max(0.0d0, min(eos_tables_cold_betaeq(:, idx_cs2), 1.0d0))
        ! apply energy shift if needed
        if (eos_tables_cold_betaeq(1, idx_eps) < 0.0d0) then
            energy_shift = -2*eos_tables_cold_betaeq(1, idx_eps)
            eos_tables_cold_betaeq(:, idx_eps) = eos_tables_cold_betaeq(:, idx_eps)+energy_shift
        endif
        ! init limits of eos
        eos_rhomax = logrho_table(nrho)
        eos_rhomin = logrho_table(1)
        eos_epsmax = eos_tables_cold_betaeq(nrho, idx_eps)
        eos_epsmin = eos_tables_cold_betaeq(1, idx_eps)
        eos_hmin = temp_eosh(1)
        ! take log rho and check equidistant
        logrho_table = dlog10(logrho_table) ! take log at this step
        call check_equidistant(logrho_table, is_equidistant)
        if (is_equidistant) dlog_rho = (logrho_table(nrho)-logrho_table(1))/(nrho-1)
        if (mype==0) then
            if (is_equidistant) then
                write(*, *) "Equi-distant table in log rho!"
            else
                write(*, *) "Non equi-distant table in log rho!"
            endif
        endif
        deallocate(temp_eosh)
    end subroutine


!==================================================================
!===================== hybrid polytrope ===========================
!==================================================================


    subroutine hybrid_poly_get_pressure_one_grid(prs, rho, eps, temp, ye)
  
        use mod_eos
        include 'amrvacdef.f'
        
        double precision, intent(inout) :: prs
        double precision, intent(in) :: rho
        double precision, intent(in) :: eps
        double precision, intent(in), optional :: ye
        double precision, intent(inout), optional :: temp
        ! local variable below
        double precision             :: eps_cold, eps_th
    
        if (.not. old_bhac_safety) then
            if (rho<small_rho_thr) then
                call atmo_get_pressure_one_grid(prs,rho,eps)
                return
            end if
        endif
    
        eps_cold = (eos_adiab*rho**(eos_gamma-1.0d0))/(eos_gamma-1.0d0)
        !here may be different with FIL where they change eps to be eps_cold if eps<eps_cold
        eps_th = max(eps-eps_cold, 0.0d0)
        prs = (eos_gamma_th-1.0d0)*rho*eps_th+eos_adiab*rho**eos_gamma ! the same for the two gammas?
  
    end subroutine hybrid_poly_get_pressure_one_grid


    subroutine hybrid_poly_get_eps_one_grid(prs, rho, eps, temp, ye)
  
        use mod_eos
        include 'amrvacdef.f'
        
        double precision, intent(in) :: prs
        double precision, intent(in) :: rho
        double precision, intent(inout) :: eps
        double precision, intent(in), optional :: temp, ye
        ! local variables
        double precision             :: prs_cold, eps_cold
        double precision             :: prs_th
    
        if (.not. old_bhac_safety) then
            if (rho<small_rho_thr) then
                call atmo_get_eps_one_grid(prs,rho,eps)
                return
            end if
        endif
        prs_cold = eos_adiab*rho**eos_gamma
        eps_cold = (eos_adiab*rho**(eos_gamma-1.0d0))/(eos_gamma-1.0d0)
        !here may be different with FIL where they change prs to be prs_cold if prs<prs_cold
        prs_th = max(prs-prs_cold, 0.0d0)
        if (abs(eos_gamma_th-1.0d0)<2*tiny(1.0d0)) then
            eps = eps_cold
        else
            eps = eps_cold+prs_th/(rho*(eos_gamma_th-1.0d0))
        endif
  
    end subroutine hybrid_poly_get_eps_one_grid


    subroutine hybrid_poly_get_cs2_one_grid(cs2, rho, eps, temp, ye)
        use mod_eos
        include 'amrvacdef.f'
    
        double precision, intent(inout) :: cs2
        double precision, intent(in) :: rho
        double precision, intent(in) :: eps
        double precision, intent(in), optional :: ye
        double precision, intent(inout), optional :: temp
        ! local variables
        double precision             :: eps_cold, eps_th
        double precision             :: prs
        double precision             :: dpdrho_cold
        double precision             :: enthalpy

        if (.not. old_bhac_safety) then
            if (rho<small_rho_thr) then
                call atmo_get_cs2_one_grid(cs2,rho,eps)
                return
            end if
        endif

        eps_cold = (eos_adiab*rho**(eos_gamma-1.0d0))/(eos_gamma-1.0d0)
        !here may be different with FIL where they change eps to be eps_cold if eps<eps_cold
        eps_th = max(eps-eps_cold, 0.0d0)
        prs = (eos_gamma_th-1.0d0)*rho*eps_th+eos_adiab*rho**eos_gamma
        dpdrho_cold = eos_gamma*eos_adiab*rho**(eos_gamma-1.0d0)
        enthalpy = 1.0d0+eps+prs/rho
        !(See Miguel Alcubierre 2007, Eq.7.6.44, Page 262; and Luciano 2013, Eq.2.173, Page 108)
        cs2 = (dpdrho_cold+eos_gamma_th*(eos_gamma_th-1.0d0)*eps_th)/enthalpy
        cs2 = min(max(0.0, cs2), 0.9999d0)
    end subroutine hybrid_poly_get_cs2_one_grid


    subroutine hybrid_poly_get_temp_one_grid(rho, eps, temp, ye)
        use mod_eos
        double precision, intent(in)    :: rho 
        double precision, intent(in)    :: ye
        double precision, intent(inout) :: temp, eps
        ! local
        double precision                :: eps_cold, eps_th

        eps_cold = (eos_adiab*rho**(eos_gamma-1.0d0))/(eos_gamma-1.0d0)
        !here may be different with FIL where they change eps to be eps_cold if eps<eps_cold
        eps_th = max(eps-eps_cold, 0.0d0)
        temp = eps_th*(eos_gamma_th-1.0d0)
    end subroutine hybrid_poly_get_temp_one_grid


    subroutine hybrid_poly_get_eps_range(rho, eps_min, eps_max, ye)
        use mod_eos
        double precision, intent(in) :: rho
        double precision, intent(in), optional :: ye
        double precision, intent(out) :: eps_max, eps_min

        !eps_min = small_eps
        eps_min = (eos_adiab*rho**(eos_gamma-1.0d0))/(eos_gamma-1.0d0)
        eps_max = bigdouble_eos
    end subroutine


!==================================================================
!===================== hybrid tabulated ===========================
!==================================================================


    subroutine search_and_interpolate_cold_eos_table(intp_res, rho, idx_var, log_interp)
        ! cold table in the order of ye, eps, prs, cs^2, no log envolved
        use mod_eos
        use mod_interpolate, only : lagrange_interpolation
        double precision, intent(out)         :: intp_res
        double precision, intent(in)          :: rho
        integer, intent(in)                   :: idx_var
        logical, intent(in)                   :: log_interp
        ! local
        integer                               :: low_idx, high_idx
        double precision                      :: log_rho
        logical                               :: local_log_interp

        local_log_interp = log_interp
        log_rho = dlog10(rho)
        if (log_rho<logrho_table(1) .or. log_rho>logrho_table(nrho)) then 
            write(*, *) "logrho: ", log_rho, ", range: ", logrho_table(1), " to ", logrho_table(nrho)
            call mpistop("logrho outside the allowed interpolation range.")
        endif
        if (is_equidistant) then
            low_idx = INT((log_rho-logrho_table(1))/dlog_rho)+1
            if (low_idx>nrho) call mpistop("Could not find rho in eos table: rho larger than all value in the table")
            do while(logrho_table(low_idx)>log_rho)
                !write(*, *) 'ddebug', logrho_table(1), logrho_table(nrho), dlog_rho, log_rho, low_idx
                low_idx = low_idx-1
                if (low_idx<1) call mpistop("Couldn't find rho in eos table.")
            enddo
        else
            call binary_search_bounds(logrho_table, nrho, log_rho, low_idx, high_idx)
            if (low_idx==-1) call mpistop("Could not find rho in eos table: rho less than all value in the table")
            if (high_idx==-1) call mpistop("Could not find rho in eos table: rho larger than all value in the table")
        endif
        high_idx = low_idx+eos_table_interp_order ! high_idx is changed anyway
        ! debug use
        if (debug_mod_eos_hybrid) then
            write(*, *) "rho_start, rho_end, dlog_rho, rho_aim, find index, interp var index"
            write(*, *) logrho_table(1), logrho_table(nrho), dlog_rho, log_rho, low_idx, idx_var
            write(*, *) "your logrho: ", log_rho, ", find rho index: ", low_idx, ", corresponding logrho: ", logrho_table(low_idx)
            write(*, *) "    neighbor eps: ", eos_tables_cold_betaeq(low_idx, 2), ", neighbor pressure: ", eos_tables_cold_betaeq(low_idx, 3)
            write(*, *) "    neighbor eps: ", eos_tables_cold_betaeq(high_idx, 2), ", neighbor pressure: ", eos_tables_cold_betaeq(high_idx, 3)
        endif
        if (high_idx>nrho) call mpistop("Try to interpolate at rho larger than the largest table value.")
        ! avoid log interp for values <= 0
        if (eos_tables_cold_betaeq(low_idx, idx_var)<=tiny(1.0d0)) local_log_interp=.false.
        if (local_log_interp) then
            call lagrange_interpolation(logrho_table(low_idx:high_idx), &
                dlog10(eos_tables_cold_betaeq(low_idx:high_idx, idx_var)), log_rho, intp_res)
            intp_res = 10.0d0**intp_res
        else
            call lagrange_interpolation(logrho_table(low_idx:high_idx), &
                eos_tables_cold_betaeq(low_idx:high_idx, idx_var), log_rho, intp_res)
        endif
        if ((idx_var==idx_eps) .and. (energy_shift>tiny(1.0d0))) intp_res = intp_res-energy_shift
    end subroutine search_and_interpolate_cold_eos_table 


    subroutine get_zero_temp_betaeqm_sequence()
        use mod_eos
        use mod_eos_tabulated, only: eos_tabulated_read_params, &
            tabulated_get_all_beta_eqm_one_grid, tabulated_get_cs2_one_grid
        include 'amrvacdef.f'
        double precision           :: zero_temp, rho
        integer                    :: i

        call eos_tabulated_read_params()
        zero_temp = 10.0d0**logtemp_table(1)
        allocate(eos_tables_cold_betaeq(1:nrho, 1:4)) ! in the order of ye, eps, prs, cs^2, no log involved
        do i=1, nrho
            rho = 10.0d0**logrho_table(i)
            call tabulated_get_all_beta_eqm_one_grid(rho, zero_temp, &
                eos_tables_cold_betaeq(i, idx_ye), eps=eos_tables_cold_betaeq(i, idx_eps), &
                prs=eos_tables_cold_betaeq(i, idx_prs)) ! return true eps
            call tabulated_get_cs2_one_grid(eos_tables_cold_betaeq(i, idx_cs2), rho, &
                eos_tables_cold_betaeq(i, idx_eps), ye=eos_tables_cold_betaeq(i, idx_ye), temp=zero_temp) ! use true eps
        enddo
        ! energy_shift already assigned by eos reader, we should do the shift again because we get the true eps now which could <0
        eos_tables_cold_betaeq(i, idx_eps) = eos_tables_cold_betaeq(i, idx_eps)+energy_shift
        call check_equidistant(logrho_table, is_equidistant)
        if (is_equidistant) dlog_rho = (logrho_table(nrho)-logrho_table(1))/(nrho-1)
        if (mype==0) then
            if (is_equidistant) then
                write(*, *) "Equi-distant table in log rho!"
            else
                write(*, *) "Non equi-distant table in log rho!"
            endif
        endif
        deallocate(eos_tables)
    end subroutine get_zero_temp_betaeqm_sequence


    subroutine hybrid_tab_get_pressure_one_grid(prs, rho, eps, temp, ye)
  
        use mod_eos
        include 'amrvacdef.f'
        
        double precision, intent(inout) :: prs
        double precision, intent(in) :: rho
        double precision, intent(in) :: eps
        double precision, intent(in), optional :: ye
        double precision, intent(inout), optional :: temp
        ! local variable below
        double precision             :: eps_cold, eps_th, prs_cold
    
        ! check boundary 
        if (.not. old_bhac_safety) then
            if (rho<small_rho_thr) then
                call atmo_get_pressure_one_grid(prs,rho,eps)
                return
            endif
        endif
        if (rho>eos_rhomax .or. rho<eos_rhomin) then
           write(*, '(A4, ES12.4, A56, ES12.4, A2, ES12.4, A1)') 'rho ', rho, &
                ' in hybrid_tab_get_pressure_one_grid outside the bound (', eos_rhomin, ', ', eos_rhomax, ')'
           call mpistop("hybrid_eos: rho outside the bound allowed by the cold beta_eqm table")
        endif
        ! calculate through interpolation
        call search_and_interpolate_cold_eos_table(eps_cold, rho, idx_eps, .true.)
        call search_and_interpolate_cold_eos_table(prs_cold, rho, idx_prs, .true.)
        !if (eps_cold>=eps) write(*, *) "Eps_cold larger than eps! set pressure=pressure_cold."
        !here may be different with FIL where they change prs to be prs_cold if prs<prs_cold
        if (debug_mod_eos_hybrid) then
            write(*, *) "cold eps: ", eps_cold, ", cold pressure: ", prs_cold
        endif
        eps_th = max(eps-eps_cold, 0.0d0)
        prs = (eos_gamma_th-1.0d0)*rho*eps_th+prs_cold
  
    end subroutine hybrid_tab_get_pressure_one_grid


    subroutine hybrid_tab_get_eps_one_grid(prs, rho, eps, temp, ye)
  
        use mod_eos
        include 'amrvacdef.f'
        
        double precision, intent(in) :: prs
        double precision, intent(in) :: rho
        double precision, intent(inout) :: eps
        double precision, intent(in), optional :: temp, ye
        ! local variables
        double precision             :: prs_cold, eps_cold
        double precision             :: prs_th
    
        ! check atmo and boundary
        if (.not. old_bhac_safety) then
            if (rho<small_rho_thr) then
                call atmo_get_eps_one_grid(prs,rho,eps)
                return
            end if
        endif
        if (rho>eos_rhomax .or. rho<eos_rhomin) then
           write(*, '(A4, ES12.4, A51, ES12.4, A2, ES12.4, A1)') 'rho ', rho, &
                ' in hybrid_tab_get_eps_one_grid outside the bound (', eos_rhomin, ', ', eos_rhomax, ')'
           call mpistop("hybrid_eos: rho outside the bound allowed by the cold beta_eqm table")
        endif
        ! calucate through intepolation
        call search_and_interpolate_cold_eos_table(eps_cold, rho, idx_eps, .true.)
        call search_and_interpolate_cold_eos_table(prs_cold, rho, idx_prs, .true.)
        ! here may be different with FIL where they change prs to be prs_cold if prs<prs_cold
        prs_th = max(prs-prs_cold, 0.0d0)
        !if (prs_cold>=prs) write(*, *) "Pressure_cold larger than total pressure! set eps=eps_cold."
        if (abs(eos_gamma_th-1)<2*tiny(1.0d0)) then
            eps = eps_cold
        else
            eps = eps_cold+prs_th/(rho*(eos_gamma_th-1.0d0))
        endif
  
    end subroutine hybrid_tab_get_eps_one_grid


    subroutine hybrid_tab_get_cs2_one_grid(cs2, rho, eps, temp, ye)
        use mod_eos
        include 'amrvacdef.f'
        double precision, intent(inout) :: cs2
        double precision, intent(in) :: rho
        double precision, intent(in) :: eps
        double precision, intent(in), optional :: ye
        double precision, intent(inout), optional :: temp
        ! local variables
        double precision             :: prs, eps_th, enthalpy
        double precision             :: cs2_cold, eps_cold, enthalpy_cold, prs_cold
    
        ! check boundary
        if (.not. old_bhac_safety) then
            if (rho<small_rho_thr) then
                call atmo_get_cs2_one_grid(cs2,rho,eps)
                return
            end if
        endif
        if (rho>eos_rhomax .or. rho<eos_rhomin) then
           write(*, '(A4, ES12.4, A51, ES12.4, A2, ES12.4, A1)') 'rho ', rho, &
                ' in hybrid_tab_get_cs2_one_grid outside the bound (', eos_rhomin, ', ', eos_rhomax, ')'
           call mpistop("hybrid_eos: rho outside the bound allowed by the cold beta_eqm table")
        endif
        ! calculate by interpolate
        call search_and_interpolate_cold_eos_table(eps_cold, rho, idx_eps, .true.)
        call search_and_interpolate_cold_eos_table(prs_cold, rho, idx_prs, .true.)
        call search_and_interpolate_cold_eos_table(cs2_cold, rho, idx_cs2, .true.)
        ! here may be different with FIL where they change eps to be eps_cold if eps<eps_cold
        eps_th = max(eps-eps_cold, 0.0d0)
        prs = (eos_gamma_th-1.0d0)*rho*eps_th+prs_cold
        enthalpy_cold = 1.0d0+eps_cold+prs_cold/rho
        enthalpy = 1.0d0+eps+prs/rho
        !(See Miguel Alcubierre 2007, Eq.7.6.44, Page 262; and Luciano 2013, Eq.2.173, Page 108)
        cs2 = (cs2_cold*enthalpy_cold+eos_gamma_th*(eos_gamma_th-1.0d0)*eps_th)/enthalpy
        cs2 = min(max(0.0, cs2), 0.99999999d0)

    end subroutine hybrid_tab_get_cs2_one_grid


    subroutine hybrid_tab_get_eps_range(rho, eps_min, eps_max, ye)
        use mod_eos
        double precision, intent(in) :: rho
        double precision, intent(in), optional :: ye
        double precision, intent(out) :: eps_max, eps_min

        !eps_min = small_eps
        call search_and_interpolate_cold_eos_table(eps_min, rho, idx_eps, .true.)
        eps_max = bigdouble_eos
    end subroutine


    subroutine hybrid_tab_get_temp_one_grid(rho, eps, temp, ye)
        use mod_eos
        double precision, intent(in) :: rho 
        double precision, intent(in) :: ye
        double precision, intent(inout) :: temp, eps
        ! local
        double precision             :: eps_cold, eps_th

        if (rho>eos_rhomax .or. rho<eos_rhomin) then
           write(*, '(A4, ES12.4, A52, ES12.4, A2, ES12.4, A1)') 'rho ', rho, &
                ' in hybrid_tab_get_temp_one_grid outside the bound (', eos_rhomin, ', ', eos_rhomax, ')'
           call mpistop("hybrid_eos: rho outside the bound allowed by the cold beta_eqm table")
        endif
        call search_and_interpolate_cold_eos_table(eps_cold, rho, idx_eps, .true.)
        ! here may be different with FIL where they change eps to be eps_cold if eps<eps_cold
        eps_th = max(eps-eps_cold, 0.0d0)
        temp = eps_th*(eos_gamma_th-1.)
    end subroutine hybrid_tab_get_temp_one_grid


    subroutine hybrid_tab_get_all_beta_eqm_one_grid(rho, temp, ye, eps, prs)
        use mod_eos
        double precision, intent(in)              :: rho, temp
        double precision, intent(inout)           :: ye
        double precision, intent(inout), optional :: eps, prs
    
        ! given rho and temp, get Ye which satisfies beta equilibrium
        if (rho>eos_rhomax .or. rho<eos_rhomin) then
           write(*, '(A4, ES12.4, A60, ES12.4, A2, ES12.4, A1)') 'rho ', rho, &
                ' in hybrid_tab_get_all_beta_eqm_one_grid outside the bound (', eos_rhomin, ', ', eos_rhomax, ')'
           call mpistop("hybrid_eos: rho outside the bound allowed by the cold beta_eqm table")
        endif
        call search_and_interpolate_cold_eos_table(ye, rho, idx_ye, .false.)
        if (present(eps)) then 
            call search_and_interpolate_cold_eos_table(eps, rho, idx_eps, .true.)
        endif
        if (present(prs)) call search_and_interpolate_cold_eos_table(prs, rho, idx_prs, .true.)
        
    end subroutine hybrid_tab_get_all_beta_eqm_one_grid

end module mod_eos_hybrid
