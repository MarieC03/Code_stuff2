!> Module for eos
module mod_eos_tabulated
  use mod_eos
  use mod_eos_tabulated_parameters
  use mod_eos_readtable_scollapse
  use mod_eos_readtable_compose
  public

contains

  !> Read this module's parameters from a file
  subroutine eos_tabulated_read_params()
    include 'amrvacdef.f'


   
    
    select case (table_type_input)
    case ('scollapse')
      table_type = scollapse
      call activate_tablereader_scollapse
    case ('compose')
      table_type = compose
      call activate_tablereader_compose
    case default
      call mpistop('No this table type for tabeos, pls specific one!')
    end select
  end subroutine eos_tabulated_read_params

  subroutine eos_tabulated_activate()
     use mod_eos
     include 'amrvacdef.f'

    eos_type = tabulated

    eos_get_pressure_one_grid         => tabulated_get_pressure_one_grid
    eos_get_eps_one_grid              => tabulated_get_eps_one_grid
    eos_get_cs2_one_grid              => tabulated_get_cs2_one_grid

    eos_get_all_beta_eqm_one_grid     => tabulated_get_all_beta_eqm_one_grid
    eos_get_temp_one_grid             => tabulated_get_temp_one_grid
    eos_eps_get_all_one_grid          => tabulated_eps_get_all_one_grid
    eos_temp_get_all_one_grid         => tabulated_temp_get_all_one_grid

    eos_get_eps_range                 => tabulated_get_eps_range

    if (use_realistic_mp_table) then
       nvars = 25
       if (mype == 0 ) then 
          print*, 'Using realistic mp table for neutrino microphysics'
       endif
    else
       nvars = 6  ! p, eps, cs2, mue, mup, mun
       if (mype == 0 ) then 
          print*, 'Not using realistic mp table for saving memory'
       endif
    endif

    call eos_tabulated_read_params()
    if (mype == 0 ) then 
        print*, '****************************************'
        print*, '-------Tabulated EOS activated----------'
        print*, '****************************************'
    endif
  end subroutine eos_tabulated_activate

  !> Get Pressure based on rho, eps, ye 
  subroutine tabulated_get_pressure_one_grid(prs,rho,eps,temp,ye)
    use mod_eos_interpolation
    
    double precision, intent(inout) :: prs
    double precision, intent(in)    :: rho
    double precision, intent(in)    :: eps
    double precision, intent(inout), optional :: temp
    double precision, intent(in), optional    :: ye

    double precision                :: log_rho, log_temp, log_eps, temp_in, eps_max, eps_min
    double precision                :: pnu(1:3)

    if (rho<=small_rho_thr) then
       call atmo_get_pressure_one_grid(prs,rho,eps)
       if (present(temp)) temp = small_temp
       return
    end if


    if (.not.present(ye))   call mpistop("nuc_eos get_press: need input ye!")

    if (rho>eos_rhomax)     call mpistop("nuc_eos get_press: rho > rhomax")
    if (rho<eos_rhomin)     call mpistop("nuc_eos get_press: rho < rhomin")
    if (ye>eos_yemax)       call mpistop("nuc_eos get_press: ye > yemax")
    if (ye<eos_yemin)       call mpistop("nuc_eos get_press: ye < yemin")


    if (present(temp)) then
      temp_in = temp
      if (temp>eos_tempmax)     call mpistop("nuc_eos get_press: temp > tempmax")
      if (temp<eos_tempmin)     call mpistop("nuc_eos get_press: temp < tempmin")
      log_temp = dlog10(temp_in)
    else
      if (eps>eos_epsmax)     call mpistop("nuc_eos get_press: eps > epsmax")
      if (eps<eos_epsmin)     then
            write(*,*) eps, 'eps in get_press'
            call mpistop("nuc_eos get_press: eps < epsmin")
      endif
    endif

    log_rho  = dlog10(rho)
    log_eps  = dlog10(eps + energy_shift)

    if (.not. present(temp)) then
      if (fil_eps_range_bound) then
        call eos_get_eps_range(rho, eps_min, eps_max, ye)
        if (eps .le. eps_min .or. eps .ge. eps_max) then
           if (eps .le. eps_min) then
             log_temp = dlog10(eos_tempmin)
           endif
           if (eps .ge. eps_max) then
             log_temp = dlog10(eos_tempmax)
           endif
        else
          call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
        endif
      else
        call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
      endif
    endif

    call intep3d(log_rho, log_temp, ye, &
           prs, eos_tables(:,:,:,i_logpress), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

    prs = 10.0d0**prs
  end subroutine tabulated_get_pressure_one_grid

  !> Get specific internal energy from (rho, temp, ye)
  subroutine tabulated_get_eps_one_grid(prs,rho,eps,temp,ye)
    !use mod_eos
    use mod_eos_interpolation

    
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in), optional :: temp, ye
    double precision, intent(inout) :: eps

    double precision             :: log_rho, log_temp

    if (.not.present(temp)) call mpistop("nuc_eos get_eps: need input temp!")
    if (.not.present(ye))   call mpistop("nuc_eos get_eps: need input ye!")

    if (rho<=small_rho_thr) then
       call atmo_get_eps_one_grid(prs,rho,eps)
       return
    end if

    log_rho  = dlog10(rho)
    log_temp = dlog10(temp)

    if (rho>eos_rhomax) then
       write(*,*) rho, 'rho in get_eps'
       call mpistop("nuc_eos get_eps: rho > rhomax")
    endif
    if (rho<eos_rhomin)     call mpistop("nuc_eos get_eps: rho < rhomin")
    if (temp>eos_tempmax)   call mpistop("nuc_eos get_eps: temp > tempmax")
    if (temp<eos_tempmin)   call mpistop("nuc_eos get_eps: temp < tempmin")
    if (ye>eos_yemax)       call mpistop("nuc_eos get_eps: ye > yemax")
    if (ye<eos_yemin)       call mpistop("nuc_eos get_eps: ye < yemin")

    call intep3d(log_rho, log_temp, ye, &
           eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

    eps = 10.0d0**eps - energy_shift

  end subroutine tabulated_get_eps_one_grid

  !> Get cs2 from (rho, eps, ye)
  subroutine tabulated_get_cs2_one_grid(cs2,rho,eps,temp,ye)
    !use mod_eos
    use mod_eos_interpolation
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: ye
    double precision, intent(inout), optional :: temp
    double precision             :: log_rho, log_eps, log_temp, temp_in, eps_max, eps_min

    if (.not.present(ye))   call mpistop("nuc_eos get_cs2: need input ye!")

    if (rho<=small_rho_thr) then
       call atmo_get_cs2_one_grid(cs2,rho,eps)
       if (present(temp)) temp = small_temp
       return
    end if

    if (present(temp)) then
      temp_in = temp
      log_temp = dlog10(temp_in)
    else
      if (eps>eos_epsmax)     call mpistop("nuc_eos get_cs2: eps > epsmax")
      if (eps<eos_epsmin)     call mpistop("nuc_eos get_cs2: eps < epsmin")
    endif

    log_rho  = dlog10(rho)

    if (rho>eos_rhomax)     call mpistop("nuc_eos get_cs2: rho > rhomax")
    if (rho<eos_rhomin)     call mpistop("nuc_eos get_cs2: rho < rhomin")
!    if (temp>eos_tempmax)   call mpistop("nuc_eos get_cs2: temp > tempmax")
!    if (temp<eos_tempmin) log_temp = dlog10(eos_tempmin)
    if (ye>eos_yemax)       call mpistop("nuc_eos get_cs2: ye > yemax")
    if (ye<eos_yemin)       call mpistop("nuc_eos get_cs2: ye < yemin")

    log_eps  = dlog10(eps + energy_shift)

    if (.not. present(temp)) then
      if (fil_eps_range_bound) then
        call eos_get_eps_range(rho, eps_min, eps_max, ye)
        if (eps .le. eps_min .or. eps .ge. eps_max) then
           if (eps .le. eps_min) then
             log_temp = dlog10(eos_tempmin)
           endif
           if (eps .ge. eps_max) then
             log_temp = dlog10(eos_tempmax)
           endif
        else
          call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
        endif
      else
        call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
      endif
    endif

    call intep3d(log_rho, log_temp, ye, &
           cs2, eos_tables(:,:,:,i_cs2), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
    cs2  = max(min(cs2, 1.0d0), 0.0d0)

    cs2  = max(min(cs2, 1.0d0), 0.0d0)

  end subroutine tabulated_get_cs2_one_grid

  !> Get Temperature from (rho, eps, ye)
  subroutine tabulated_get_temp_one_grid(rho,eps,temp,ye)
    !use mod_eos
    use mod_eos_interpolation

    
    double precision, intent(in)    :: rho, ye
    double precision, intent(inout) :: temp, eps

    double precision                :: log_rho, log_temp, log_eps, eps_min, eps_max

    if (rho<=small_rho_thr) then
       temp = max(small_temp, eos_tempmin)
       return
    end if

    if (rho>eos_rhomax)   call mpistop("nuc_eos get_temp: rho > rhomax"  ) 
    if (rho<eos_rhomin)   call mpistop("nuc_eos get_temp: rho < rhomin"  ) 
    if (ye>eos_yemax)     call mpistop("nuc_eos get_temp: ye > yemax"    ) 
    if (ye<eos_yemin)     call mpistop("nuc_eos get_temp: ye < yemin"    ) 
    if (eps>eos_epsmax)   call mpistop("nuc_eos get_temp: eps > epsmax"  ) 
    if (eps<eos_epsmin)   call mpistop("nuc_eos get_temp: eps < epsmin"  ) 

    log_rho = dlog10(rho)
    log_eps = dlog10(eps + energy_shift)

    if (fil_eps_range_bound) then
      call eos_get_eps_range(rho, eps_min, eps_max, ye)
      if (eps .le. eps_min .or. eps .ge. eps_max) then
         if (eps .le. eps_min) then
           log_temp = dlog10(eos_tempmin)
           eps = eps_min
         endif
         if (eps .ge. eps_max) then
           log_temp = dlog10(eos_tempmax)
           eps = eps_max
         endif
      else
        call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
      endif
    else
      call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
    endif

    temp = 10.0d0**log_temp

    ! check the results
    if (temp>eos_tempmax) call mpistop("nuc_eos get_temp: temp > tempmax") 
    if (temp .lt. eos_tempmin) then
        write(*,*) 'rho, eps, temp, ye'
        write(*,*) rho, eps, temp, ye
        call mpistop("nuc_eos get_temp: temp .le. tempmin")
     endif

  end subroutine tabulated_get_temp_one_grid

  !> Get all tab eos var from (rho, eps, ye)
  subroutine tabulated_eps_get_all_one_grid(rho,eps,ye,temp,prs,ent,cs2,dedt,&
                 dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)
    use mod_eos_interpolation

    double precision, intent(inout) :: rho, eps, ye
    double precision, intent(inout), optional :: ent, prs, temp, cs2, dedt
    double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
    double precision, intent(inout), optional :: mu_e,mu_n,mu_p,muhat,munu

    double precision             :: log_rho, log_temp, log_eps, temp_in
    double precision             :: pnu(1:3), epsnu(1:3)
    double precision             :: ffx(nvars), eps_min, eps_max

    if (rho<=small_rho_thr) then
       rho = small_rho
       eps = small_eps
       ye = big_ye
       ! fixme: should include temp and rest of the variables as well
       if (present(prs))    call atmo_get_pressure_one_grid(prs,rho,eps)
       if (present(cs2))    call atmo_get_cs2_one_grid(cs2,rho,eps)
       if (present(ent))    ent  = 4.0d0
       if (present(temp))   temp = small_temp
       return
    end if

    if (rho>eos_rhomax) then
              write(*,*)  rho, 'rho'
            call mpistop("nuc_eos eps_get_all: rho > rhomax")
          
    endif
    if (ye>eos_yemax)   call mpistop("nuc_eos eps_get_all: ye  > yemax ")
    if (ye<eos_yemin)   then 
       write(*,*) 'ye, eos_yemin'
       write(*,*) ye, eos_yemin
       call mpistop("nuc_eos eps_get_all: ye  < yemin ")
    endif
    if (eps>eos_epsmax) call mpistop("nuc_eos eps_get_all: eps > eos_epsmax")
    if (eps<eos_epsmin) then
       write(*,*) 'eps, eos_epsmin'
       write(*,*) eps, eos_epsmin
       call mpistop("nuc_eos eps_get_all: eps < eos_epsmin")
    endif
    
    !if (present(temp)) temp_in = temp

    !log_temp = dlog10(temp_in) 
    log_rho  = dlog10(rho)

    log_eps  = dlog10(eps + energy_shift)

    if (fil_eps_range_bound) then
      call eos_get_eps_range(rho, eps_min, eps_max, ye)
      if (eps .le. eps_min .or. eps .ge. eps_max) then
         if (eps .le. eps_min) then
           eps = eps_min
           log_temp = dlog10(eos_tempmin)
         endif
         if (eps .ge. eps_max) then
           eps = eps_max
           log_temp = dlog10(eos_tempmax)
         endif
      else
        call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
      endif
    else
      call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
      !change eps to consistent with T
      eps = 10.d0**log_eps - energy_shift
    endif

    ! check temp first b4 find the others
    if (present(temp))   temp = 10.d0**log_temp
       if (temp>eos_tempmax) call mpistop("nuc_eos eps_get_all: temp > tempmax") 
       if (temp .lt. eos_tempmin) then
         write(*,*) 'rho, eps, temp, ye,  nuc eos eps get all: temp .lt. tempmin'
         write(*,*) rho, eps, temp, ye
         !temp = max(small_temp,eos_tempmin)
         !log_temp = max(log10(small_temp),log10(eos_tempmin))
         call mpistop("nuc_eos eps_get_all: temp .lt. tempmin") 
       endif

    call intep3d_many(log_rho, log_temp, ye, &
                      ffx, eos_tables,  &
                      nrho, ntemp, nye, nvars,     &
                      logrho_table, logtemp_table, ye_table)

    if (present(prs))  then
          prs  = 10.d0**ffx(i_logpress)
    endif

    if (present(cs2)) then 
       cs2  = max(min(ffx(i_cs2), 1.0d0), 0.0d0)
      ! if (ffx(i_cs2)> 1.00000001d0) then 
      !     write(*,*) 'rho, eps, ye in nuc_eos eps_get_all:'
      !     write(*,*) rho, eps, ye  
      !     write(*,*) cs2, 'cs2'
      !     call mpistop("nuc_eos eps_get_all: cs2>1.0d0")
      ! endif
      ! if (ffx(i_cs2)< 0.0d0) then 
      !     write(*,*) 'rho, eps, ye in nuc_eos eps_get_all:'
      !     write(*,*) rho, eps, ye  
      !     write(*,*) cs2, 'cs2'
      !     call mpistop("nuc_eos eps_get_all: cs2<0.0d0")
      ! endif
    endif
    if (present(mu_e))     mu_e    = ffx(i_mu_e)
    if (present(mu_p))     mu_p    = ffx(i_mu_p)
    if (present(mu_n))     mu_n    = ffx(i_mu_n)

    if (use_realistic_mp_table) then
      if (present(ent))      ent  = ffx(i_entropy)
      if (present(munu))     munu = ffx(i_munu)
      !  derivatives
      if (present(dedt))     dedt    = ffx(i_dedT)
      if (present(dpdrhoe))  dpdrhoe = ffx(i_dpdrhoe)
      if (present(dpderho))  dpderho = ffx(i_dpderho)
      !  chemical potentials
      if (present(muhat))    muhat   = ffx(i_muhat)
      !  compositions
      if (present(xa))       xa      = ffx(i_xa)
      if (present(xh))       xh      = ffx(i_xh)
      if (present(xn))       xn      = ffx(i_xn)
      if (present(xp))       xp      = ffx(i_xp)
      if (present(abar))     abar    = ffx(i_abar)
      if (present(zbar))     zbar    = ffx(i_zbar)
    endif 
  end subroutine tabulated_eps_get_all_one_grid

  !> Get all eos tab var from (rho, temp, ye)
  subroutine tabulated_temp_get_all_one_grid(rho,temp,ye,eps,prs,ent,cs2,dedt,&
                       dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)
    use mod_eos_interpolation

    double precision, intent(in)    :: ye
    double precision, intent(inout) :: eps, rho, temp
    double precision, intent(inout), optional :: ent, prs, cs2, dedt
    double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
    double precision, intent(inout), optional :: mu_e,mu_n,mu_p,muhat,munu

    double precision             :: log_rho, log_temp, log_eps
    double precision             :: ffx(nvars)
    double precision             :: pnu(1:3)

    ! why not modify ye?
    if (rho<=small_rho_thr) then
       rho = small_rho
       eps = small_eps
       temp= small_temp
       ! fixme: should include rest of the variables as well
       if (present(prs))    call atmo_get_pressure_one_grid(prs,rho,eps)
       if (present(cs2))    call atmo_get_cs2_one_grid(cs2,rho,eps)
       if (present(ent)) ent = 4.0d0
       return
    end if

    log_temp = dlog10(temp)   
    log_rho  = dlog10(rho)

    if (rho>eos_rhomax) call mpistop("nuc_eos temp_get_all: rho > rhomax")
    if (ye>eos_yemax)   call mpistop("nuc_eos temp_get_all: ye  > yemax ")
    if (ye<eos_yemin)   call mpistop("nuc_eos temp_get_all: ye  < yemin ")
    if (temp>eos_tempmax)   call mpistop("nuc_eos temp_get_all: temp > tempmax")
    if (temp<eos_tempmin)   call mpistop("nuc_eos temp_get_all: temp < tempmin")
    !if (temp<eos_tempmin) log_temp = dlog10(eos_tempmin)


    call intep3d_many(log_rho, log_temp, ye, &
                      ffx, eos_tables,  &
                      nrho, ntemp, nye, nvars,     &
                      logrho_table, logtemp_table, ye_table)

    eps  = 10.0d0**ffx(i_logenergy) - energy_shift
    eps  = min(eos_epsmax, max(eps, eos_epsmin))
    if (present(prs))  then
       prs  = 10.d0**ffx(i_logpress)
    endif

    if (present(cs2))  then 
        cs2  = max(min(ffx(i_cs2), 1.0d0), 0.0d0)
        !if (ffx(i_cs2)> 1.0d0) stop "nuc_eos in temp get all:cs2>1.0d0" !stop "nuc_eos: cs2 > 1.0d0"
        !if (ffx(i_cs2)< 0.0d0) then
        !    write(*,*) 'rho, eps, ye in temp get all'
        !    write(*,*) rho, eps, ye  
        !    write(*,*) cs2, 'cs2'
        !    stop "nuc_eos in temp get all :cs2<0.0d0" !stop "nuc_eos: cs2 < 0.0d0"
        !endif
    endif
    if (present(mu_e))     mu_e    = ffx(i_mu_e)
    if (present(mu_p))     mu_p    = ffx(i_mu_p)
    if (present(mu_n))     mu_n    = ffx(i_mu_n)

    if (use_realistic_mp_table) then
      if (present(ent))      ent  = ffx(i_entropy)
      if (present(munu))     munu = ffx(i_munu)
      !  derivatives
      if (present(dedt))     dedt    = ffx(i_dedT)
      if (present(dpdrhoe))  dpdrhoe = ffx(i_dpdrhoe)
      if (present(dpderho))  dpderho = ffx(i_dpderho)
      !  chemical potentials
      if (present(muhat))    muhat   = ffx(i_muhat)
      !  compositions
      if (present(xa))       xa      = ffx(i_xa)
      if (present(xh))       xh      = ffx(i_xh)
      if (present(xn))       xn      = ffx(i_xn)
      if (present(xp))       xp      = ffx(i_xp)
      if (present(abar))     abar    = ffx(i_abar)
      if (present(zbar))     zbar    = ffx(i_zbar)
    endif
  end subroutine tabulated_temp_get_all_one_grid

  subroutine tabulated_get_eps_range(rho, eps_min, eps_max, ye)
    use mod_eos_interpolation

    double precision, intent(in) :: rho
    double precision, intent(in), optional :: ye
    double precision, intent(out) :: eps_max, eps_min
    double precision              :: log_rho, log_temp, eps

    if (.not.present(ye))   call mpistop("get_eps_range: Tabulate EOS need ye!")

    !find eps_min    
    log_rho  = dlog10(rho)
    log_temp = dlog10(eos_tempmin)
    
!code test 
!EOS_call_iteration = EOS_call_iteration + 1
    call intep3d(log_rho, log_temp, ye, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
    eps_min = 10.d0**eps - energy_shift

    !find eps_max
    log_rho  = dlog10(rho)
    log_temp = dlog10(eos_tempmax)
    
!code test 
!EOS_call_iteration = EOS_call_iteration + 1
    call intep3d(log_rho, log_temp, ye, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
    eps_max = 10.d0**eps - energy_shift

    ! Debug
    !eps_min = eos_epsmin
    !eps_max = eos_epsmax

  end subroutine tabulated_get_eps_range  

!  subroutine tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
!    use mod_eos_interpolation
!    use mod_rootfinding
!    include 'amrvacdef.f'
!
!    
!    double precision, intent(in)    :: log_rho, log_eps, ye
!    double precision, intent(inout) :: log_temp
!    double precision, save          :: log_rho0, ye0, log_eps_local
!
!    double precision                :: logtemp_min, logtemp_max
!    double precision                :: log_eps_min, log_eps_max
!    ! for root finding
!    integer                         :: error_code
!
!    log_rho0 = log_rho
!    ye0 = ye
!
!    error_code = -1
!    log_eps_local = log_eps
!
!    ! range of the root
!    logtemp_min = logtemp_table(1)
!    logtemp_max = logtemp_table(ntemp)
!    !log_temp = max(logtemp_min, min(logtemp_max, log_temp))
!
!    ! check if the initial guess is close enough
!    !if (dabs(func_eps_of_temp(log_temp)) < eos_precision * dabs(log_eps_local) ) return
!
!    ! get logtemp from logeps
!    !call rootfinding_illinois(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
!    call rootfinding_brent(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
!
!    select case (error_code)
!    !case (0) ! root is found
!    !   write(*,*) "z= ", z
!    case (-1) ! nothing happened
!       call mpistop("have you ever attemp to find the root in tabulated eos?")
!    case (2) ! z is NaN
!       call mpistop("NaN")
!    case (1,3) ! if log_temp cannot be found or is not bracketed
!       ! adjust the range of log_eps, find the root again
!       !write(*,*) 'code test bounded eps'
!       error_code = -1
!       call eos_get_eps_range(10.0d0**log_rho, log_eps_min, log_eps_max, ye)
!       log_eps_max = dlog10(log_eps_max+energy_shift)
!       log_eps_min = dlog10(log_eps_min+energy_shift)
!       log_eps_local = max(log_eps_min, min(log_eps_max, log_eps_local))
!       ! using the bound values of logeps to find logtemp again
!       !call rootfinding_illinois(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
!       call rootfinding_brent(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
!       select case (error_code)
!       !case (0) ! root is found
!       !   write(*,*) "z= ", z
!       case (-1) ! nothing happened
!          call mpistop("forgot to find the root in tabulated eos after adjusting eps range?")
!       case (1) ! failed to find the root
!          call mpistop("Fail to find the root in tabulated eos after adjusting eps range")
!       case (2) ! z is NaN
!          call mpistop("NaN after adjusting eps range")
!       case (3) ! z is not bracketed
!          !write(*,*)  "log_eps, logeps_min, logeps_max, rho, ye "
!          !write(*,*)  log_eps_local, log_eps_min, log_eps_max, 10**log_rho, ye
!          !write(*,*) log_temp, logtemp_min, logtemp_max, func_eps_of_temp(logtemp_min), func_eps_of_temp(logtemp_max)
!          !write(*,*)  "---------------------------"
!          ! If still not bracketed, probably because they are out of the range
!          ! fixme: maybe there are better ways to handle this case
!          !if (log_eps_local .le. log_eps_min) then
!          !   log_temp = logtemp_min
!          !else if (log_eps_local >= log_eps_max) then
!          !   log_temp = logtemp_max
!          !else
!
!               !write(*,*)  'eps,   rho    ye'
!               !write(*,*) 10.d0**log_eps_local - energy_shift, 10.0d0**log_rho, ye
!             if ((10.d0**log_eps_local-energy_shift) .lt. 0.0d0 .and. ye .gt. 0.3d0) then
!                 log_temp = max(log10(eos_tempmin), log10(small_temp))
!                 if (mype == 0) then
!                   !write(*,*) log_temp, 'after artificial treatment inside rootfinding of logeps to logtemp'
!                 endif
!                 return
!             else
!                 ! the table has problem in high density pts, e.g. LS220(stellarcollapse)
!                 log_temp = max(log10(eos_tempmin), log10(small_temp))
!                 return 
!                 ! write(*,*) 'the root is not bracketed in tabulated eos after adjusting eps range &
!                 !          and also hot fixing of temp = tempmin with loweps,highye'
!             !call mpistop("the root is not bracketed in tabulated eos after adjusting eps range &
!             !              and also hot fixing of temp = tempmin with loweps,highye")
!             endif
!          !end if
!          !do error_code = 1,ntemp
!          !   write(*,*) logtemp_table(error_code), func_eps_of_temp(logtemp_table(error_code))
!          !end do
!          !call mpistop("the root is not bracketed in tabulated eos after adjusting eps range")
!       end select
!    end select
!
!    contains
!       double precision function func_eps_of_temp(log_temp_in)
!         double precision, intent(in) :: log_temp_in
!!code test 
!!EOS_call_iteration = EOS_call_iteration + 1
!         call intep3d(log_rho0, log_temp_in, ye0, & 
!                  func_eps_of_temp, eos_tables(:,:,:,i_logenergy), &
!                   nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
!         func_eps_of_temp = 1.0d0 - func_eps_of_temp / log_eps_local 
!       end function
!  end subroutine tabulated_logeps_to_logtemp

  ! fil version
  subroutine tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
    use mod_eos_interpolation
    use mod_rootfinding
    include 'amrvacdef.f'

    
    double precision, intent(in)    :: log_rho, ye
    double precision, intent(inout) :: log_temp, log_eps
    double precision                :: log_rho0, ye0, log_eps_local

    double precision                :: logtemp_min, logtemp_max
    double precision                :: eps, eps_max, eps_min
    ! for root finding
    integer                         :: error_code

    log_rho0 = log_rho
    ye0 = ye

    error_code = -1
    log_eps_local = log_eps

    ! range of the root
    logtemp_min = logtemp_table(1)
    logtemp_max = logtemp_table(ntemp)

    ! get logtemp from logeps
    !call rootfinding_illinois_old(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
    !call rootfinding_illinois(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
    call rootfinding_brent(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
    ! FIL one
    !call rootfinding_zero_brent(log_temp, logtemp_min, logtemp_max, 1.0d-14, func_eps_of_temp)


    select case (error_code)
    !case (0) ! root is found
    !   write(*,*) "z= ", z
    case (-1) ! nothing happened
       call mpistop("have you ever attemp to find the root in tabulated eos?")
    case (2) ! z is NaN
       call mpistop("NaN")
    case (1,3) ! if log_temp cannot be found or is not bracketed
       call eos_get_eps_range(10.0d0**log_rho0, eps_min, eps_max, ye0)
       eps = 10.0d0**(log_eps_local) - energy_shift
       if (eps .le. eps_min .or. eps .ge. eps_max) then
          if (eps .le. eps_min) then
            log_temp = dlog10(eos_tempmin)
            ! return consistent eps as well
            log_eps = dlog10(eps_min + energy_shift)
          endif
          if (eps .ge. eps_max) then
            log_temp = dlog10(eos_tempmax)
            ! return consistent eps as well
            log_eps = dlog10(eps_max + energy_shift)
          endif
       endif
    end select

    contains
       double precision function func_eps_of_temp(log_temp_in)
         double precision, intent(in) :: log_temp_in
         double precision             :: logeps_hat
!code test 
!EOS_call_iteration = EOS_call_iteration + 1
         call intep3d(log_rho0, log_temp_in, ye0, & 
                  logeps_hat, eos_tables(:,:,:,i_logenergy), &
                   nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

         !func_eps_of_temp = 1.0d0 - func_eps_of_temp / log_eps_local 
         !func_eps_of_temp = 1.0d0 - logeps_hat/log_eps_local
         func_eps_of_temp = log_eps_local - logeps_hat

       end function
  end subroutine tabulated_logeps_to_logtemp


  !subroutine tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
  subroutine tabulated_logeps_to_logtemp_old(log_rho,log_eps,log_temp,ye)
      use mod_eos_interpolation
      use mod_rootfinding

      double precision, intent(in) :: log_rho, log_eps, ye
      double precision, intent(inout) :: log_temp
      integer                      :: irho, itemp, iye
      double precision             :: tmp
      double precision             :: logtemp_min, logtemp_max
      ! for root finding
      integer                      :: iter_max = 40  !ntemp
      integer                      :: error_code = -1
      double precision             :: fp, fm, f
      double precision             :: rho_local, eps_min, eps_max, ye_local, eps_local, log_eps_local
      !integer                      :: it_root
      !integer                      :: side=0
      !double precision             :: zp, zm, z

!      call intep3d(log_rho, log_temp, ye, &
!                  tmp, eos_tables(:,:,:,i_logenergy), &
!             nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, irho, itemp, iye)
!
!      f  = tmp - log_eps
!      if (abs(f) < eos_precision*abs(log_eps) ) return


      log_eps_local = log_eps
      rho_local = 10**log_rho
      ye_local  = ye

      logtemp_min = logtemp_table(1)
      logtemp_max = logtemp_table(ntemp)

!      fm = func_eps_of_temp(logtemp_min)
      ! log_eps < eos_tab(log_rho, logtemp_min, ye)
!      if ( fm.ge. 0.d0 ) then
!          log_temp = logtemp_min
!          return
!      endif
!      fp = func_eps_of_temp(logtemp_max)

!      if ( fm*fp .ge. 0.d0 ) call mpistop("the root is not bracketed in tabulated eos (fm*fp > 0)")
!      if ( f *fm .ge. 0.d0 ) logtemp_min = log_temp
!      if ( f *fp .ge. 0.d0 ) logtemp_max = log_temp
!         call eos_get_eps_range(rho_local,eps_min,eps_max,ye_local)
!         if      (10**log_eps .gt. 10**eps_max+energy_shift) then
!              write(*,*) "10**log_eps > 10**eps_max+energy_shift)"
!              log_eps_local = log10(eps_max+energy_shift)
!              write(*,*)  "log_eps, log_eps_local"
!              write(*,*)  log_eps, log_eps_local
!         else if (10**log_eps .lt. 10**eps_min+energy_shift) then
!              write(*,*) "10**log_eps < 10**eps_min+energy_shift)"
!              log_eps_local = log10(eps_min+energy_shift)
!              write(*,*)  "log_eps, log_eps_local"
!              write(*,*)  log_eps, log_eps_local
!         endif

      ! get eps from temperature
      call rootfinding_illinois_old(log_temp, logtemp_min, logtemp_max, eos_precision, iter_max, error_code, func_eps_of_temp)

!      does not bracketed
      if (error_code == 3 .or. error_code == 1) then
         error_code = -1
         call eos_get_eps_range(rho_local,eps_min,eps_max,ye_local)
         if      (10**log_eps .gt. eps_max+energy_shift) then
       !  if      (10**log_eps .gt. 10**eps_max+energy_shift) then
          !    write(*,*) "10**log_eps > eps_max+energy_shift)"
              log_eps_local = log10(eps_max+energy_shift)
          !    write(*,*)  "log_eps, log_eps_local"
          !    write(*,*)  log_eps, log_eps_local
         else if (10**log_eps .lt. eps_min+energy_shift) then
       !  else if (10**log_eps .lt. 10**eps_min+energy_shift) then
          !    write(*,*) "10**log_eps < eps_min+energy_shift)"
              log_eps_local = log10(eps_min+energy_shift)
          !    write(*,*)  "log_eps, log_eps_local"
          !    write(*,*)  log_eps, log_eps_local
         endif
      call rootfinding_illinois_old(log_temp, logtemp_min, logtemp_max, eos_precision, iter_max, error_code, func_eps_of_temp)
      endif

      select case (error_code)
      !case (0) ! root is found
      !   write(*,*) "z= ", z
      case (-1) ! nothing happened
         call mpistop("have you ever attemp to find the root in tabulated eos?")
      case (1) ! failed to find the root
         write(*,*)  "eps, rho, ye "
         write(*,*)  10**log_eps - energy_shift, 10**log_rho, ye
         call mpistop("Fail to find the root in tabulated eos")
      case (2) ! z is NaN
         write(*,*)  "eps, rho, ye "
         write(*,*)  10**log_eps - energy_shift, 10**log_rho, ye
         call mpistop("NaN")
      case (3) ! z is not bracketed
         write(*,*)  "eps, rho, ye "
         write(*,*)  10**log_eps - energy_shift, 10**log_rho, ye
         call mpistop("the root is not bracketed in tabulated eos")
      end select
      contains
         double precision function func_eps_of_temp(log_temp_in)
           double precision, intent(in) :: log_temp_in
           call intep3d(log_rho, log_temp_in, ye, &
                    func_eps_of_temp, eos_tables(:,:,:,i_logenergy), &
                     nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
           func_eps_of_temp = func_eps_of_temp - log_eps_local
         end function
   ! end subroutine tabulated_logeps_to_logtemp
    end subroutine tabulated_logeps_to_logtemp_old
                                                                                                       
  ! Beta equilibrium requires that \mu_n =  \mu_p +\mu_e
  ! Since we assume that (rest mass difference)/T between n and p should be negligible we demand
  ! \mu_n-\mu_p-\mu_e = \mu_hat =0
  subroutine tabulated_get_all_beta_eqm_one_grid(rho,temp,ye,eps,prs)
    use mod_eos_interpolation
    use mod_rootfinding
    !use mod_eos_tabulated_parameters

    double precision, intent(in)              :: rho, temp
    double precision, intent(inout)           :: ye
    double precision, intent(inout), optional :: eps, prs

    double precision, save                    :: log_rho, log_temp
    double precision                          :: log_eps, ye_min, ye_max, log_prs
    integer                                   :: error_code = -1

    log_rho  = dlog10(rho)
    log_temp = dlog10(temp)

    ! given rho and temp, get Ye which satisfies beta equilibrium
    if (rho>eos_rhomax)     call mpistop("nuc_eos beta_eqm: rho > rhomax")
    if (rho<eos_rhomin)     call mpistop("nuc_eos beta_eqm: rho < rhomin")
    if (temp>eos_tempmax)   call mpistop("nuc_eos beta_eqm: temp > tempmax")
    ! if (temp<eos_tempmin)   call mpistop("nuc_eos beta_eqm: temp < tempmin")
    if (temp<eos_tempmin) log_temp = dlog10(eos_tempmin)

    ye_min = ye_table(1)
    ye_max = ye_table(nye)

!code test
!    write(*,*) 'rho, temp'
!    write(*,*) rho, temp
!    write(*,*) ye_min , ye_max, 'ye min max'
!    write(*,*) "################################################"
!    write(6,*) "Not through set_eps_max_min"
!    write(6,*) "Done reading eos tables ,  all in code unit"
!    write(*,*)  eos_yemin, eos_yemax, "eos_yemin  and max"
!    write(*,*)  eos_tempmin, eos_tempmax, "eos_tempmin  and max"
!    write(*,*)  eos_rhomin, eos_rhomax, "eos_rhomin  and max"
!    write(*,*)  eos_epsmin, eos_epsmax, "eos_epsmin  and max"
!    write(*,*)  eos_hmin,            'eos_hmin'
!    write(*,*) "################################################"


    call rootfinding_brent(ye, ye_min, ye_max, eos_precision, eos_iter_max, error_code, func_munu_of_ye)
!  stop 'beta eqm'
    select case (error_code)
    !case (0) ! root is found
    !   write(*,*) "z= ", z
    case (-1) ! nothing happened
       call mpistop("have you ever attemp to find the root in tabulated eos?")
    case (1) ! failed to find the root
       call mpistop("Fail to find the root in beta_eq")
    case (2) ! z is NaN
       call mpistop("NaN in beta_eq")
    case (3) ! z is not bracketed
       call mpistop("the root is not bracketed in beta_eq")
    end select
    ! after got the correct ye, check the EOS range 
    if (ye>eos_yemax)       call mpistop("nuc_eos beta_eqm: ye > yemax")
    if (ye<eos_yemin)       call mpistop("nuc_eos beta_eqm: ye < yemin")

    if (present(eps)) then
           call intep3d(log_rho, log_temp, ye, &
                    log_eps, eos_tables(:,:,:,i_logenergy), &
                     nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
           eps = 10.0**log_eps - energy_shift
    endif

    if (present(prs)) then
           call intep3d(log_rho, log_temp, ye, &
                    log_prs, eos_tables(:,:,:,i_logpress), &
                     nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
           prs = 10.0**log_prs
    endif



    contains
       double precision function func_munu_of_ye(ye_in)
         double precision, intent(in) :: ye_in
         double precision             :: ffx(nvars)
 
         call intep3d_many(log_rho, log_temp, ye_in, &
                      ffx, eos_tables,  &
                      nrho, ntemp, nye, nvars,     &
                      logrho_table, logtemp_table, ye_table)
         !f(Ye) =  mu_e(Ye) + mu_p(Ye) - mu_n(Ye)  ! included the mass difference Qnp (we should)
          !write(*,*) ffx(i_mu_e), ffx(i_mu_p), ffx(i_mu_n) , ye_in
         ! seems adding the later term do no match with the FIL result. FIL also wrong, so we do not add Qnp
         func_munu_of_ye = ffx(i_mu_e) + ffx(i_mu_p) - ffx(i_mu_n) - Qnp   !KEN
       end function
  end subroutine tabulated_get_all_beta_eqm_one_grid

end module mod_eos_tabulated
