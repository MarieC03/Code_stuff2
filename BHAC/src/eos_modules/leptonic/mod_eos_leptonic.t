module mod_eos_leptonic
  use mod_eos
  use mod_eos_leptonic_parameters
  use mod_eos_leptonic_interpolation
  use mod_rootfinding
  implicit none
  public

contains

  subroutine eos_leptonic_activate()
    use mod_eos_readtable_leptonic_scollapse
    use mod_eos_readtable_leptonic
    include 'amrvacdef.f'

    select case (trim(baryon_table_type))
    case ('scollapse')
      call activate_baryon_tablereader_scollapse()
    case default
      call mpistop('eos_leptonic_activate: unsupported baryon table type')
    end select

    call activate_tablereader_leptonic()

    eos_type    = leptonic
    eos_rhomin  = lep_eos_rhomin
    eos_rhomax  = lep_eos_rhomax
    eos_tempmin = lep_eos_tempmin
    eos_tempmax = lep_eos_tempmax
    eos_yemin   = lep_eos_ylemin
    eos_yemax   = lep_eos_ylemax
    eos_ymumin  = lep_eos_ymumin
    eos_ymumax  = lep_eos_ymumax
    eos_epsmin  = lep_eos_epsmin
    eos_epsmax  = lep_eos_epsmax
    eos_hmin    = lep_eos_hmin

    eos_get_pressure_one_grid     => leptonic_get_pressure_one_grid
    eos_get_eps_one_grid          => leptonic_get_eps_one_grid
    eos_get_cs2_one_grid          => leptonic_get_cs2_one_grid
    eos_get_temp_one_grid         => leptonic_get_temp_one_grid
    eos_get_eps_range             => leptonic_get_eps_range
    eos_get_all_beta_eqm_one_grid => leptonic_get_all_beta_eqm_one_grid
    eos_eps_get_all_one_grid      => leptonic_eps_get_all_one_grid
    eos_temp_get_all_one_grid     => leptonic_temp_get_all_one_grid

    call eos_check

    if (mype == 0) then
      write(*,*) '****************************************'
      write(*,*) '-------- Leptonic EOS activated --------'
      write(*,*) ' baryon table   = ', trim(baryon_table_name)
      write(*,*) ' leptonic table = ', trim(leptonic_table_name)
      write(*,*) ' use muons      = ', lep_use_muons
      write(*,*) '****************************************'
    end if
  end subroutine eos_leptonic_activate

  pure subroutine lep_bound_ye_ymu(ye_in, ymu_in, yp_out, ye_out, ymu_out, logymu_out)
    double precision, intent(in)  :: ye_in, ymu_in
    double precision, intent(out) :: yp_out, ye_out, ymu_out, logymu_out

    ye_out = max(lep_eos_ylemin, min(lep_eos_ylemax, ye_in))
    ymu_out = max(lep_eos_ymumin, min(lep_eos_ymumax, ymu_in))

    if (lep_fix_ymu_high_yp) then
      if (ye_out + ymu_out > lep_eos_yemax) ymu_out = lep_eos_ymumin
    end if

    yp_out = max(lep_eos_yemin, min(lep_eos_yemax, ye_out + ymu_out))
    logymu_out = log(max(lep_eos_ymumin, ymu_out))
  end subroutine lep_bound_ye_ymu

  pure subroutine lep_bound_rho_temp(rho_in, temp_in, rho_out, temp_out, lrho_out, ltemp_out)
    double precision, intent(in)  :: rho_in, temp_in
    double precision, intent(out) :: rho_out, temp_out, lrho_out, ltemp_out

    rho_out = max(lep_eos_rhomin, min(lep_eos_rhomax, rho_in))
    temp_out = max(lep_eos_tempmin, min(lep_eos_tempmax, temp_in))
    lrho_out = log(rho_out)
    ltemp_out = log(temp_out)
  end subroutine lep_bound_rho_temp

  subroutine leptonic_eval_total_from_logs(lrho, ltemp, ye_in, ymu_in, eps, prs, ent, cs2, &
       mu_e, mu_n, mu_p, mu_mu, xa, xh, xn, xp, abar, zbar, muhat, munu)
    double precision, intent(in) :: lrho, ltemp, ye_in, ymu_in
    double precision, intent(out) :: eps, prs, ent, cs2
    double precision, intent(out) :: mu_e, mu_n, mu_p, mu_mu
    double precision, intent(out) :: xa, xh, xn, xp, abar, zbar
    double precision, intent(out) :: muhat, munu

    double precision :: yp, ye_local, ymu_local, logymu
    double precision :: ffx_b(lep_nvars_baryon)
    double precision :: ffx_e(lep_nvars_ele)
    double precision :: ffx_m(lep_nvars_muon)
    double precision :: press_ele, eps_ele, ent_ele
    double precision :: press_mu, eps_mu, ent_mu

    call lep_bound_ye_ymu(ye_in, ymu_in, yp, ye_local, ymu_local, logymu)

    call intep3d_lep_many(lrho, ltemp, yp, ffx_b, lep_tables_baryon, &
         lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon, &
         lep_logrho_table, lep_logtemp_table, lep_ye_table)

    press_ele = 0.0d0
    eps_ele = 0.0d0
    ent_ele = 0.0d0
    if (lep_add_ele_contribution) then
      call intep3d_lep_many(lrho, ltemp, ye_local, ffx_e, lep_tables_ele, &
           lep_nrho, lep_ntemp, lep_nyle, lep_nvars_ele, &
           lep_logrho_table, lep_logtemp_table, lep_yle_table)
      press_ele = ffx_e(i_lep_press_eminus) + ffx_e(i_lep_press_eplus)
      eps_ele = ffx_e(i_lep_eps_eminus) + ffx_e(i_lep_eps_eplus)
      ent_ele = ffx_e(i_lep_s_eminus) + ffx_e(i_lep_s_eplus)
      mu_e = ffx_e(i_lep_mu_ele)
    else
      mu_e = ffx_b(i_lep_mu_e)
    end if

    press_mu = 0.0d0
    eps_mu = 0.0d0
    ent_mu = 0.0d0
    if (lep_use_muons) then
      call intep3d_lep_many(lrho, ltemp, logymu, ffx_m, lep_tables_muon, &
           lep_nrho, lep_ntemp, lep_nymu, lep_nvars_muon, &
           lep_logrho_table, lep_logtemp_table, lep_logymu_table)
      press_mu = ffx_m(i_lep_press_muminus) + ffx_m(i_lep_press_muplus)
      eps_mu = ffx_m(i_lep_eps_muminus) + ffx_m(i_lep_eps_muplus)
      ent_mu = ffx_m(i_lep_s_muminus) + ffx_m(i_lep_s_muplus)
      mu_mu = ffx_m(i_lep_mu_mu)
    else
      mu_mu = 105.6583755d0
    end if

    prs = exp(ffx_b(i_lep_logpress)) + press_ele + press_mu
    eps = exp(ffx_b(i_lep_logenergy)) - lep_energy_shift + eps_ele + eps_mu
    ent = ffx_b(i_lep_entropy) + ent_ele + ent_mu
    cs2 = max(0.0d0, min(0.9999999d0, ffx_b(i_lep_cs2)))

    mu_p = ffx_b(i_lep_mu_p)
    mu_n = ffx_b(i_lep_mu_n)
    xa = ffx_b(i_lep_xa)
    xh = ffx_b(i_lep_xh)
    xn = ffx_b(i_lep_xn)
    xp = ffx_b(i_lep_xp)
    abar = ffx_b(i_lep_abar)
    zbar = ffx_b(i_lep_zbar)
    muhat = mu_n - mu_p - mu_e + Qnp
    munu = 0.0d0
  end subroutine leptonic_eval_total_from_logs

  subroutine leptonic_get_pressure_one_grid(prs, rho, eps, temp, ye, ymu)
    double precision, intent(inout) :: prs
    double precision, intent(in) :: rho, eps
    double precision, intent(inout), optional :: temp
    double precision, intent(in), optional :: ye, ymu

    double precision :: rho_local, eps_local, temp_local, ye_local

    ye_local = lep_eos_ylemin
    if (present(ye)) ye_local = ye
    eps_local = eps
    if (present(temp)) then
      temp_local = temp
    else
      temp_local = lep_eos_tempmin
      call leptonic_get_temp_one_grid(rho, eps_local, temp_local, ye_local, ymu)
    end if
    rho_local = rho

    call leptonic_temp_get_all_one_grid(rho_local, temp_local, ye_local, eps_local, prs=prs, ymu=ymu)
    if (present(temp)) temp = temp_local
  end subroutine leptonic_get_pressure_one_grid

  subroutine leptonic_get_eps_one_grid(prs, rho, eps, temp, ye, ymu)
    double precision, intent(in) :: prs, rho
    double precision, intent(inout) :: eps
    double precision, intent(in), optional :: temp, ye, ymu

    double precision :: rho_local, temp_local, ye_local

    if (.not. present(temp)) call mpistop('leptonic_get_eps_one_grid: need temp input')
    if (.not. present(ye)) call mpistop('leptonic_get_eps_one_grid: need ye input')

    rho_local = rho
    temp_local = temp
    ye_local = ye
    call leptonic_temp_get_all_one_grid(rho_local, temp_local, ye_local, eps, ymu=ymu)
  end subroutine leptonic_get_eps_one_grid

  subroutine leptonic_get_cs2_one_grid(cs2, rho, eps, temp, ye, ymu)
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho, eps
    double precision, intent(in), optional :: ye, ymu
    double precision, intent(inout), optional :: temp

    double precision :: rho_local, eps_local, temp_local, ye_local

    ye_local = lep_eos_ylemin
    if (present(ye)) ye_local = ye
    eps_local = eps
    if (present(temp)) then
      temp_local = temp
    else
      temp_local = lep_eos_tempmin
      call leptonic_get_temp_one_grid(rho, eps_local, temp_local, ye_local, ymu)
    end if
    rho_local = rho

    call leptonic_temp_get_all_one_grid(rho_local, temp_local, ye_local, eps_local, cs2=cs2, ymu=ymu)
    if (present(temp)) temp = temp_local
  end subroutine leptonic_get_cs2_one_grid

  subroutine leptonic_get_temp_one_grid(rho, eps, temp, ye, ymu)
    double precision, intent(in) :: rho, ye
    double precision, intent(inout) :: eps, temp
    double precision, intent(in), optional :: ymu

    double precision :: rho_local, temp_local, lrho, ltemp
    double precision :: ymu_local, eps_min, eps_max
    integer :: error_code

    ymu_local = lep_eos_ymumin
    if (present(ymu)) ymu_local = ymu

    call leptonic_get_eps_range(rho, eps_min, eps_max, ye, ymu_local)
    if (eps <= eps_min) then
      eps = eps_min
      temp = lep_eos_tempmin
      return
    end if
    if (eps >= eps_max) then
      eps = eps_max
      temp = lep_eos_tempmax
      return
    end if

    call lep_bound_rho_temp(rho, max(temp, lep_eos_tempmin), rho_local, temp_local, lrho, ltemp)
    ltemp = 0.5d0 * (lep_logtemp_table(1) + lep_logtemp_table(lep_ntemp))
    error_code = -1
    call rootfinding_brent(ltemp, lep_logtemp_table(1), lep_logtemp_table(lep_ntemp), &
         lep_eos_precision, lep_eos_iter_max, error_code, func_eps_of_temp)

    if (error_code == 2) call mpistop('leptonic_get_temp_one_grid: NaN in rootfinding')
    temp = exp(ltemp)

  contains
    double precision function func_eps_of_temp(ltemp_in)
      double precision, intent(in) :: ltemp_in
      double precision :: prs_tmp, ent_tmp, cs2_tmp
      double precision :: mu_e_tmp, mu_n_tmp, mu_p_tmp, mu_mu_tmp
      double precision :: xa_tmp, xh_tmp, xn_tmp, xp_tmp, abar_tmp, zbar_tmp
      double precision :: muhat_tmp, munu_tmp

      call leptonic_eval_total_from_logs(lrho, ltemp_in, ye, ymu_local, func_eps_of_temp, prs_tmp, &
           ent_tmp, cs2_tmp, mu_e_tmp, mu_n_tmp, mu_p_tmp, mu_mu_tmp, xa_tmp, xh_tmp, xn_tmp, xp_tmp, &
           abar_tmp, zbar_tmp, muhat_tmp, munu_tmp)
      func_eps_of_temp = eps - func_eps_of_temp
    end function func_eps_of_temp
  end subroutine leptonic_get_temp_one_grid

  subroutine leptonic_get_eps_range(rho, epsmin, epsmax, ye, ymu)
    double precision, intent(in) :: rho
    double precision, intent(out) :: epsmin, epsmax
    double precision, intent(in), optional :: ye, ymu

    double precision :: rho_local, temp_local, lrho, ltemp
    double precision :: ye_local, ymu_local
    double precision :: prs_tmp, ent_tmp, cs2_tmp
    double precision :: mu_e_tmp, mu_n_tmp, mu_p_tmp, mu_mu_tmp
    double precision :: xa_tmp, xh_tmp, xn_tmp, xp_tmp, abar_tmp, zbar_tmp
    double precision :: muhat_tmp, munu_tmp

    ye_local = lep_eos_ylemin
    ymu_local = lep_eos_ymumin
    if (present(ye)) ye_local = ye
    if (present(ymu)) ymu_local = ymu

    call lep_bound_rho_temp(rho, lep_eos_tempmin, rho_local, temp_local, lrho, ltemp)
    call leptonic_eval_total_from_logs(lrho, lep_logtemp_table(1), ye_local, ymu_local, epsmin, prs_tmp, &
         ent_tmp, cs2_tmp, mu_e_tmp, mu_n_tmp, mu_p_tmp, mu_mu_tmp, xa_tmp, xh_tmp, xn_tmp, xp_tmp, &
         abar_tmp, zbar_tmp, muhat_tmp, munu_tmp)
    call leptonic_eval_total_from_logs(lrho, lep_logtemp_table(lep_ntemp), ye_local, ymu_local, epsmax, prs_tmp, &
         ent_tmp, cs2_tmp, mu_e_tmp, mu_n_tmp, mu_p_tmp, mu_mu_tmp, xa_tmp, xh_tmp, xn_tmp, xp_tmp, &
         abar_tmp, zbar_tmp, muhat_tmp, munu_tmp)
  end subroutine leptonic_get_eps_range

  subroutine leptonic_temp_get_all_one_grid(rho,temp,ye,eps,prs,ent,cs2,dedt,&
       dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,mu_mu,muhat,munu,ymu)
    double precision, intent(inout) :: rho, eps, temp
    double precision, intent(in) :: ye
    double precision, intent(in), optional :: ymu
    double precision, intent(inout), optional :: prs,ent,cs2,dedt
    double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
    double precision, intent(inout), optional :: mu_e,mu_n,mu_p,mu_mu,muhat,munu

    double precision :: rho_local, temp_local, lrho, ltemp
    double precision :: ymu_local
    double precision :: prs_local, ent_local, cs2_local
    double precision :: mu_e_local, mu_n_local, mu_p_local, mu_mu_local
    double precision :: xa_local, xh_local, xn_local, xp_local, abar_local, zbar_local
    double precision :: muhat_local, munu_local

    if (rho <= small_rho_thr) then
      rho = small_rho
      eps = small_eps
      temp = max(small_temp, lep_eos_tempmin)
      if (present(prs)) call atmo_get_pressure_one_grid(prs, rho, eps)
      if (present(cs2)) call atmo_get_cs2_one_grid(cs2, rho, eps)
      if (present(ent)) ent = 0.0d0
      if (present(xa)) xa = 0.0d0
      if (present(xh)) xh = 0.0d0
      if (present(xn)) xn = 0.0d0
      if (present(xp)) xp = 0.0d0
      if (present(abar)) abar = 0.0d0
      if (present(zbar)) zbar = 0.0d0
      if (present(mu_e)) mu_e = 0.0d0
      if (present(mu_n)) mu_n = 0.0d0
      if (present(mu_p)) mu_p = 0.0d0
      if (present(mu_mu)) mu_mu = 105.6583755d0
      if (present(muhat)) muhat = 0.0d0
      if (present(munu)) munu = 0.0d0
      return
    end if

    ymu_local = lep_eos_ymumin
    if (present(ymu)) ymu_local = ymu

    call lep_bound_rho_temp(rho, temp, rho_local, temp_local, lrho, ltemp)
    call leptonic_eval_total_from_logs(lrho, ltemp, ye, ymu_local, eps, prs_local, ent_local, cs2_local, &
         mu_e_local, mu_n_local, mu_p_local, mu_mu_local, xa_local, xh_local, xn_local, xp_local, &
         abar_local, zbar_local, muhat_local, munu_local)

    rho = rho_local
    temp = temp_local

    if (present(prs)) prs = prs_local
    if (present(ent)) ent = ent_local
    if (present(cs2)) cs2 = cs2_local
    if (present(dedt)) dedt = 0.0d0
    if (present(dpderho)) dpderho = 0.0d0
    if (present(dpdrhoe)) dpdrhoe = 0.0d0
    if (present(xa)) xa = xa_local
    if (present(xh)) xh = xh_local
    if (present(xn)) xn = xn_local
    if (present(xp)) xp = xp_local
    if (present(abar)) abar = abar_local
    if (present(zbar)) zbar = zbar_local
    if (present(mu_e)) mu_e = mu_e_local
    if (present(mu_n)) mu_n = mu_n_local
    if (present(mu_p)) mu_p = mu_p_local
    if (present(mu_mu)) mu_mu = mu_mu_local
    if (present(muhat)) muhat = muhat_local
    if (present(munu)) munu = munu_local
  end subroutine leptonic_temp_get_all_one_grid

  subroutine leptonic_eps_get_all_one_grid(rho,eps,ye,temp,prs,ent,cs2,dedt,&
       dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,mu_mu,muhat,munu,ymu)
    double precision, intent(inout) :: rho, eps, ye
    double precision, intent(in), optional :: ymu
    double precision, intent(inout), optional :: ent, prs, temp, cs2, dedt
    double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
    double precision, intent(inout), optional :: mu_e,mu_n,mu_p,mu_mu,muhat,munu

    double precision :: temp_local

    if (.not. present(temp)) then
      temp_local = lep_eos_tempmin
    else
      temp_local = temp
    end if

    call leptonic_get_temp_one_grid(rho, eps, temp_local, ye, ymu)
    call leptonic_temp_get_all_one_grid(rho, temp_local, ye, eps, prs=prs, ent=ent, cs2=cs2, &
         dedt=dedt, dpderho=dpderho, dpdrhoe=dpdrhoe, xa=xa, xh=xh, xn=xn, xp=xp, abar=abar, zbar=zbar, &
         mu_e=mu_e, mu_n=mu_n, mu_p=mu_p, mu_mu=mu_mu, muhat=muhat, munu=munu, ymu=ymu)

    if (present(temp)) temp = temp_local
  end subroutine leptonic_eps_get_all_one_grid

  double precision function leptonic_find_ye_of_mumu(lrho, ltemp, mu_mu)
    double precision, intent(in) :: lrho, ltemp, mu_mu
    integer :: error_code

    leptonic_find_ye_of_mumu = 0.5d0 * (lep_eos_ylemin + lep_eos_ylemax)
    error_code = -1
    call rootfinding_brent(leptonic_find_ye_of_mumu, lep_eos_ylemin, lep_eos_ylemax, &
         lep_eos_precision, lep_eos_iter_max, error_code, func_yle)
    if (error_code /= 0) then
      leptonic_find_ye_of_mumu = max(lep_eos_ylemin, min(lep_eos_ylemax, leptonic_find_ye_of_mumu))
    end if

  contains
    double precision function func_yle(ye_try)
      double precision, intent(in) :: ye_try
      double precision :: ffx_e(lep_nvars_ele)

      call intep3d_lep_many(lrho, ltemp, ye_try, ffx_e, lep_tables_ele, &
           lep_nrho, lep_ntemp, lep_nyle, lep_nvars_ele, &
           lep_logrho_table, lep_logtemp_table, lep_yle_table)
      func_yle = ffx_e(i_lep_mu_ele) - mu_mu
    end function func_yle
  end function leptonic_find_ye_of_mumu

  double precision function leptonic_find_betaeq_ye(lrho, ltemp)
    double precision, intent(in) :: lrho, ltemp
    integer :: error_code

    leptonic_find_betaeq_ye = 0.5d0 * (lep_eos_yemin + lep_eos_yemax)
    error_code = -1
    call rootfinding_brent(leptonic_find_betaeq_ye, lep_eos_yemin, lep_eos_yemax, &
         lep_eos_precision, lep_eos_iter_max, error_code, func_betaeq_ye)
    if (error_code /= 0) then
      leptonic_find_betaeq_ye = max(lep_eos_yemin, min(lep_eos_yemax, leptonic_find_betaeq_ye))
    end if

  contains
    double precision function func_betaeq_ye(ye_try)
      double precision, intent(in) :: ye_try
      double precision :: ffx_b(lep_nvars_baryon)

      call intep3d_lep_many(lrho, ltemp, ye_try, ffx_b, lep_tables_baryon, &
           lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon, &
           lep_logrho_table, lep_logtemp_table, lep_ye_table)
      func_betaeq_ye = ffx_b(i_lep_mu_e) + ffx_b(i_lep_mu_p) - ffx_b(i_lep_mu_n) - Qnp
    end function func_betaeq_ye
  end function leptonic_find_betaeq_ye

  subroutine leptonic_get_all_beta_eqm_one_grid(rho, temp, ye, eps, prs, ymu)
    double precision, intent(in) :: rho, temp
    double precision, intent(inout) :: ye
    double precision, intent(inout), optional :: eps, prs
    double precision, intent(inout), optional :: ymu

    double precision :: rho_local, temp_local, lrho, ltemp
    double precision :: ymu_local, ye_local
    double precision :: eps_local, prs_local
    double precision :: mu_mu_root
    integer :: error_code

    call lep_bound_rho_temp(rho, temp, rho_local, temp_local, lrho, ltemp)

    if (.not. lep_use_muons .or. .not. present(ymu)) then
      ye = leptonic_find_betaeq_ye(lrho, ltemp)
      if (present(ymu)) ymu = lep_eos_ymumin
    else
      ymu_local = 0.5d0 * (log(lep_eos_ymumin) + log(lep_eos_ymumax))
      error_code = -1
      call rootfinding_brent(ymu_local, log(lep_eos_ymumin), log(lep_eos_ymumax), &
           lep_eos_precision, lep_eos_iter_max, error_code, func_betaeq_ymu)

      if (ymu_local <= log(lep_eos_ymumin) + 1.0d-10 .or. error_code /= 0) then
        ymu_local = lep_eos_ymumin
        ye_local = leptonic_find_betaeq_ye(lrho, ltemp)
      else
        ymu_local = exp(ymu_local)
        ye_local = leptonic_find_ye_of_mumu(lrho, ltemp, mu_mu_root)
      end if

      ye = ye_local
      ymu = ymu_local
    end if

    if (present(eps) .or. present(prs)) then
      ye_local = ye
      if (present(ymu)) then
        ymu_local = ymu
      else
        ymu_local = lep_eos_ymumin
      end if
      call leptonic_temp_get_all_one_grid(rho_local, temp_local, ye_local, eps_local, prs=prs_local, ymu=ymu_local)
      if (present(eps)) eps = eps_local
      if (present(prs)) prs = prs_local
    end if

  contains
    double precision function func_betaeq_ymu(logymu_try)
      double precision, intent(in) :: logymu_try
      double precision :: ffx_b(lep_nvars_baryon)
      double precision :: ffx_m(lep_nvars_muon)
      double precision :: yp_try

      call intep3d_lep_many(lrho, ltemp, logymu_try, ffx_m, lep_tables_muon, &
           lep_nrho, lep_ntemp, lep_nymu, lep_nvars_muon, &
           lep_logrho_table, lep_logtemp_table, lep_logymu_table)
      mu_mu_root = ffx_m(i_lep_mu_mu)
      ye_local = leptonic_find_ye_of_mumu(lrho, ltemp, mu_mu_root)
      yp_try = max(lep_eos_yemin, min(lep_eos_yemax, ye_local + exp(logymu_try)))
      call intep3d_lep_many(lrho, ltemp, yp_try, ffx_b, lep_tables_baryon, &
           lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon, &
           lep_logrho_table, lep_logtemp_table, lep_ye_table)
      func_betaeq_ymu = mu_mu_root + ffx_b(i_lep_mu_p) - ffx_b(i_lep_mu_n) - Qnp
    end function func_betaeq_ymu
  end subroutine leptonic_get_all_beta_eqm_one_grid

end module mod_eos_leptonic
