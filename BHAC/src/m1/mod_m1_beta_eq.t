module mod_m1_beta_eq
  use mod_m1_internal
  use mod_m1_eas_param
  use mod_m1_fermi
  use mod_m1_constants
  use mod_eos, only: eos_rhomin, eos_rhomax, eos_yemin, eos_yemax, eos_tempmin, eos_tempmax, &
       eos_ymumin, eos_temp_get_all_one_grid, eos_has_ymu
  use mod_rootfinding, only: rootfinding_global_multid_newton_raphson, numerical_jacobian

  implicit none
  !KEN changed this completely

  ! ====================================================================
  ! CORRECTED CONSTANTS - Must match mod_m1_constants and C++ values
  ! ====================================================================
  double precision :: rho0 = 1.0d+11           ! in g/cm^3 (density suppression threshold)
  double precision :: mb = MNUC_CGS            ! nucleon mass from mod_m1_constants
  double precision :: mb_msun = MNUC_MSUN            ! nucleon mass from mod_m1_constants
  double precision :: kB = 8.61738568d-11              ! Boltzmann constant from mod_m1_constants
  double precision :: Pi_const = PI            ! Pi from mod_m1_constants
  double precision :: hc_const = HC_MEVCM      ! ℏc from mod_m1_constants
  
  integer, parameter :: m1_beta_ndim = 2

  type m1_beta_eq_helpers
    double precision :: eta_nu_e, eta_nu_ebar, eta_nu_x, eta_nu_mu, eta_nu_mubar
    double precision :: J_nu_e, J_nu_ebar, J_nu_x, J_nu_mu, J_nu_mubar
    double precision :: fluid_prim_rho, fluid_prim_ye, fluid_prim_ymu, fluid_prim_T, fluid_prim_eps
    double precision :: N_KSP1,N_KSP2
  end type m1_beta_eq_helpers
  
  type(m1_beta_eq_helpers) :: stateEQ
  
contains
  !> Calculate the correct beta-equilibrium including neutrinos and update T_fluid, Y_e_fluid
  subroutine m1_get_beta_equilibrium(wrad_KSP,N_KSP,fluid_Prim,Jrad_closure,Gamma_closure,ix^D,T_eq,Ye_eq)
    integer, intent(in) :: ix^D
    double precision, intent(in) :: fluid_Prim(1:fluid_vars)
    double precision, intent(in) :: wrad_KSP(1:m1_numvars_internal,^NS)
    double precision, intent(in) :: N_KSP(^NS),Jrad_closure(^NS),Gamma_closure(^NS)  
    double precision, intent(out) :: T_eq,Ye_eq
    !internal:
    integer :: error
    double precision :: z_vec(1:2)
    double precision :: eps,prs,ent,cs2,dedt
    double precision :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,mu_mu,muhat,munu
    double precision :: T_fluid, Ye_fluid, Ymu_fluid, rho_fluid

    stateEQ%fluid_prim_rho = fluid_Prim(idx_rho)
    stateEQ%fluid_prim_Ye = fluid_Prim(idx_ye)
    stateEQ%fluid_prim_T = fluid_Prim(idx_T)
    stateEQ%fluid_prim_ymu = eos_ymumin
    if (m1_i_mu > 0) stateEQ%fluid_prim_ymu = fluid_Prim(idx_ymu)

    rho_fluid = stateEQ%fluid_prim_rho
    Ye_fluid = stateEQ%fluid_prim_Ye
    T_fluid = stateEQ%fluid_prim_T
    Ymu_fluid = stateEQ%fluid_prim_ymu
    
    if(Ye_fluid > eos_yemax) write(*,*) "how"
    if (eos_has_ymu()) then
      call eos_temp_get_all_one_grid(rho_fluid,T_fluid,Ye_fluid,eps,prs=prs,ent=ent,&
      cs2=cs2,dedt=dedt,dpderho=dpderho,dpdrhoe=dpdrhoe,xa=xa,xh=xh,xn=xn,xp=xp,&
      abar=abar,zbar=zbar,mu_e=mu_e,mu_n=mu_n,mu_p=mu_p,mu_mu=mu_mu,muhat=muhat,munu=munu,&
      ymu=Ymu_fluid)
    else
      call eos_temp_get_all_one_grid(rho_fluid,T_fluid,Ye_fluid,eps,prs=prs,ent=ent,&
      cs2=cs2,dedt=dedt,dpderho=dpderho,dpdrhoe=dpdrhoe,xa=xa,xh=xh,xn=xn,xp=xp,&
      abar=abar,zbar=zbar,mu_e=mu_e,mu_n=mu_n,mu_p=mu_p,muhat=muhat,munu=munu)
      mu_mu = 0.0d0
    end if

    stateEQ%fluid_prim_eps = eps 

    !> fluid frame J
    stateEQ%J_nu_e = 0.0d0
    stateEQ%J_nu_ebar = 0.0d0
    stateEQ%J_nu_x = 0.0d0
    stateEQ%J_nu_mu = 0.0d0
    stateEQ%J_nu_mubar = 0.0d0
    if (m1_i_nue > 0) stateEQ%J_nu_e = Jrad_closure(m1_i_nue)
    if (m1_i_nuebar > 0) stateEQ%J_nu_ebar = Jrad_closure(m1_i_nuebar)
    if (m1_i_nux > 0) stateEQ%J_nu_x = Jrad_closure(m1_i_nux)
    if (m1_i_mu > 0) stateEQ%J_nu_mu = Jrad_closure(m1_i_mu)
    if (m1_i_mubar > 0) stateEQ%J_nu_mubar = Jrad_closure(m1_i_mubar)

    stateEQ%N_KSP1 = 0.0d0
    stateEQ%N_KSP2 = 0.0d0
    if (m1_i_nue > 0) stateEQ%N_KSP1 = N_KSP(m1_i_nue)
    if (m1_i_nuebar > 0) stateEQ%N_KSP2 = N_KSP(m1_i_nuebar)

    !> first guess:
    z_vec(1) = stateEQ%fluid_prim_Ye !Ye_eq
    z_vec(2) = stateEQ%fluid_prim_T  !T_eq

    call betaeq_rootfinding_bounded_newton_raphson(m1_beta_ndim, z_vec, 1.0d-11, 100, &
    error, m1_beta_equilib_function, m1_beta_equilib_jacobian, m1_beta_equilib_min )
    
    if(error .ne. 0) then
      write(*,*) "beta-equilibrium: root not found at ix^D, error:",error
      write(*,*) "  rho=", rho_fluid, " Ye=", Ye_fluid, " T=", T_fluid
        
      ! NEW: Evaluate residuals at final point
      block
        double precision :: final_residuals(2)
        final_residuals = m1_beta_equilib_function(z_vec)
        write(*,*) "  Final residuals: f1=", final_residuals(1), " f2=", final_residuals(2)
      end block
      
      !write(*,*) "  Falling back to current state"
      ! Fallback to current fluid state
      Ye_eq = stateEQ%fluid_prim_Ye
      T_eq = stateEQ%fluid_prim_T
      return
      else 
          !write(*,*) "We passed"
          !write(*,*) "  rho=", rho_fluid, " Ye=", Ye_fluid, " T=", T_fluid
        block
        double precision :: final_residuals(2)
        final_residuals = m1_beta_equilib_function(z_vec)
        !write(*,*) "  Final residuals: f1=", final_residuals(1), " f2=", final_residuals(2)
      end block

    end if


    !> new equilibrium values
    Ye_eq = z_vec(1)
    T_eq = z_vec(2)

    ! Final safety clamp (should rarely be needed with bounded NR)
    if (Ye_eq .ge. eos_Yemax) then
      Ye_eq=eos_Yemax*0.99999999
   else if (Ye_eq .le. eos_Yemin) then
      Ye_eq=eos_Yemin*1.000000001
   endif
   if (T_eq .ge. eos_tempmax) then
      T_eq=eos_tempmax*0.99999999
   else if (T_eq .le. eos_tempmin) then
      T_eq=eos_tempmin*1.00000001
   endif 

  contains
  !> Equations for finding root T_Eq, Ye_Eq
  function m1_beta_equilib_function(z_vec)
    implicit none
    double precision, dimension(:), intent(in) :: z_vec
    double precision, dimension(1:size(z_vec)) :: m1_beta_equilib_function
    ! internal
    double precision, parameter :: Qnp = 1.2935d0 ! MeV - rest-mass energy difference n&p
    double precision, parameter :: PENALTY = 1.0d30
    double precision :: Yl,Y_nu,Y_nue,Y_nue_bar,Y_rhs
    double precision :: e,u,eps_internal,chem_pot
    double precision :: Z_nu,Z_nue,Z_nue_bar,Z_nux,Z_numu,Z_numubar
    double precision :: eps,prs,ent,cs2,dedt
    double precision :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,mu_mu,muhat,munu
    double precision :: eta_nu(^NS),Tequil,Yeequil
    double precision :: rho_cgs

    Tequil = z_vec(2)
    Yeequil = z_vec(1)

    !> Check if out of bounds and return penalty
    if (Yeequil .le. eos_Yemin .or. Yeequil .ge. eos_Yemax .or. &
        Tequil .le. eos_tempmin .or. Tequil .ge. eos_tempmax) then
        m1_beta_equilib_function(1) = PENALTY
        m1_beta_equilib_function(2) = PENALTY
        return
    endif

    if (eos_has_ymu()) then
      call eos_temp_get_all_one_grid(stateEQ%fluid_prim_rho,Tequil,Yeequil,eps,prs=prs,ent=ent,&
      cs2=cs2,dedt=dedt,dpderho=dpderho,dpdrhoe=dpdrhoe,xa=xa,xh=xh,xn=xn,xp=xp,&
      abar=abar,zbar=zbar,mu_e=mu_e,mu_n=mu_n,mu_p=mu_p,mu_mu=mu_mu,muhat=muhat,munu=munu,&
      ymu=stateEQ%fluid_prim_ymu)
    else
      call eos_temp_get_all_one_grid(stateEQ%fluid_prim_rho,Tequil,Yeequil,eps,prs=prs,ent=ent,&
      cs2=cs2,dedt=dedt,dpderho=dpderho,dpdrhoe=dpdrhoe,xa=xa,xh=xh,xn=xn,xp=xp,&
      abar=abar,zbar=zbar,mu_e=mu_e,mu_n=mu_n,mu_p=mu_p,muhat=muhat,munu=munu)
      mu_mu = 0.0d0
    end if

    ! Compute eta_nu from chemical potentials
    {^KSP&
    if(^KSP .eq. m1_i_nue) then
      chem_pot = mu_e + mu_p - mu_n - Qnp
    else if(^KSP .eq. m1_i_nuebar) then
      chem_pot = -1.0d0*(mu_e + mu_p - mu_n - Qnp)
    else if(^KSP .eq. m1_i_nux) then
      chem_pot = 0.0d0
    else if(^KSP .eq. m1_i_mu) then
      chem_pot = mu_mu + mu_p - mu_n - Qnp
    else if(^KSP .eq. m1_i_mubar)then
      chem_pot = -1.0d0*(mu_mu + mu_p - mu_n - Qnp)
    else
      chem_pot = 0.0d0
    end if 
    eta_nu(^KSP) = chem_pot/Tequil
    \}

    stateEQ%eta_nu_e = 0.0d0
    stateEQ%eta_nu_ebar = 0.0d0
    stateEQ%eta_nu_x = 0.0d0
    stateEQ%eta_nu_mu = 0.0d0
    stateEQ%eta_nu_mubar = 0.0d0
    if (m1_i_nue > 0) stateEQ%eta_nu_e = eta_nu(m1_i_nue)
    if (m1_i_nuebar > 0) stateEQ%eta_nu_ebar = eta_nu(m1_i_nuebar)
    if (m1_i_nux > 0) stateEQ%eta_nu_x = eta_nu(m1_i_nux)
    if (m1_i_mu > 0) stateEQ%eta_nu_mu = eta_nu(m1_i_mu)
    if (m1_i_mubar > 0) stateEQ%eta_nu_mubar = eta_nu(m1_i_mubar)

    ! ====================================================================
    ! CRITICAL: Convert density from geometrized units to CGS
    ! ====================================================================
    ! stateEQ%fluid_prim_rho is in GEOMETRIZED units (code units)
    ! We need CGS for the neutrino density calculations
    rho_cgs = stateEQ%fluid_prim_rho * INVRHOGF
    
    ! ====================================================================
    ! Equation (91): Lepton number conservation
    ! Y_l = Y_e + (m_b/ρ) × [N_νe - N_ν̄e]
    ! ====================================================================
    Yl = stateEQ%fluid_prim_ye + mb_msun / stateEQ%fluid_prim_rho * (stateEQ%N_KSP1 - stateEQ%N_KSP2)
    
    ! Neutrino number density prefactor (dimensionless, per baryon)
    ! Y_ν = 4π/(ℏc)³ × (m_b/ρ) × (k_B T)³ × exp(-ρ₀/ρ)
    ! This matches C++ implementation:
    ! Ynu = 4*pi / hc_mevcm³ * temp³ * FD[2+n](eta) * exp(-rho_lim/F.rho) * mnuc_cgs/F.rho
    Y_nu = 4.0d0 * Pi_const / (hc_const**3) * mb_msun / stateEQ%fluid_prim_rho * &
           (Tequil)**3 * exp(-rho0 / rho_cgs)
    
    if(m1_i_nue > 0) then
      Y_nue = Y_nu * fermi_dirac(stateEQ%eta_nu_e, 2)
    else
      Y_nue = 0.0d0
    end if
    if(m1_i_nuebar > 0) then
      Y_nue_bar = Y_nu * fermi_dirac(stateEQ%eta_nu_ebar, 2)
    else 
      Y_nue_bar = 0.0d0
    end if 

    Y_rhs = Yeequil + Y_nue - Y_nue_bar

    ! Residual 1: Lepton number conservation
    m1_beta_equilib_function(1) = Y_rhs - Yl 
    
    ! ====================================================================
    ! Equation (92): Energy conservation
    ! u = e(T,Ye) + (ρ/m_b) × [Z_νe + Z_ν̄e + 4Z_νx]
    ! ====================================================================
    
    ! Current total energy density (with old temperature)
    eps_internal = stateEQ%fluid_prim_eps
    !e = (1.0d0 + stateEQ%fluid_prim_rho) * eps_internal
    e = (1.0d0 + eps_internal) * stateEQ%fluid_prim_rho
    u = e + stateEQ%J_nu_e + stateEQ%J_nu_ebar + stateEQ%J_nu_x + &
        stateEQ%J_nu_mu + stateEQ%J_nu_mubar

    ! New matter energy density (with equilibrium temperature)
    eps_internal = eps
    !e = (1.0d0 + stateEQ%fluid_prim_rho) * eps_internal
    e = (1.0d0 + eps_internal) * stateEQ%fluid_prim_rho !KEN fixed e calculation

    ! Neutrino energy density prefactor
    ! Z_ν = 4π/(ℏc)³ × (m_b/ρ) × (k_B T)⁴ × exp(-ρ₀/ρ)
    ! This matches C++ implementation:
    ! Ynu = 4*pi / hc_mevcm³ * temp⁴ * FD[3+n](eta) * exp(-rho_lim/F.rho) * mnuc_cgs/F.rho
    !Z_nu = 4.0d0 * Pi_const / (hc_const**3) * mb / rho_cgs * &
    !       (Tequil)**4 * exp(-rho0 / rho_cgs)
    Z_nu = 4.0d0 * Pi_const / (hc_const**3) * (Tequil)**4 * exp(-rho0 / rho_cgs)*MEV_TO_ERG_KEN * EPSGF * RHOGF 
    
    if(m1_i_nue > 0) then
      Z_nue = Z_nu * fermi_dirac(stateEQ%eta_nu_e, 3)
    else
      Z_nue = 0.0d0
    end if
    if(m1_i_nuebar > 0) then
      Z_nue_bar = Z_nu * fermi_dirac(stateEQ%eta_nu_ebar, 3)
    else
      Z_nue_bar = 0.0d0
    end if 
    if(m1_i_nux > 0) then
      Z_nux = Z_nu * fermi_dirac(stateEQ%eta_nu_x, 3)
    else
      Z_nux = 0.0d0
    end if
    if(m1_i_mu > 0) then
      Z_numu = Z_nu * fermi_dirac(stateEQ%eta_nu_mu, 3)
    else
      Z_numu = 0.0d0
    end if
    if(m1_i_mubar > 0) then
      Z_numubar = Z_nu * fermi_dirac(stateEQ%eta_nu_mubar, 3)
    else
      Z_numubar = 0.0d0
    end if 
    
    ! Residual 2: Energy conservation
    ! Note: ρ/m_b = n_b (baryon number density in CGS)
    !m1_beta_equilib_function(2) = -u + e + rho_cgs/mb * (Z_nue + Z_nue_bar + 4.0d0*Z_nux)
    if ((m1_i_mu > 0) .or. (m1_i_mubar > 0)) then
      m1_beta_equilib_function(2) = -u + e + (Z_nue + Z_nue_bar + Z_numu + Z_numubar + 2.0d0*Z_nux)
    else
      m1_beta_equilib_function(2) = -u + e + (Z_nue + Z_nue_bar + 4.0d0*Z_nux)
    end if

  end function m1_beta_equilib_function

  !> Jacobian
  function m1_beta_equilib_jacobian(z_vec)
    implicit none
    double precision, dimension(:), intent(in) :: z_vec
    double precision, dimension(1:size(z_vec),1:size(z_vec)) :: m1_beta_equilib_jacobian
    double precision :: jac(1:size(z_vec),1:size(z_vec)) 

    call numerical_jacobian(m1_beta_equilib_function, z_vec, jac)
    m1_beta_equilib_jacobian(:,:) = jac(:,:)

  end function m1_beta_equilib_jacobian

  !> Minimal values
  function m1_beta_equilib_min(z_vec)
    implicit none   
    double precision, dimension(:), intent(in) :: z_vec
    double precision :: m1_beta_equilib_min

    m1_beta_equilib_min = -1.0d+40

  end function m1_beta_equilib_min
 
!> Bounded Newton-Raphson with overshoot detection
subroutine betaeq_rootfinding_bounded_newton_raphson(n_dim, z, tolerance, iter_max, &
  return_code, fvec_in, fjac_in, fmin_in)
  use mod_lu
  use mod_lnsrch
  integer, intent(in)          :: n_dim
  double precision, dimension(1:n_dim), intent(inout):: z
  double precision, intent(in) :: tolerance
  integer, intent(in)          :: iter_max
  integer, intent(out) :: return_code
  interface
     function fvec_in(z_vec)
        implicit none
        double precision, dimension(:), intent(in)  :: z_vec
        double precision, dimension(1:size(z_vec))    :: fvec_in
     end function fvec_in
     function fjac_in(z_vec)
        implicit none
        double precision, dimension(:), intent(in)  :: z_vec
        double precision, dimension(1:size(z_vec),1:size(z_vec)) :: fjac_in
     end function fjac_in
     function fmin_in(z_vec)
        implicit none
        double precision :: fmin_in
        double precision, dimension(:), intent(in)  :: z_vec
     end function fmin_in
  end interface
  
  double precision, parameter :: bound_safety = 1.0d-6
  double precision, parameter :: STPMX = 1.0d3
  integer                     :: it_root, consecutive_clamps
  integer                     :: err_flag, d_flag
  integer, dimension(1:n_dim) :: indx
  double precision, dimension(1:n_dim)          :: fvec, dz, zold, g, z_before_clamp
  double precision, dimension(1:n_dim, 1:n_dim) :: jac
  double precision :: ye_min_safe, ye_max_safe, t_min_safe, t_max_safe
  double precision :: f, fold, stpmax
  logical          :: was_clamped
  
  return_code = -1
  consecutive_clamps = 0
  dz(:) = huge(0.0d0)
  stpmax = STPMX * max(dot_product(z(:),z(:)), dble(n_dim))
  
  ! Define safe boundaries
  ye_min_safe = eos_Yemin + bound_safety
  ye_max_safe = eos_Yemax - bound_safety
  t_min_safe = eos_tempmin + bound_safety
  t_max_safe = eos_tempmax - bound_safety
  
  ! Clamp initial guess
  z(1) = max(ye_min_safe, min(ye_max_safe, z(1)))
  z(2) = max(t_min_safe, min(t_max_safe, z(2)))
  
  f = fmin_in(z)
  
  do it_root = 1, iter_max
     fvec = fvec_in(z)
  
     ! Check convergence on residuals
     if (maxval(abs(fvec)) <= tolerance) then
        return_code = 0
        return
     end if
     
     ! Check if step size has become negligible (stagnation)
     if (maxval(abs(dz)) < tolerance) then
        return_code = 0
        return
     end if
  
     jac(:,:) = fjac_in(z(:))
     g = matmul(transpose(jac), fvec)
     zold(:) = z(:)
     fold = f
     dz(:) = -fvec(:)
  
     call ludcmp(jac, n_dim, indx, d_flag, err_flag)
     if (err_flag == 1) then
        return_code = 4
        return
     end if
     call lubksb(jac, n_dim, indx, dz)
     
     ! Line search to find better step
     call lnsrch(zold, fold, g, dz, z, f, stpmax, err_flag, fmin_in)
     
     if (err_flag == 1) then
        return_code = 3  ! Line search failed
        return
     else if (err_flag == 2) then
        return_code = 5  ! Roundoff problem in line search
        return
     end if
     
     ! Check if NaN before clamping
     if (any(isnan(z(:)))) then
        return_code = 2
        return
     end if
     
     ! Store state before clamping
     z_before_clamp(:) = z(:)
     was_clamped = .false.
     
     ! Clamp to bounds
     if (z(1) < ye_min_safe .or. z(1) > ye_max_safe) was_clamped = .true.
     if (z(2) < t_min_safe .or. z(2) > t_max_safe) was_clamped = .true.
     
     z(1) = max(ye_min_safe, min(ye_max_safe, z(1)))
     z(2) = max(t_min_safe, min(t_max_safe, z(2)))
     
     ! Track consecutive clamping events
     if (was_clamped) then
        consecutive_clamps = consecutive_clamps + 1
        
        ! If we've clamped twice in a row, the step becomes ~0 and we're stuck at boundary
        if (consecutive_clamps >= 2) then
           ! Check if we're actually converged despite being at boundary
           fvec = fvec_in(z)
           if (maxval(abs(fvec)) <= tolerance * 10.0d0) then
              return_code = 0  ! Good enough
              return
           else
              return_code = 6  ! Stuck at boundary, not converged
              return
           end if
        end if
     else
        consecutive_clamps = 0  ! Reset counter if no clamping
     end if
     
     ! Update function value after clamping (if clamped)
     if (was_clamped) then
        f = fmin_in(z)
     end if
  
  end do
  
  return_code = 1  ! Max iterations reached

end subroutine betaeq_rootfinding_bounded_newton_raphson
  
end subroutine m1_get_beta_equilibrium

end module mod_m1_beta_eq
