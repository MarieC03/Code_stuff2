module mod_m1_eas
  use mod_m1_eas_param
  use mod_m1_metric_interface
  use mod_m1_closure
  implicit none

  public

  double precision :: timeCurrent
  procedure(sub_m1_get_eas), pointer :: m1_get_eas => null()
  
  abstract interface

     subroutine sub_m1_get_eas(wrad,speciesKSP,ix^D,eas_eq,fluid_Prim) !T_fluid,rho_fluid,ye_fluid)
      use mod_m1_internal
      use mod_m1_eas_param
      {#IFDEF UNIT_TESTS
      use mod_m1_tests
      }
      {#IFNDEF UNIT_TESTS
       include "amrvacdef.f"
       }
       integer, intent(in) :: ix^D, speciesKSP
       double precision, intent(in)    :: wrad(1:m1_numvars_internal)
       double precision, intent(inout) :: fluid_Prim(1:fluid_vars)
       double precision, intent(inout) :: eas_eq(1:m1_num_eas)
     end subroutine sub_m1_get_eas

  end interface
  
contains

  subroutine m1_eas_gray_activate()
    use mod_m1_eas_param
    use mod_m1_eas_microphysical_gray
    {#IFNDEF TABEOS
      {#IFNDEF UNIT_TESTS
       call mpistop("You need 3d tables for microphysical eas")
      }
      {#IFDEF UNIT_TESTS
        write(99,*) "You need 3d tables for microphysical eas"
      }
    }
    m1_get_eas => m1_get_eas_microphysical_gray_Weakhub
  
  end subroutine m1_eas_gray_activate

  subroutine m1_eas_activate()
    if( .not. associated(m1_get_eas) ) then
      {#IFDEF UNIT_TESTS
           write(99,*)"Error: no microphysics implementation provided for M1."
      }
      {#IFNDEF UNIT_TESTS
       call mpistop("Error: no microphysics implementation provided for M1.")
       }
    end if
  end subroutine m1_eas_activate

 !********************************
  subroutine m1_allocate_eas_ixD()
    
    if(allocated(m1_eas_ixD)) deallocate(m1_eas_ixD)
    allocate( m1_eas_ixD(m1_num_eas,^NS) )
    m1_eas_ixD = 0.0d0
    
  end subroutine m1_allocate_eas_ixD


  subroutine m1_deallocate_eas_ixD()
    
    if(allocated(m1_eas_ixD)) deallocate( m1_eas_ixD )
    
  end subroutine m1_deallocate_eas_ixD

  !********************************not used:
  subroutine m1_allocate_eas(ix^L)
    integer, intent(in) :: ix^L
    
    allocate( m1_eas(ix^S,m1_num_eas,^NS) )
    
  end subroutine m1_allocate_eas


  subroutine m1_deallocate_eas()
    
    deallocate( m1_eas )
    
  end subroutine m1_deallocate_eas

 !****************************************************************************
  subroutine m1_add_eas(qdt,qtC,igrid,x,wprim,wradimpl,ixI^L,calc_eqRates)
    !Note: have to make sure that wprim is prim since we need fluid vars
    use mod_m1_internal
    use mod_m1_beta_eq
    use mod_eos, only: eos_rhomin, eos_rhomax, eos_yemin, eos_yemax, eos_tempmin, eos_tempmax
    use, intrinsic :: ieee_arithmetic  ! NEW: For NaN checking
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    use amrvacdef
    }
    {#IFNDEF UNIT_TESTS
      include 'amrvacdef.f' 
    }
    integer, intent(in) :: igrid,ixI^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout)    :: wprim(ixI^S,1:nw) ! wprim is not output
    double precision, intent(inout) :: wradimpl(ixI^S,1:nm1rad_eas)
    double precision, intent(in) :: qdt,qtC
    logical, intent(in) :: calc_eqRates
    !internal
    integer :: ix^D, ixO^L
    double precision :: tau_beta_min, tau_beta
    double precision :: T_eq,Ye_eq !> correct equilibrium fluid vars

    double precision :: T_old, Ye_old                    ! Old fluid values
    double precision :: beta_equil_tscale                ! Dimensionless timescale parameter
    double precision :: fac                              ! Interpolation factor

    double precision, dimension(1:3) :: fluid_Prim
    double precision, dimension(1:^NS) :: N_KSP
    double precision, dimension(1:m1_numvars_internal, ^NS) :: wrad_KSP
    double precision, dimension(1:^NC) :: vel_
    double precision, dimension(1:^NS) :: eta_nu
    double precision, dimension(1:^NC, 1:^NS) :: vel_closure
    double precision, dimension(1:^NS) :: Gamma_closure,Wlor_closure,Jrad_closure
    integer, dimension(1:^NS) :: signum
    logical :: calc_eqRates_here   !> whether or not to calculates equilib. kappas
    type(m1_metric_helper) :: metricM1_help 
    type(m1_closure_helpers) :: stateM1_help 	
    !double precision :: tset

    !tset = 0.3d0

    ixO^L=ixI^L^LSUBdixB;

    {#IFDEF M1_TAU_OPT
    calc_eqRates_here = calc_eqRates !false if not using tau_path
    }
    {#IFNDEF M1_TAU_OPT
    calc_eqRates_here = .true.
    }

    tau_beta_min = 1.0d+40
    timeCurrent = qtC

    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1_help)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1_help%fill_metric(wprim,x,ixI^L,ixO^L)
    }

    !allocate space for m1_eas
    {#IFNDEF UNIT_TESTS
    call m1_allocate_eas_ixD()
    \}
  
    !*********************************************************************      
    {do ix^D=ixOmin^D,ixOmax^D \}

      fluid_Prim(idx_rho) = wprim(ix^D,rho_) 
      fluid_Prim(idx_T) = wprim(ix^D,T_eps_) 
      fluid_Prim(idx_Ye) = wprim(ix^D,Ye_) 

      !KEN this if statements prevents eas calculations were it would be zero.
      if(fluid_Prim(idx_rho) .le. small_rho_thr) then
        {^KSP&
       wradimpl(ix^D,Q_er^KSP_) = 1.0d-30
       wradimpl(ix^D,kappa_a^KSP_) = 1.0d-30
       wradimpl(ix^D,kappa_s^KSP_) = 1.0d-30
       wradimpl(ix^D,Q_nr^KSP_) = 1.0d-30
       wradimpl(ix^D,kappa_nr^KSP_) = 1.0d-30
      \}
        cycle
      endif

      {^KSP&  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        N_KSP(^KSP) = wprim(ix^D,nrad^KSP_)/ metricM1_help%sqrtg(ix^D)
        wrad_KSP(m1_energy_,^KSP) = wprim(ix^D,erad^KSP_)/ metricM1_help%sqrtg(ix^D)
        {^C& wrad_KSP(m1_flux^C_,^KSP) = wprim(ix^D,frad^KSP^C_)/ metricM1_help%sqrtg(ix^D) \}
        {^C& vel_(^C) = wprim(ix^D,u0_+^C) \}
        !----------------------------------
	      stateM1_help%E = wrad_KSP(m1_energy_,^KSP)
        {^C& stateM1_help%F_low(^C) = wrad_KSP(m1_flux^C_,^KSP) \} 
        {^C& stateM1_help%vel(^C)   = vel_(^C) \} 
	      
        !> update closure
	      call m1_update_closure_ixD(stateM1_help,metricM1_help,ix^D,.true.,get_vel_impl=.true.)
        !> Now J, H_i are up to date, fill wrad with update values
        wrad_KSP(m1_energy_,^KSP) = stateM1_help%E
        {^C& wrad_KSP(m1_flux^C_,^KSP) = stateM1_help%F_low(^C) \}
        Gamma_closure(^KSP) = stateM1_help%Gamma
        {^C& vel_closure(^C,^KSP) = stateM1_help%vel(^C) \}
        Wlor_closure(^KSP) = stateM1_help%W
        Jrad_closure(^KSP) = stateM1_help%J
        !-----------------------------------
        {#IFDEF M1_TAU_OPT
           !> the equilibium kappas were calculated already for tau_path 
           m1_eas_ixD(Q_ems,^KSP) = wradimpl(ix^D,Q_er^KSP_) 
           m1_eas_ixD(k_a,^KSP) = wradimpl(ix^D,kappa_a^KSP_)
           m1_eas_ixD(k_s,^KSP) = wradimpl(ix^D,kappa_s^KSP_)
           m1_eas_ixD(Q_ems_n,^KSP) = wradimpl(ix^D,Q_nr^KSP_)
           m1_eas_ixD(k_n,^KSP) = wradimpl(ix^D,kappa_nr^KSP_)   
           m1_eas_ixD(tau_path,^KSP) = wradimpl(ix^D,tau_path^KSP_)            
        }
        {#IFDEF M1_RATES
          !> ----- get the rates: 
        call calculate_m1_eas_ixD(wrad_KSP(:,^KSP),N_KSP(^KSP),ix^D,Gamma_closure(^KSP), &
                vel_closure(:,^KSP),Wlor_closure(^KSP),Jrad_closure(^KSP),ixI^L,x,m1_eas_ixD(:,^KSP),^KSP,&
                fluid_Prim,eta_nu(^KSP),calc_eqRates_here,qdt)
        }
    
      \} ! end KSP !++++++++++++++++++++++++++++++++++++++++++++++++++++
       
      !------------------------------------------
!KEN changed it nearly completely
      !------------------------------------------
      {#IFDEF M1_BETA_EQ

        ! ====================================================================
        ! Calculate beta-equilibration timescale for all species
        ! ====================================================================
        tau_beta_min = huge(1.0d0)  ! Initialize to large value
        
        {^KSP& 
        ! tau_beta = 1 / sqrt(k_a * (k_a + k_s))
        tau_beta = 1.0d0 / sqrt(m1_eas_ixD(k_a,^KSP) * &
                   (m1_eas_ixD(k_a,^KSP) + m1_eas_ixD(k_s,^KSP)) + 1.0d-45)
        
        ! Get minimum across all species
        tau_beta_min = min(tau_beta_min, tau_beta)
        \}
        
        ! ====================================================================
        ! Compute dimensionless timescale parameter
        ! beta_equil_tscale = tau_beta / dt
        ! ====================================================================
        if (tau_beta_min > 0.0d0) then
            beta_equil_tscale = tau_beta_min / qdt
        else
            beta_equil_tscale = 0.0d0
        end if
        
        ! ====================================================================
        ! Only apply beta-equilibrium correction if tau_beta < dt
        ! (i.e., beta_equil_tscale < 1.0)
        ! ====================================================================
        if (beta_equil_tscale < 1.0d0 .and. N_KSP(1) > 1.d-16 .and. N_KSP(2) > 1.d-16 .and. N_KSP(3) > 1.d-16 ) then
        
            ! Save old values
            T_old  = fluid_Prim(idx_T)
            Ye_old = fluid_Prim(idx_Ye)
            
            ! Solve for beta-equilibrium
            call m1_get_beta_equilibrium(wrad_KSP, N_KSP, fluid_Prim, &
                 Jrad_closure, Gamma_closure, ix^D, T_eq, Ye_eq)
            
            ! ----------------------------------------------------------------
            ! Interpolation logic (exactly like C++)
            ! ----------------------------------------------------------------
            if (beta_equil_tscale < 0.5d0) then
                ! Fast equilibration: use full equilibrium values
                fluid_Prim(idx_T)  = T_eq
                fluid_Prim(idx_Ye) = Ye_eq
                
            else
                ! Intermediate regime: linear interpolation
                ! fac = 2 * (1 - beta_equil_tscale)
                ! Maps tscale from [0.5, 1.0] to fac from [1.0, 0.0]
                fac = 2.0d0 * (1.0d0 - beta_equil_tscale)
                fac = max(0.0d0, min(1.0d0, fac))  ! Safety clamp
                
                fluid_Prim(idx_T)  = fac * T_eq  + (1.0d0 - fac) * T_old
                fluid_Prim(idx_Ye) = fac * Ye_eq + (1.0d0 - fac) * Ye_old
            end if

            if (.not. ieee_is_finite(fluid_Prim(idx_T)) .or. &
                .not. ieee_is_finite(fluid_Prim(idx_Ye))) then
                
                write(*,*) "================================"
                write(*,*) "NaN DETECTED in beta-equilibrium at "
                write(*,*) "================================"
                write(*,*) "  T_old  =", T_old
                write(*,*) "  Ye_old =", Ye_old
                write(*,*) "  T_eq   =", T_eq
                write(*,*) "  Ye_eq  =", Ye_eq
                write(*,*) "  fac    =", fac
                write(*,*) "  beta_equil_tscale =", beta_equil_tscale
                write(*,*) "  tau_beta_min      =", tau_beta_min
                write(*,*) "  qdt               =", qdt
                write(*,*) "  T_new  =", fluid_Prim(idx_T)
                write(*,*) "  Ye_new =", fluid_Prim(idx_Ye)
                write(*,*) "  rho    =", fluid_Prim(idx_rho)
                write(*,*) "================================"
                
                ! Fallback: keep old values
                fluid_Prim(idx_T)  = T_old
                fluid_Prim(idx_Ye) = Ye_old
                
                write(*,*) "  Falling back to old values"
                write(*,*) "================================"
            end if

            
            calc_eqRates_here = .true.
            
            {^KSP&
            call calculate_m1_eas_ixD(wrad_KSP(:,^KSP), N_KSP(^KSP), ix^D, &
                    Gamma_closure(^KSP), vel_closure(:,^KSP), &
                    Wlor_closure(^KSP), Jrad_closure(^KSP), ixI^L, x, &
                    m1_eas_ixD(:,^KSP), ^KSP, fluid_Prim, eta_nu(^KSP), &
                    calc_eqRates_here, qdt)
            \}
            
        end if  ! beta_equil_tscale < 1.0
        
      } ! endif M1_BETA_EQ
      !------------------------------------------
     !> update rates vars: 
     {^KSP&
       wradimpl(ix^D,Q_er^KSP_) = m1_eas_ixD(Q_ems,^KSP) 
       wradimpl(ix^D,kappa_a^KSP_) = m1_eas_ixD(k_a,^KSP)
       wradimpl(ix^D,kappa_s^KSP_) = m1_eas_ixD(k_s,^KSP) 
       wradimpl(ix^D,Q_nr^KSP_) = m1_eas_ixD(Q_ems_n,^KSP)
       wradimpl(ix^D,kappa_nr^KSP_) = m1_eas_ixD(k_n,^KSP) 

       
        {#IFDEF DEBUG_EAS
        if(timeCurrent>2.d0 * m1_tset) then
        write(95,*)"m1_eas: Q_er",wradimpl(ix^D,Q_er^KSP_)
        write(95,*)"m1_eas: Q_nr",wradimpl(ix^D,Q_nr^KSP_)
        write(95,*)"m1_eas: k_a",wradimpl(ix^D,kappa_a^KSP_)
        write(95,*)"m1_eas: k_s",wradimpl(ix^D,kappa_s^KSP_)
        write(95,*)"m1_eas: k_nr",wradimpl(ix^D,kappa_nr^KSP_)
        end if 
        }
     \}  

  {enddo ^D&\}
  !********************************************************************* 

  call metricM1_help%destroy()

  ! deallocate m1_eas
  call m1_deallocate_eas_ixD()
  !
  end subroutine m1_add_eas

 !****************************************************************************
 
  subroutine calculate_m1_eas_ixD(wrad,N,ix^D,Gamma,vel,Wlor,Jrad,ixI^L,x,eas, speciesKSP, fluid_Prim, eta_rad,calc_eqRates_here, qdt)    
    use mod_m1_internal
    use mod_m1_fermi
    use mod_m1_constants
    use mod_m1_eas_microphysical_gray
    use mod_eos_tabulated
    use mod_eos, only: small_rho, small_temp, big_ye , eos_yemin, eos_yemax
    !............................
    use mod_Weakhub_reader !, only: logrho_IVtable,IVtemp,logTemp_IVtable,Ye_IVtable, &
      ! kappa_a_en_grey_table,kappa_s_grey_table,kappa_a_num_grey_table
    use mod_m1_eas_param
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    {#IFNDEF UNIT_TESTS
    !include "amrvacdef.f"
    }
    integer, intent(in) :: ix^D, ixI^L
    integer, intent(in) :: speciesKSP
    double precision, intent(inout) :: wrad(1:m1_numvars_internal)
    double precision, intent(inout) :: fluid_Prim(1:fluid_vars) !> fluid variables
    double precision, intent(inout) :: eas(1:m1_num_eas)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(in)    :: Gamma, Wlor, Jrad
    double precision, intent(in)    :: vel(1:^NC) 
    double precision, intent(inout) :: N      !> (nrad*Gamma)/sqrtg
    double precision, intent(inout) :: eta_rad !> degeneracy parameter of neutrinos
    logical, intent(in) :: calc_eqRates_here  !> if to calculate equilibrium opacities
    double precision, intent(in)    :: qdt    !> timestep
    ! internal
    double precision, parameter :: Qnp1 = 1.2935 ! MeV !> rest-mass energy difference n&p (Ruffert 95)
    double precision :: eas_eq(1:m1_num_eas)
    double precision :: T_fluid, rho_fluid, ye_fluid 
    double precision :: tau_opt    !> optical depth
    double precision :: T_nu       !> temperature of neutrinos
    double precision :: av_energy  !> average energy of neutrinos
    double precision :: eta_e !> degeneracy parameter of electrons
    double precision :: Black_er   !> Blackbody energy rates
    double precision :: Black_nr   !> Blackbody energy rates
    double precision :: n_prim
    double precision :: chem_pot   !> neutrino chemical potential
    double precision :: Q_free, Q_free_n !> emissivity of free streaming alla Ruffert et.al 95
    double precision :: eps,prs,ent,cs2,dedt,dpderho,dpdrhoe
    double precision :: xa,xh,xn,xp,abar1,zbar1,mu_e,mu_n,mu_p,muhat,munu
    double precision :: UNITS_rho,UNITS_chem_pot,rho_cgs
    double precision :: Qems,Qems_n,kappa_a_analyt,kappa_n_analyt
    double precision :: Qems_tmp, Qems_n_tmp
    logical :: is_known_tau_opt = .false.
    !double precision :: tset !0.1d0 !0.5d0 !1.0d-2
    logical :: eas_correction_standard = .true. !.true.
    logical :: eas_correction_Tnu = .false. !.true.
    logical :: eas_correction_tanh = .false.
    double precision :: tanh_factor

    !...........................
    integer :: i,j,k,itemp
    double precision :: LogTemp
    double precision ATmin,ATmax,AYemin,AYemax,Arhomin,Arhomax,Atstep,Arhostep,Ayestep
    integer :: Ntmax,Nrhomax,Nyemax
    integer:: iloop  
    {#IFDEF UNIT_TESTS
    !> activate gray else called in m1_startup
       call m1_eas_gray_activate()
    }

    !tset = 0.3d0!1000000000.d0 !0.3d0 !0.4d0 !30.0*qdt !20.0*qdt !15.0 * qdt

    !n_prim = N/Gamma  ! not used

    rho_fluid  = fluid_Prim(idx_rho)
    ye_fluid   = fluid_Prim(idx_ye)
    T_fluid    = fluid_Prim(idx_T)



     {#IFDEF M1_TAU_OPT
      if(calc_eqRates_here) then
         call m1_get_eas(wrad,speciesKSP,ix^D,eas_eq,fluid_Prim)
         eas_eq(:) = eas_eq/LENGTHGF
      else 
       eas_eq = eas !since wradimpl and eas had equilib. kappas form tau_path
      end if
     }
     {#IFNDEF M1_TAU_OPT
    !> ---- obtain the equilibrium gray rates from Weakhub --------
    !> associated m1_get_eas-pointer to whatever sub_m1_get_eas is:
    !> rates from table as Weakhub if:
    !> sub_m1_get_eas() = m1_get_eas_microphysical_gray_Weakhub via pointer m1_get_eas
    !> output: eas_eq, ! equilibrium gray k_a, k_s,k_n
    call m1_get_eas(wrad,speciesKSP,ix^D,eas_eq,fluid_Prim)
    !call m1_get_eas_microphysical_gray_Weakhub(wrad,n_prim,speciesKSP,ix^D,eas_eq,fluid_Prim)
    eas_eq(:) = eas_eq/LENGTHGF
    }


    {#IFNDEF UNIT_TESTS
    !> get the chemical potentials from EoS-table
    call tabulated_temp_get_all_one_grid(rho_fluid,T_fluid,ye_fluid,eps,prs=prs,ent=ent,&
         cs2=cs2,dedt=dedt,dpderho=dpderho,dpdrhoe=dpdrhoe,xa=xa,xh=xh,xn=xn,xp=xp,&
         abar=abar1,zbar=zbar1,mu_e=mu_e,mu_n=mu_n,mu_p=mu_p,muhat=muhat,munu=munu)
    }


    UNITS_rho = INVRHOGF 
    UNITS_chem_pot = 1.0d0

    if(speciesKSP .eq. m1_i_nue) then
      chem_pot = mu_e + mu_p - mu_n - Qnp1
    else if(speciesKSP .eq. m1_i_nuebar) then
      chem_pot = -1.0d0 * (mu_e + mu_p - mu_n - Qnp1)
    else if(speciesKSP .eq. m1_i_nux) then
      chem_pot = 0.0d0
    else if(speciesKSP .eq. m1_i_mu) then
      chem_pot = 0.0d0
    else if(speciesKSP .eq. m1_i_mubar) then
      chem_pot = 0.0d0
    end if 

     eta_rad = chem_pot/T_fluid  !> degeneracy parameter
     eta_e = mu_e / T_fluid

     is_known_tau_opt = .True.
    if(is_known_tau_opt) then
      rho_cgs = rho_fluid * UNITS_rho 
      !> get optical depth from empirical formula
      tau_opt = exp(DLOG10(10.d0)*( 0.96d0 * ( DLOG10(rho_cgs) / DLOG10(10.0d0) - 11.7d0)))
      {#IFDEF M1_TAU_OPT_USE
       tau_opt = eas(tau_path)
      }
       if((speciesKSP .eq. 1) .or. (speciesKSP .eq. 2) ) then
         eta_rad = eta_rad * (1.0d0 - exp(-tau_opt))
       end if 
    end if 

    !eas(Q_ems) = 0.0d0
    !eas(Q_ems_n) = 0.0d0

    Qems = 0.0d0
    Qems_n = 0.0d0
    Qems_tmp = 0.0d0
    Qems_n_tmp = 0.0d0

    !> calculate T_eff neutrino
    !KEN had to add T_fluid
    call get_neutrino_temp_ixD(wrad,N,ix^D,Gamma,vel,Wlor,Jrad,eta_rad, speciesKSP, fluid_Prim, av_energy, T_nu, T_fluid, timeCurrent, x(ix^D,:))

    !> get Blackbody function
    call blackbody_ixD(wrad,ix^D,speciesKSP,fluid_Prim,T_fluid,eta_rad,Black_er,Black_nr)

    if((speciesKSP .eq. m1_i_nux) .or. (speciesKSP .eq. m1_i_mu) .or. (speciesKSP .eq. m1_i_mubar)) then
      call m1_eas_analytic_brems(speciesKSP,eta_rad,ye_fluid,T_fluid,rho_fluid,av_energy,T_nu,eas,Qems,Qems_n)
      Qems_tmp = Qems_tmp + Qems
      Qems_n_tmp = Qems_n_tmp + Qems_n

      call m1_eas_analytic_pair(speciesKSP,eta_rad,eta_e,ye_fluid,T_fluid,rho_fluid,eas,Qems,Qems_n)
      Qems_tmp = Qems_tmp + Qems
      Qems_n_tmp = Qems_n_tmp + Qems_n

      call m1_eas_analytic_plasmon(speciesKSP,eta_rad,eta_e,ye_fluid,T_fluid,rho_fluid,eas,Qems,Qems_n)
      Qems_tmp = Qems_tmp + Qems
      Qems_n_tmp = Qems_n_tmp + Qems_n

      Qems_tmp = Qems_tmp * MEV_TO_ERG_KEN * RHOGF * EPSGF / TIMEGF  !* LENGTHGF ! LENGTHGF only beacause of blackbody
      Qems_n_tmp = Qems_n_tmp * MNUC_CGS * RHOGF / TIMEGF !* LENGTHGF ! LENGTHGF only beacause of blackbody

      eas_eq(k_a) = eas_eq(k_a) + Qems_tmp / (Black_er + 1.0d-40) ! Teilen ist immer tricky. Deswegen liebe safe
      eas_eq(k_n) = eas_eq(k_n) + Qems_n_tmp / (Black_nr + 1.0d-40)
    end if 

    eas(Q_ems) = Black_er * eas_eq(k_a) 
    eas(Q_ems_n) = Black_nr * eas_eq(k_n)

    eas(k_a) = eas_eq(k_a) 
    eas(k_s) = eas_eq(k_s) 
    eas(k_n) = eas_eq(k_n) 
        !----------------------------------------
    !> correct kappa_a,s,n
    if((timeCurrent .ge. m1_tset) .and. (av_energy .ge. m1_E_atmo*100) ) then
            if(eas_correction_standard) then
              !
                if((speciesKSP .eq. m1_i_nue) .or. (speciesKSP .eq. m1_i_nuebar) &
                  .or. (speciesKSP == m1_i_mu) .or. (speciesKSP == m1_i_mubar)) then
                   eas(k_a) = eas_eq(k_a) * max(1.0d0, T_nu/T_fluid)**2
                   !TEST:   
                   eas(k_s) = eas_eq(k_s) * max(1.0d0, T_nu/T_fluid)**2
                   eas(k_n) = eas_eq(k_n) * max(1.0d0, T_nu/T_fluid)**2
                   eas(Q_ems) = eas(Q_ems) * max(1.0d0, T_nu/T_fluid)**2
                   eas(Q_ems_n) = eas(Q_ems_n) * max(1.0d0, T_nu/T_fluid)**2
                else if(speciesKSP .eq. m1_i_nux) then
                   eas(k_s) = eas_eq(k_s) * max(1.0d0, T_nu/T_fluid)**2      
                end if
            else if(eas_correction_Tnu) then
               !
                if(timeCurrent > 2.0 * m1_tset) then
                   if((speciesKSP .eq. m1_i_nue) .or. (speciesKSP .eq.m1_i_nuebar) &
                  .or. (speciesKSP == m1_i_mu) .or. (speciesKSP == m1_i_mubar)) then                
                      eas(k_a) = eas_eq(k_a) * max(1.0d0, T_nu/T_fluid)**2
                      eas(k_s) = eas_eq(k_s) * max(1.0d0, T_nu/T_fluid)**2
                      eas(k_n) = eas_eq(k_n) * max(1.0d0, T_nu/T_fluid)**2
                   else if(speciesKSP .eq. m1_i_nux) then
                      eas(k_s) = eas_eq(k_s) * max(1.0d0, T_nu/T_fluid)**2    
                   end if
                else    
                   if((speciesKSP .eq. m1_i_nue) .or. (speciesKSP .eq.m1_i_nuebar) &
                  .or. (speciesKSP == m1_i_mu) .or. (speciesKSP == m1_i_mubar)) then
                      eas(k_a) = eas_eq(k_a) * max(1.0d0, T_nu/T_fluid)
                      eas(k_s) = eas_eq(k_s) * max(1.0d0, T_nu/T_fluid)
                      eas(k_n) = eas_eq(k_n) * max(1.0d0, T_nu/T_fluid)
                   else if(speciesKSP .eq. m1_i_nux) then
                      eas(k_s) = eas_eq(k_s) * max(1.0d0, T_nu/T_fluid)       
                   end if
                end if 
            else if(eas_correction_tanh) then
              !
                tanh_factor = (1.0 + tanh((T_nu/T_fluid - 1.0)*10.0d0))/2.0d0
                if((speciesKSP .eq. m1_i_nue) .or. (speciesKSP .eq.m1_i_nuebar)&
                  .or. (speciesKSP == m1_i_mu) .or. (speciesKSP == m1_i_mubar)) then
                   eas(k_a) = eas_eq(k_a) * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
                   eas(k_s) = eas_eq(k_s) * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
                   eas(k_n) = eas_eq(k_n) * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
                   eas(Q_ems) = eas(Q_ems) * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
                   eas(Q_ems_n) = eas(Q_ems_n) * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)
                else if(speciesKSP .eq. m1_i_nux) then
                   eas(k_s) = eas_eq(k_s) * (1.0d0 + tanh_factor*(T_nu/T_fluid - 1.0d0)**2)   
                   eas(k_a) = eas_eq(k_a)
                   eas(k_n) = eas_eq(k_n)    
                end if
            end if 
    end if 
    
    ! ======================================================================
    eas(Q_ems)   = max(eas(Q_ems),   1.0d-30)
    eas(Q_ems_n) = max(eas(Q_ems_n), 1.0d-30)
    eas(k_a)     = max(eas(k_a),     1.0d-30)
    eas(k_s)     = max(eas(k_s),     1.0d-30)
    eas(k_n)     = max(eas(k_n),     1.0d-30)
    
  end subroutine calculate_m1_eas_ixD

  !/**************************************/
  
  !KEN did this like in FIL
subroutine m1_eas_analytic_brems(speciesKSP,eta_nu,ye_fluid,T_fluid,rho_fluid,av_energy,T_nu,eas,Qems,Qems_n)
  use mod_m1_internal
  use mod_m1_eas_param
  use mod_m1_constants
  include "amrvacdef.f"
  
  
  integer, intent(in) :: speciesKSP
  double precision, intent(in) :: eta_nu
  double precision, intent(in) :: ye_fluid,T_fluid,rho_fluid
  double precision, intent(in) :: av_energy,T_nu
  double precision, intent(inout) :: eas(1:m1_num_eas)
  double precision, intent(inout) :: Qems, Qems_n
  
  ! Internal variables
  double precision :: Y_n, Y_p
  double precision :: R_brems, Q_brems
  double precision :: rho_cgs, T_mev

  ! Calculate neutron and proton fractions
  !KEN TODO when adding muons
  Y_p = ye_fluid
  Y_n = 1.0d0 - Y_p

  Qems = 0.0d0
  Qems_n = 0.0d0
  
  ! Convert to CGS units
  rho_cgs = rho_fluid * INVRHOGF 
  T_mev = T_fluid
  
  ! Calculate bremsstrahlung emission rates
  R_brems = 0.231d0 * (2.0778d2 * ERG_TO_MEV) * 0.5d0 * &
            (Y_n**2 + Y_p**2 + 28.0d0/3.0d0 * Y_n * Y_p) * &
            rho_cgs**2 * T_fluid**4.5d0
  
  Q_brems = R_brems * T_fluid / 0.231d0 * 0.504d0
  
  ! Apply to heavy lepton neutrinos only
  if (Q_brems > 0.0d0 .and. Q_brems < 1.0d99) then
    
    if (m1_use_muons) then
      ! 5-species case: nue, nuebar, numu, numubar, nux(tau+antitau)
      if (speciesKSP == m1_i_nux) then
        ! NUX = 2 heavy species (tau + anti-tau)
        Qems = 2.0d0 * Q_brems
        Qems_n = 2.0d0 * R_brems
      else if (speciesKSP == m1_i_mu .or. speciesKSP == m1_i_mubar) then
        ! Muon neutrino or anti-muon neutrino
        Qems = Q_brems
        Qems_n = R_brems
      else
        Qems = 0.0d0
        Qems_n = 0.0d0
      end if
      
    else
      ! 3-species case: nue, nuebar, nux(mu+antimu+tau+antitau)
      if (speciesKSP == m1_i_nux) then
        ! NUX = 4 heavy species
        Qems = 4.0d0 * Q_brems
        Qems_n = 4.0d0 * R_brems
      else
        Qems = 0.0d0
        Qems_n = 0.0d0
      end if
    end if
    
  else
    Qems = 0.0d0
    Qems_n = 0.0d0
  end if

end subroutine m1_eas_analytic_brems

!KEN pair process emission like in FIL
subroutine m1_eas_analytic_pair(speciesKSP,eta_nu,eta_e,ye_fluid,T_fluid,rho_fluid,eas,Qems,Qems_n)
  use mod_m1_internal
  use mod_m1_eas_param
  use mod_m1_constants
  use mod_m1_fermi
  include "amrvacdef.f"
  
  integer, intent(in) :: speciesKSP
  double precision, intent(in) :: eta_nu, eta_e
  double precision, intent(in) :: ye_fluid, T_fluid, rho_fluid
  double precision, intent(inout) :: eas(1:m1_num_eas)
  double precision, intent(inout) :: Qems, Qems_n
  
  ! Internal variables
  double precision :: block_factor_nue, block_factor_nuebar, block_factor_nux
  double precision :: block_factor_numu, block_factor_numubar
  double precision :: eps_const, eps_m, eps_p, eps_fraction
  double precision :: pair_const, R_pair, R_pair_x
  double precision :: R_pair_numu_numubar
  double precision :: T_mev
  double precision :: FD4ep, FD4em
  double precision :: FD3ep, FD3em
  double precision :: FDRp, FDRm
  
  ! Constants from C++ code
  double precision, parameter :: me_mev = 0.510998910d0  ! MeV

  Qems = 0.0d0
  Qems_n = 0.0d0
  
  ! Early exit for NUE and NUEBAR (not used - overproducing)
  if((speciesKSP .eq. m1_i_nue) .or. (speciesKSP .eq. m1_i_nuebar)) then
    return
  endif
  
  T_mev = T_fluid

  FD4ep = fermi_dirac(eta_e, 4)
  FD4em = fermi_dirac(-eta_e, 4)
  FD3ep = fermi_dirac(eta_e, 3)
  FD3em = fermi_dirac(-eta_e, 3)

  if (FD3ep .gt. 1.0d-15) then
    FDRp = FD4ep/FD3ep
  else 
    FDRp = 4
  endif
  if (FD3em .gt. 1.0d-15) then
    FDRm = FD4em/FD3em
  else 
    FDRm = 4
  endif
  
  ! Calculate blocking factor for heavy leptons (nux)
  block_factor_nux = 1.0d0 + exp(0.0d0 - 0.5d0 * &  ! eta_nux = 0
                     (FDRp + FDRm))
  
  ! eps_const = 8π / (ℏc)³
  eps_const = 8.0d0 * PI / HC_MEVCM**3
  
  ! Ruffert et al. (B5)
  eps_m = eps_const * T_mev**4 * FD3ep
  eps_p = eps_const * T_mev**4 * FD3em
  
  ! Ruffert et al. (B16) - average energy fraction
  eps_fraction = 0.5d0 * T_mev * &
                 (FDRp + FDRm)
  
  ! Pair constant - Ruffert (B8)
  pair_const = sigma_0 * CLITE_CM / me_mev**2 * eps_m * eps_p
  
  ! ================================================================
  ! NUMU, NUMU_BAR, NUX
  ! ================================================================
  if (m1_use_muons) then
    ! 5-species case
    block_factor_numu = 1.0d0 + exp(eta_nu - 0.5d0 * &  ! eta_nux = 0
                     (FDRp + FDRm))
    block_factor_numubar = 1.0d0 + exp(-eta_nu - 0.5d0 * &  ! eta_nux = 0
                     (FDRp + FDRm))
    
    ! Ruffert et al. (B10) - for numu and numubar
    R_pair_numu_numubar = (1.0d0 / 36.0d0) * pair_const * &
                          ((Cv - Ca)**2 + (Cv + Ca - 2.0d0)**2) / &
                          (block_factor_numu * block_factor_numubar)
    
    if (speciesKSP == m1_i_mu .or. speciesKSP == m1_i_mubar) then
      if (R_pair_numu_numubar > 0.0d0 .and. R_pair_numu_numubar < 1.0d99) then
        Qems = R_pair_numu_numubar * eps_fraction
        Qems_n = R_pair_numu_numubar
      end if
      return
    end if
    
    ! NUX (tau + antitau only, factor 2)
    R_pair_x = (1.0d0 / 18.0d0) * pair_const * &
               ((Cv - Ca)**2 + (Cv + Ca - 2.0d0)**2) / &
               block_factor_nux**2
    
  else
    ! 3-species case: NUX includes mu, antimu, tau, antitau (factor 4)
    R_pair_x = (1.0d0 / 9.0d0) * pair_const * &
               ((Cv - Ca)**2 + (Cv + Ca - 2.0d0)**2) / &
               block_factor_nux**2
  end if
  
  ! Apply to NUX
  if (speciesKSP == m1_i_nux) then
    if (R_pair_x > 0.0d0 .and. R_pair_x < 1.0d99) then
      Qems = R_pair_x * eps_fraction
      Qems_n = R_pair_x
    end if
  end if

end subroutine m1_eas_analytic_pair

!KEN plasmon decay emission like in FIL
subroutine m1_eas_analytic_plasmon(speciesKSP,eta_nu,eta_e,ye_fluid,T_fluid,rho_fluid,eas,Qems,Qems_n)
  use mod_m1_internal
  use mod_m1_eas_param
  use mod_m1_constants
  use mod_m1_fermi
  include "amrvacdef.f"
  
  integer, intent(in) :: speciesKSP
  double precision, intent(in) :: eta_nu, eta_e
  double precision, intent(in) :: ye_fluid, T_fluid, rho_fluid
  double precision, intent(inout) :: eas(1:m1_num_eas)
  double precision, intent(inout) :: Qems, Qems_n
  
  ! Internal variables
  double precision :: gamma, gamma_const
  double precision :: block_factor_nue, block_factor_nuebar, block_factor_nux
  double precision :: block_factor_numu, block_factor_numubar
  double precision :: R_gamma, R_gamma_numu, R_gamma_numu_bar, R_gamma_x
  double precision :: Q_gamma
  double precision :: T_mev

  double precision :: FD4ep, FD4em
  double precision :: FD3ep, FD3em
  double precision :: FDRp, FDRm


  ! Constants from C++ code
  double precision, parameter :: me_mev = 0.510998910d0  ! MeV
  double precision :: fsc = 0.0072973525205055605d0 ! Fine structure constant
  
  Qems = 0.0d0
  Qems_n = 0.0d0
  
  ! NUE and NUE_BAR - early return (overproducing)
  if((speciesKSP .eq. m1_i_nue) .or. (speciesKSP .eq. m1_i_nuebar)) then
    return
  endif
  
  T_mev = T_fluid

  FD4ep = fermi_dirac(eta_e, 4)
  FD4em = fermi_dirac(-eta_e, 4)
  FD3ep = fermi_dirac(eta_e, 3)
  FD3em = fermi_dirac(-eta_e, 3)

  if (FD3ep .gt. 1.0d-15) then
    FDRp = FD4ep/FD3ep
  else 
    FDRp = 4
  endif
  if (FD3em .gt. 1.0d-15) then
    FDRm = FD4em/FD3em
  else 
    FDRm = 4
  endif
  
  ! Ruffert et al. see comment after (B12)
  ! gamma = gamma_0 * sqrt(1/3 * (pi^2 + 3*eta_e^2))
  gamma = GAMMA_0 * sqrt((1.0d0 / 3.0d0) * (PI**2 + 3.0d0 * eta_e**2))
  
  ! Ruffert et al. (B11)
  ! gamma_const = pi^3 * sigma_0 * c / (me^2 * 3 * fsc * (hc)^6) * 
  !               gamma^6 * exp(-gamma) * (1+gamma) * T^8
  gamma_const = PI**3 * SIGMA_0 * CLITE_CM / &
                (ME_MEV**2 * 3.0d0 * FSC * HC_MEVCM) * &
                gamma**6 * exp(-gamma) * (1.0d0 + gamma) * T_mev**8
  
  gamma_const = gamma_const / HC_MEVCM**1
  gamma_const = gamma_const / HC_MEVCM**1
  gamma_const = gamma_const / HC_MEVCM**1
  gamma_const = gamma_const / HC_MEVCM**1
  gamma_const = gamma_const / HC_MEVCM**1

  
  ! Calculate blocking factor for heavy leptons (nux)
  ! block_factor[i] = 1 + exp(eta[i] - (1 + 0.5*gamma^2/(1+gamma)))
  block_factor_nux = 1.0d0 + exp(0.0d0 - &  ! eta_nux = 0
                     (1.0d0 + 0.5d0 * gamma**2 / (1.0d0 + gamma)))
  
  ! Q_gamma - average energy factor (Ruffert et al. B11)
  Q_gamma = T_mev * 0.5d0 * (2.0d0 + gamma**2 / (1.0d0 + gamma))
  
  ! ================================================================
  ! NUMU, NUMU_BAR, NUX
  ! Ruffert et al. (B12)
  ! ================================================================
  if (m1_use_muons) then
    ! 5-species case
    ! KEN: Use pair-process style blocking factors for muons
    block_factor_numu = 1.0d0 + exp(eta_nu - 0.5d0 * &
                        (FDRp + FDRm))
    
    block_factor_numubar = 1.0d0 + exp(-eta_nu - 0.5d0 * &
                           (FDRp + FDRm))
    
    ! R_gamma for numu and numubar
    R_gamma_numu = (Cv - 1.0d0)**2 * gamma_const / &
                   (block_factor_numu * block_factor_numubar)
    
    R_gamma_numu_bar = (Cv - 1.0d0)**2 * gamma_const / &
                       (block_factor_numu * block_factor_numubar)
    
    if (speciesKSP == m1_i_mu) then
      if (R_gamma_numu > 0.0d0 .and. R_gamma_numu < 1.0d99) then
        Qems = Q_gamma * R_gamma_numu
        Qems_n = R_gamma_numu
      end if
      return
    end if
    
    if (speciesKSP == m1_i_mubar) then
      if (R_gamma_numu_bar > 0.0d0 .and. R_gamma_numu_bar < 1.0d99) then
        Qems = Q_gamma * R_gamma_numu_bar
        Qems_n = R_gamma_numu_bar
      end if
      return
    end if
    
    ! NUX (tau + antitau only, factor 2)
    R_gamma_x = 2.0d0 * (Cv - 1.0d0)**2 * gamma_const / &
                block_factor_nux**2
    
  else
    ! 3-species case: NUX includes mu, antimu, tau, antitau (factor 4)
    R_gamma_x = 4.0d0 * (Cv - 1.0d0)**2 * gamma_const / &
                block_factor_nux**2
  end if
  
  ! Apply to NUX
  if (speciesKSP == m1_i_nux) then
    if (R_gamma_x > 0.0d0 .and. R_gamma_x < 1.0d99) then
      Qems = Q_gamma * R_gamma_x
      Qems_n = R_gamma_x
    end if
  end if

end subroutine m1_eas_analytic_plasmon

  subroutine m1_eas_analytic_kappa_AN(wrad,speciesKSP,eta_nu,mu_e,mu_n,mu_p,ye_fluid,T_fluid,rho_fluid,av_energy,T_nu,eas,kappa_a_analyt,kappa_n_analyt)
    use mod_m1_internal
    use mod_m1_fermi
    use mod_m1_constants
    integer, intent(in) :: speciesKSP
    double precision, intent(in) :: eta_nu,mu_e,mu_n,mu_p
    double precision, intent(in) :: ye_fluid,T_fluid,rho_fluid
    double precision, intent(in) :: av_energy,T_nu
    double precision, intent(in) :: wrad(1:m1_numvars_internal)
    double precision, intent(inout) :: eas(1:m1_num_eas)
    double precision, intent(inout) :: kappa_a_analyt, kappa_n_analyt
    ! internal:
    double precision, parameter :: m_el = 0.510998910d0    ! mass of electron in MeV
    double precision, parameter :: alpha_const = 1.23d0    ! no dim
    double precision, parameter :: sigma0 = 1.76d-44       ! in units of cm^2
    double precision, parameter :: AA = 6.0221367d23       !> Avogadro constant
    double precision, parameter :: clight = 2.99792458d+10 !> speed of light
    double precision, parameter :: hc_const = 1.23984172d-10 ! MeV cm
    double precision, parameter :: Pi_const = 3.14159265358979 !> Pi

    double precision :: xi_pn,xi_np,Y_n,Y_p,Ap
    double precision :: eta_n, eta_p, eta_e, Ynp, Ypn
    double precision :: fermi_fac, fermi_block
    double precision :: energy_av    

    Y_p = ye_fluid
    Y_n = Y_p - 1.0d0

    ! fuagcities
    eta_p = mu_p/T_fluid
    eta_n = mu_n/T_fluid
    eta_e = mu_e/T_fluid

    ! baryion number density
    ! nb = rho/mb = Na*rho
    !Ap = rho_fluid / MNUC_CU
    Ap = rho_fluid * INVRHOGF / MNUC_CGS

    xi_np = Ap * (Y_p - Y_n)/(exp(eta_p - eta_n) - 1.0d0)
    xi_pn = Ap * (Y_n - Y_p)/(exp(eta_n - eta_p) - 1.0d0)

    if(speciesKSP .eq. m1_i_nue) then
      ! average energy positrons
      energy_av = T_fluid * fermi_dirac(eta_nu,5)/fermi_dirac(eta_nu,4)
      ! fermi blocking factor:
      fermi_block = 1.0d0 /(1.0d0 + exp(-1.0d0 *(energy_av/T_fluid - eta_e)))
      ! ------ kappa_n nue -----------
      fermi_fac = fermi_dirac(eta_nu,4)/fermi_dirac(eta_nu,2)
      ! kappa_n
      !kappa_n_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_np * sigma0 * (T_fluid/m_el/clight**2.0)**2.0d0 * fermi_fac * fermi_block
      kappa_n_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_np * sigma0 * (T_fluid/m_el)**2.0d0 * fermi_fac * fermi_block
      ! ------ kappa_a nue_bar -----------
      fermi_fac = fermi_dirac(eta_nu,5)/fermi_dirac(eta_nu,3)
      ! kappa_a 
      !kappa_a_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_np * sigma0 * (T_fluid/m_el/clight**2.0)**2.0d0 * fermi_fac * fermi_block
      kappa_a_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_np * sigma0 * (T_fluid/m_el )**2.0d0 * fermi_fac * fermi_block

    else if(speciesKSP .eq. m1_i_nuebar) then
      ! average energy positrons
      energy_av = T_fluid * fermi_dirac(eta_nu,5)/fermi_dirac(eta_nu,4)
      ! fermi blocking factor:
      fermi_block = 1.0d0 /(1.0d0 + exp(-1.0d0 *(energy_av/T_fluid + eta_e)))
      ! ------ kappa_n -----------
      fermi_fac = fermi_dirac(eta_nu,4)/fermi_dirac(eta_nu,2)
      ! kappa_a
      !kappa_n_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_pn * sigma0 * (T_fluid/m_el/clight**2.0)**2.0d0 * fermi_fac * fermi_block
      kappa_n_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_pn * sigma0 * (T_fluid/m_el)**2.0d0 * fermi_fac * fermi_block
      ! ------ kappa_a -----------
      fermi_fac = fermi_dirac(eta_nu,5)/fermi_dirac(eta_nu,3)
      ! kappa_a
      !kappa_a_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_pn * sigma0 * (T_fluid/m_el/clight**2.0d0)**2.0d0 * fermi_fac * fermi_block
      kappa_a_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_pn * sigma0 * (T_fluid/m_el)**2.0d0 * fermi_fac * fermi_block

      
      !eas(k_n)  = kappa_n_analyt
      !eas(k_a)  = kappa_a_analyt
    end if 

    !eas(k_n)  = kappa_n_analyt
    !eas(k_a)  = kappa_a_analyt

  end subroutine m1_eas_analytic_kappa_AN


  !----------------------------------
  function f_tau(tau_opt)
     double precision :: f_tau, tau_opt
     if(tau_opt .gt. 2.0d0) then
      f_tau = 1.0d0
     else if(tau_opt .lt. 2.0d0/3.0d0) then
      f_tau = 0.0d0
     else 
      f_tau = (tau_opt - 2.0d0/3.0d0)/(4.0d0/3.0d0)
     end if 
    end function f_tau
  !----------------------------------
  ! this routine is not used and undone!!
  subroutine Q_Ruffert(wrad,speciesKSP,eta_nu,mu_e,mu_n,mu_p,ye_fluid,T_fluid,rho_fluid,av_energy,T_nu,eas,Q_free,Q_free_n)
    use mod_m1_internal
    use mod_m1_fermi
    integer, intent(in) :: speciesKSP
    double precision, intent(in) :: eta_nu,mu_e,mu_n,mu_p
    double precision, intent(in) :: ye_fluid,T_fluid,rho_fluid
    double precision, intent(in) :: av_energy,T_nu
    double precision, intent(in) :: wrad(1:m1_numvars_internal)
    double precision, intent(inout) :: eas(1:m1_num_eas)
    double precision, intent(out) :: Q_free, Q_free_n
    ! internal:
    double precision, parameter :: m_el = 0.510998910d0    ! mass of electron in MeV
    double precision, parameter :: alpha_const = 1.23d0    ! no dim
    double precision, parameter :: sigma0 = 1.76d-44       ! in units of cm^2
    double precision, parameter :: AA = 6.0221367d23       !> Avogadro constant
    double precision, parameter :: clight = 2.99792458d+10 !> speed of light
    double precision, parameter :: hc_const = 1.23984172d-10 ! MeV cm
    double precision, parameter :: Pi_const = 3.14159265358979 !> Pi
    ! TODO:
    double precision, parameter :: UnitsR = 1.0d0
    double precision, parameter :: UnitsQ = 1.0d0
    double precision, parameter :: UnitsQ_n = 1.0d0

    double precision :: eta_n, eta_p, eta_e, Ynp, Ypn
    double precision :: fermi_fac, R_free

    eta_p = mu_p/T_fluid
    eta_n = mu_n/T_fluid
    eta_e = mu_e/T_fluid

    Ynp = (2.0d0 *ye_fluid -1.0d0)/(exp(eta_p - eta_n)-1.0d0)
    Ypn = exp(eta_p - eta_n) * Ynp

    if(speciesKSP .eq. m1_i_nue) then
      !> electron antineutrino : beta process
      fermi_fac = 1.0d0 / (1.0d0 +  exp(-1.0d0*(fermi_dirac(eta_e,5)/fermi_dirac(eta_e,4)- eta_nu)))
      
      R_free = (1.0d0 + alpha_const**2)/8.0d0 * sigma0 * clight/(m_el * clight**2)**2 
      R_free = Q_free * AA * rho_fluid * Ypn * fermi_fac * UnitsR !* energy_electron 
     
      Q_free = R_free * 8.0d0 * Pi_const/(hc_const)**3 * T_fluid**6 * fermi_dirac(eta_e,5) * UnitsQ
      Q_free_n = R_free * 8.0d0 * Pi_const/(hc_const)**3 * T_fluid**5 * fermi_dirac(eta_e,4) * UnitsQ_n
    
    else if(speciesKSP .eq. m1_i_nuebar) then 
      !> electron antineutrino : beta process
      fermi_fac = 1.0d0 / (1.0d0 +  exp(-1.0d0*(fermi_dirac(-eta_e,5)/fermi_dirac(-eta_e,4)- eta_nu)))
      
      R_free = (1.0d0 + alpha_const**2)/8.0d0 * sigma0 * clight/(m_el * clight**2)**2 
      R_free = Q_free * AA * rho_fluid * Ypn * fermi_fac * UnitsR !* energy_positron
      
      Q_free = R_free * 8.0d0 * Pi_const/(hc_const)**3 * T_fluid**6 * fermi_dirac(-eta_e,5) * UnitsQ
      Q_free_n = R_free * 8.0d0 * Pi_const/(hc_const)**3 * T_fluid**5 * fermi_dirac(-eta_e,4) * UnitsQ_n
    end if 
  end subroutine Q_Ruffert
  !------------------------------------------------
    !***************************************************
  subroutine m1_eas_analytic(wrad,speciesKSP,eta_nu,mu_e,mu_n,mu_p,ye_fluid,T_fluid,rho_fluid,av_energy,T_nu,eas)
    use mod_m1_internal
    use mod_m1_fermi
    use mod_m1_constants
    integer, intent(in) :: speciesKSP
    double precision, intent(in) :: eta_nu,mu_e,mu_n,mu_p
    double precision, intent(in) :: ye_fluid,T_fluid,rho_fluid
    double precision, intent(in) :: av_energy,T_nu
    double precision, intent(in) :: wrad(1:m1_numvars_internal)
    double precision, intent(inout) :: eas(1:m1_num_eas)
    ! internal:
    double precision, parameter :: m_el = 0.510998910d0    ! mass of electron in MeV
    double precision, parameter :: alpha_const = 1.23d0    ! no dim
    double precision, parameter :: sigma0 = 1.76d-44       ! in units of cm^2
    double precision, parameter :: AA = 6.0221367d23       !> Avogadro constant
    double precision, parameter :: clight = 2.99792458d+10 !> speed of light
    double precision, parameter :: hc_const = 1.23984172d-10 ! MeV cm
    double precision, parameter :: Pi_const = 3.14159265358979 !> Pi
    ! TODO:
    double precision, parameter :: UnitsR = 1.0d0
    double precision, parameter :: UnitsQ = 1.0d0
    double precision, parameter :: UnitsQ_n = 1.0d0

    double precision :: eta_n, eta_p, eta_e, Ynp, Ypn
    double precision :: xi_pn,xi_np,Yn,Yp,Ap
    double precision :: fermi_fac, fermi_block, R_free
    double precision :: Q_free, Q_free_n
    double precision :: kappa_a_analyt, kappa_n_analyt
    double precision :: energy_av    
    double precision :: Y_n,Y_p

    Y_p = ye_fluid
    Y_n = 1.0d0 - Y_p

    eta_p = mu_p/T_fluid
    eta_n = mu_n/T_fluid
    eta_e = mu_e/T_fluid

    Ynp = (2.0d0 *ye_fluid -1.0d0)/(exp(eta_p - eta_n)-1.0d0)
    Ypn = exp(eta_p - eta_n) * Ynp
    
    ! baryion number density
    ! nb = rho/mb = Na*rho
    !Ap = rho_fluid / MNUC_CU
    Ap = rho_fluid * INVRHOGF / MNUC_CGS

    xi_np = Ap * (Y_p - Y_n)/(exp(eta_p - eta_n) - 1.0d0)
    xi_pn = Ap * (Y_n - Y_p)/(exp(eta_n - eta_p) - 1.0d0)


    if(speciesKSP .eq. m1_i_nue) then
      !> electron antineutrino : beta process
      fermi_fac = 1.0d0 / (1.0d0 +  exp(-1.0d0*(fermi_dirac(eta_e,5)/fermi_dirac(eta_e,4)- eta_nu)))
      
      R_free = (1.0d0 + alpha_const**2)/8.0d0 * sigma0 * clight/(m_el * clight**2)**2 
      R_free = Q_free * AA * rho_fluid * Ypn * fermi_fac * UnitsR !* energy_electron 
     
      Q_free = R_free * 8.0d0 * Pi_const/(hc_const)**3 * T_fluid**6 * fermi_dirac(eta_e,5) * UnitsQ
      Q_free_n = R_free * 8.0d0 * Pi_const/(hc_const)**3 * T_fluid**5 * fermi_dirac(eta_e,4) * UnitsQ_n
    
    else if(speciesKSP .eq. m1_i_nuebar) then 
      !> electron antineutrino : beta process
      fermi_fac = 1.0d0 / (1.0d0 +  exp(-1.0d0*(fermi_dirac(-eta_e,5)/fermi_dirac(-eta_e,4)- eta_nu)))
      
      R_free = (1.0d0 + alpha_const**2)/8.0d0 * sigma0 * clight/(m_el * clight**2)**2 
      R_free = Q_free * AA * rho_fluid * Ypn * fermi_fac * UnitsR !* energy_positron
      
      Q_free = R_free * 8.0d0 * Pi_const/(hc_const)**3 * T_fluid**6 * fermi_dirac(-eta_e,5) * UnitsQ
      Q_free_n = R_free * 8.0d0 * Pi_const/(hc_const)**3 * T_fluid**5 * fermi_dirac(-eta_e,4) * UnitsQ_n
    end if 

    if(speciesKSP .eq. m1_i_nue) then
      ! average energy positrons
      energy_av = T_fluid * fermi_dirac(eta_nu,5)/fermi_dirac(eta_nu,4)
      ! fermi blocking factor:
      fermi_block = 1.0d0 /(1.0d0 + exp(-1.0d0 *(energy_av/T_fluid - eta_e)))
      ! ------ kappa_n nue -----------
      fermi_fac = fermi_dirac(eta_nu,4)/fermi_dirac(eta_nu,2)
      ! kappa_n
      kappa_n_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_np * sigma0 * (T_fluid/m_el)**2.0d0 * fermi_fac * fermi_block
      ! ------ kappa_a nue_bar -----------
      fermi_fac = fermi_dirac(eta_nu,5)/fermi_dirac(eta_nu,3)
      ! kappa_a 
      kappa_a_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_np * sigma0 * (T_fluid/m_el)**2.0d0 * fermi_fac * fermi_block

    else if(speciesKSP .eq. m1_i_nuebar) then
      ! average energy positrons
      energy_av = T_fluid * fermi_dirac(eta_nu,5)/fermi_dirac(eta_nu,4)
      ! fermi blocking factor:
      fermi_block = 1.0d0 /(1.0d0 + exp(-1.0d0 *(energy_av/T_fluid + eta_e)))
      ! ------ kappa_n -----------
      fermi_fac = fermi_dirac(eta_nu,4)/fermi_dirac(eta_nu,2)
      ! kappa_a
      kappa_n_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_pn * sigma0 * (T_fluid/m_el)**2.0d0 * fermi_fac * fermi_block
      ! ------ kappa_a -----------
      fermi_fac = fermi_dirac(eta_nu,5)/fermi_dirac(eta_nu,3)
      ! kappa_a
      kappa_a_analyt = (1.0d0 + 3.0d0*ALPHA**2.0)/4.0d0 * xi_pn * sigma0 * (T_fluid/m_el)**2.0d0 * fermi_fac * fermi_block

      
      !eas(k_n)  = kappa_n_analyt
      !eas(k_a)  = kappa_a_analyt
    end if 

    eas(k_a)  = kappa_a_analyt
    eas(k_n)  = kappa_n_analyt

  end subroutine m1_eas_analytic


end module mod_m1_eas
