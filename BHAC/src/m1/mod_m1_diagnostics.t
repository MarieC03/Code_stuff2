module mod_m1_diagnostics
    
   !************************

    contains
      
    {#IFNDEF UNIT_TESTS_EXPLICIT
    {#IFNDEF M1_EXPLICIT
    subroutine m1_compute_diagnostics(qdt,qtC,ixO^L,ixI^L,x,wprim,ynue,ynue_bar,mu_nue, mu_nue_bar, mu_delta_npe, mu_delta_npenu)
  
      use mod_m1_internal
      use mod_m1_constants
      use mod_eos_tabulated
      use mod_eos, only: small_rho, small_temp, big_ye , eos_yemin, eos_yemax
      use mod_m1_eas_param
      use mod_m1_closure
      
     {#IFNDEF UNIT_TESTS
     use mod_m1_metric_interface
     include "amrvacdef.f"
     }
     integer, intent(in) :: ixO^L, ixI^L
     double precision, intent(in) :: qdt,qtC
     double precision, intent(inout) :: wprim(ixI^S,1:nw)
     double precision, intent(in)  :: x(ixI^S, 1:ndim)
     double precision, intent(inout) :: mu_delta_npe(ixI^S), mu_delta_npenu(ixI^S)
     double precision, intent(inout)  :: mu_nue(ixI^S), mu_nue_bar(ixI^S)
     double precision, intent(inout)  :: ynue(ixI^S), ynue_bar(ixI^S)

     integer :: ix^D
     double precision, parameter :: Qnp1 = 1.2935 ! MeV !> rest-mass energy difference n&p (Ruffert 95)
    
     double precision :: T_fluid, rho, rho_cgs, y_e, y_mu, y_p
     double precision :: eps,prs,ent,cs2,dedt,dpderho,dpdrhoe
     double precision :: xa,xh,xn,xp,abar1,zbar1,mu_e,mu_n,mu_p,muhat,mu_nu, munu

     double precision, dimension(1:^NS) :: Gamma_closure,Wlor_closure,Jrad_closure
     double precision, dimension(1:^NS) :: N_KSP, nprimitive
     double precision, dimension(1:m1_numvars_internal, ^NS) :: wrad_KSP
     double precision, dimension(1:^NC) :: vel_
     double precision, dimension(1:^NC, 1:^NS) :: vel_closure

     type(m1_metric_helper) :: metricM1_help 
     type(m1_closure_helpers) :: stateM1_help 	

     {#IFNDEF UNIT_TESTS2
	   call metricM1_help%fill_metric(wprim,x,ixI^L,ixO^L)
     }
  
     
     {^KSP&
     {do ix^D=ixOmin^D,ixOmax^D \}

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

        nprimitive(^KSP) = N_KSP(^KSP)/Gamma_closure(^KSP)


        rho = wprim(ix^D,rho_) 
        T_fluid = wprim(ix^D,T_eps_) 
        y_e = wprim(ix^D,Ye_) 
        y_mu = 0.0
   
        {#IFDEF MUONS
          y_mu = wprim(ix^D,Ymu_) 
        \}
   
        y_p = y_e + y_mu
        rho_cgs = rho/RHOGF
   
        !! TODO if MUONS blabla
   
        {#IFNDEF UNIT_TESTS
        !> get the chemical potentials from EoS-table
        call tabulated_temp_get_all_one_grid(rho,T_fluid,y_e,eps,prs=prs,ent=ent,&
           cs2=cs2,dedt=dedt,dpderho=dpderho,dpdrhoe=dpdrhoe,xa=xa,xh=xh,xn=xn,xp=xp,&
           abar=abar1,zbar=zbar1,mu_e=mu_e,mu_n=mu_n,mu_p=mu_p,muhat=muhat,munu=munu)
        } 
     
       
        mu_delta_npe(ix^D) = mu_n - mu_p - mu_e + Qnp1
 
        if (rho_cgs > 1.0e10) then
           if(^KSP == m1_i_nue ) then
             ynue(ix^D)     = nprimitive(m1_i_nue)/ (rho)
             call m1_diagnostics_mu(rho, T_fluid, ynue(ix^D), mu_nue(ix^D))           
           else if(^KSP == m1_i_nuebar) then
             ynue_bar(ix^D) = nprimitive(m1_i_nuebar)/ (rho)
             call m1_diagnostics_mu(rho, T_fluid, ynue_bar(ix^D), mu_nue_bar(ix^D))
           endif 
        else
           ynue(ix^D) = 0.0d0
           ynue_bar(ix^D) = 0.0d0 
           mu_nue(ix^D) = 0.0d0
           mu_nue_bar(ix^D) = 0.0d0
        endif
        mu_delta_npenu(ix^D) = mu_n - mu_p - mu_e - mu_nue_bar(ix^D) + Qnp1

     {enddo ^D&\}
      \}
    end subroutine m1_compute_diagnostics

    subroutine m1_diagnostics_mu(rho,T_fluid,ynu,mu_out)
      use mod_m1_constants
      use mod_rootfinding, only: rootfinding_brent
      double precision :: rho,T_fluid,ynu
      double precision :: mu_out
       ! internal
      double precision :: etanu
      double precision :: F2_local
      integer :: flag

      F2_local = rho / RHOGF * (HC_MEVCM)**3 * ynu / (4.0d0 * PI * MNUC_CGS * T_fluid**3)

      if (ynu < 1e-20) then 
         mu_out = 0.0d0
      else
         call rootfinding_brent(etanu, -10.0d0, 10.0d0, 1.0d-10, 100, flag, Root_Fermi2)
         mu_out = etanu * T_fluid
      end if 

      contains

      function Root_Fermi2(etanu) result(res)
        use mod_m1_fermi
        implicit none
        double precision, intent(in) :: etanu
        double precision :: res
        double precision :: F2_etanu
      
        F2_etanu = fermi_dirac(etanu, 2)
        res = F2_local - F2_etanu
      end function Root_Fermi2

  end subroutine

    }
    }
 
end module mod_m1_diagnostics
!

