module mod_m1_backreaction
    double precision :: t_back
   !************************

    contains
      
    {#IFNDEF UNIT_TESTS_EXPLICIT
    {#IFNDEF M1_EXPLICIT
    subroutine m1_add_backreaction(dtfactor,qdt,qtC,t_backreact,species,x,wcons,wrad,ix^D,ixI^L,eas_ixD,sources,N,vel_,metricM1,energy_good,number_good)
  
      use mod_m1_internal
      use mod_m1_eas
      use mod_eos, only: small_rho, small_temp, big_ye, eos_yemin 
      use mod_Weakhub_reader, only: logtemp_min_IV, logtemp_max_IV,ye_min_IV, ye_max_IV,logrho_min_IV,logrho_max_IV
      use mod_m1_eas_param
      
     {#IFNDEF UNIT_TESTS
     use mod_m1_metric_interface
     include "amrvacdef.f"
     }
  
     integer, intent(in) :: ix^D, ixI^L
     integer, intent(in) :: species   
     double precision, intent(in) :: qdt,qtC
     double precision, intent(in) :: t_backreact
     double precision, intent(in)    :: x(ixI^S, 1:ndim)
     double precision, intent(in)    :: dtfactor   !< Timestep factor 
     double precision, intent(inout) :: wcons(ixI^S,1:nw)
     double precision, intent(inout) :: wrad(1:m1_numvars_internal)
     double precision, intent(in)    :: eas_ixD(1:m1_num_eas)
     double precision, intent(inout) :: sources(1:m1_numvars_internal+1)   
     double precision, intent(inout) :: N
     double precision, intent(in)    :: vel_(1:^NC)
     type(m1_metric_helper), intent(in)  :: metricM1 
     logical, intent(inout) :: energy_good
     logical, intent(inout) :: number_good
     ! internal
     !double precision :: couple_fluid = 0.02d0
     double precision :: couple_tau = 1.0d-3
     double precision :: couple_momenta = 1.0d-3
     double precision :: couple_ye = 1.0d-4
     double precision :: couple_ymu = 1.0d-4
     double precision :: dE_source, dS_source1, dS_source2, dS_source3, dS_norm, dYeC_source, dYmuC_source
     double precision:: tiny, Snorm, fracE, fracS, rel_dye_den, fracYe_rel, absYe_cap, fracYe_abs
     double precision :: rel_dymu_den,fracYmu_rel, absYmu_cap, fracYmu_abs
     double precision:: capE, capS, capYe, capYmu, factorA
     double precision :: deltaErad, deltaN
     double precision :: Er,Fmag,facFrad
     double precision :: fmaxfact = 0.999d0 !0.99999d0 !0.999d0
     double precision :: deltaFrad(1:^NC)
     double precision :: tau_tmp, ye_tmp, dye_tmp, error_Ye, Delta_Ye
     double precision :: rel_error_Ye = 1.0d-3
     double precision :: damping
     double precision :: stepsBetweenDamp = 10.d0 !5.0d0 
     logical :: m1_apply_limits = .false. !whether to limit source changes and limit wrad
     logical :: m1_Iterative_Damping = .false. !.true.
     integer, dimension(1:6) :: signum
     integer :: i
  
     if(m1_use_neutrinos) then
      signum(m1_i_nue) = 1
      signum(m1_i_nuebar) = -1
      signum(m1_i_nux) = 0
     end if 
     if(m1_use_muons) then
      signum(m1_i_mu) = 1 !M1_muons_TODO
      signum(m1_i_mubar) = -1  !M1_muons_TODO
     end if 
     if(m1_use_photons)then
      signum(m1_i_photon) = 0  !M1_photons_TODO
     end if     
  
  
       if(wcons(ix^D,T_eps_) .gt. 1000.d0 ) then
         write(88,*)"------------------------"
         write(88,*)"DYSP ixD",ix^D
         write(88,*)"x",x(ix^D,1),x(ix^D,2),x(ix^D,3)
         write(88,*)"in correct kappa: at t", qtC
         write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
         write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
         write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
         write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
         write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
         do i=1,m1_numvars_internal+1
           write(88,*) "i sources", i, sources(i)
         end do 
         if(wcons(ix^D,T_eps_) .gt. 1000.d0 ) then
           write(88,*)"DANGER: fluid T too large!",wcons(ix^D,T_eps_)
         end if 
         {^C& write(88,*)"vel_i",vel_(^C) \}
         write(88,*)"wrad: N",N/metricM1%sqrtg(ix^D)
         write(88,*)"wrad: E",wrad(m1_energy_)/metricM1%sqrtg(ix^D)
         {^C& write(88,*)"wrad: Fi",wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) \}
         write(88,*)"wprim: T",wcons(ix^D,T_eps_)
         write(88,*)"wprim: Ye",wcons(ix^D,ye_)
         write(88,*)"wprim: Rho",wcons(ix^D,rho_)
         write(88,*)"wcons: tau",wcons(ix^D,tau_),wcons(ix^D,tau_) - qdt*dtfactor*sources(1)
         {^C& write(88,*)"wcons: sC",wcons(ix^D,s^C_) ,wcons(ix^D,s^C_) - qdt*dtfactor*sources(1+^C) \}
         write(88,*)"wcons: dye",wcons(ix^D,tau_),wcons(ix^D,dye_) - qdt*dtfactor * signum(species) * sources(1+^NC+1)   
       end if 
  
  
       do i=1,m1_numvars_internal+1
       if(sources(i) .ne. sources(i)) then
         write(88,*) "SOURCES ARE NAN"
         write(88,*) "i sources", i, sources(i)
         sources(i) = 0.0d0
       end if 
       end do 
  
  
       ! propose raw updates (fluid gets -sources)
       dE_source   = - qdt*dtfactor * sources(1)
       dS_source1  = - qdt*dtfactor * sources(1+1)
       dS_source2  = - qdt*dtfactor * sources(1+2)
       dS_source3  = - qdt*dtfactor * sources(1+3)   ! if 3D
       if ((species .eq. m1_i_nue) .or. (species .eq. m1_i_nuebar)) then
         dYeC_source = - qdt*dtfactor * signum(species) * sources(1+^NC+1)
       else
         dYeC_source = 0.0d0
       endif
  
       !if(m1_use_muons) then
       !     if((species .eq. m1_i_mu) .or. (species .eq. m1_i_mubar)) then
       !      dYmuC_source = - qdt*dtfactor * signum(species) * sources(1+^NC+2)
       !     else
       !       dYmuC_source = 0.0d0
       !     end if 
       !end if 
  
       damping = 1.0d0
       if(M1_Iterative_Damping) then
           if((qtC > t_backreact) .and. (qtC < t_backreact + stepsBetweenDamp*qdt)) then
             damping = 0.005d0
           else if((qtC > t_backreact+ stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 2.0 * stepsBetweenDamp*qdt)) then
             damping = 0.01d0
           else if((qtC > t_backreact+ 2.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 3.0 * stepsBetweenDamp*qdt)) then
             damping = 0.02d0
           else if((qtC > t_backreact+ 3.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 4.0 * stepsBetweenDamp*qdt)) then
             damping = 0.05d0
           else if((qtC > t_backreact+ 4.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 5.0*stepsBetweenDamp*qdt)) then
             damping = 0.1d0
           else if((qtC > t_backreact+ 5.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 6.0*stepsBetweenDamp*qdt)) then
             damping = 0.2d0
           else if((qtC > t_backreact+ 6.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 7.0*stepsBetweenDamp*qdt)) then
             damping = 0.3d0
           else if((qtC > t_backreact+ 7.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 8.0*stepsBetweenDamp*qdt)) then
             damping = 0.4d0
           else if((qtC > t_backreact+ 8.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 9.0*stepsBetweenDamp*qdt)) then
             damping = 0.5d0
           else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
             damping = 0.6d0
           else if((qtC > t_backreact+ 10.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 11.0*stepsBetweenDamp*qdt)) then
             damping = 0.7d0
           else if((qtC > t_backreact+ 11.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 12.0*stepsBetweenDamp*qdt)) then
             damping = 0.8d0
           else if((qtC > t_backreact+ 12.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 13.0*stepsBetweenDamp*qdt)) then
             damping = 0.9d0
           else if((qtC > t_backreact+ 13.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 14.0*stepsBetweenDamp*qdt)) then
             damping = 1.0d0
           end if 
       end if 
  
       !-------------------------------------------
       if(m1_apply_limits) then
          ! fractional changes (use vector norm for momentum)
          tiny = 1.0d-60
          Snorm    = sqrt(wcons(ix^D,s1_)**2 + wcons(ix^D,s2_)**2 + wcons(ix^D,s3_)**2)
          dS_norm  = sqrt(dS_source1**2 + dS_source2**2 + dS_source3**2)
          
          fracE  = dabs(dE_source) / max(tiny, dabs(wcons(ix^D,tau_)))
          fracS  = dS_norm / max(tiny, Snorm)
          fracS = 0.d0
          fracS = max(fracS, dabs(dS_source1) / max(tiny, dabs(wcons(ix^D,s1_))))
          fracS = max(fracS, dabs(dS_source2) / max(tiny, dabs(wcons(ix^D,s2_))))
          fracS = max(fracS, dabs(dS_source3) / max(tiny, dabs(wcons(ix^D,s3_)))) 
  
          
          ! for lepton: combine a relative cap (if dye is a conserved density)
          ! and an absolute cap on deltaYe per step to protect EOS
          {#IFDEF TABEOS
          rel_dye_den = max(tiny, dabs(wcons(ix^D,dye_)))
          fracYe_rel  = dabs(dYeC_source) / rel_dye_den
          absYe_cap   = 2.0d-3     ! e.g. limit absolute deltaYe per substep
          fracYe_abs  = 0.0d0
          if ((species .eq. m1_i_nue) .or. (species .eq. m1_i_nuebar)) then
            fracYe_abs = dabs(dYeC_source) / max(tiny, absYe_cap)  ! treat like a fraction vs cap
          endif      
          !if(m1_use_muons) then
          !   ! for lepton: combine a relative cap (if dye is a conserved density)
          !   ! and an absolute cap on deltaYmu per step to protect EOS
          !   rel_dymu_den = max(tiny, dabs(wcons(ix^D,dymu_)))
          !   fracYmu_rel  = dabs(dYmuC_source) / rel_dymu_den
          !   absYmu_cap   = 2.0d-3     ! e.g. limit absolute deltaYmu per substep
          !   fracYmu_abs  = 0.0d0
          !   if ((species .eq. m1_i_mu) .or. (species .eq. m1_i_mu)) then
          !     fracYmu_abs = dabs(dYmuC_source) / max(tiny, absYu_cap)  ! treat like a fraction vs cap
          !   endif
          !end if 
          }
          
          ! user-set caps
          capE  = couple_tau      ! e.g. 0.02
          capS  = couple_momenta    ! e.g. 0.02
          capYe = couple_ye         ! e.g. 0.02 (relative), plus abs cap above
          !capYmu = couple_ymu      ! e.g. 0.02 (relative), plus abs cap above
          
          factorA = 1.0d0
          factorA = min(factorA, capE  / max(tiny, fracE))
          factorA = min(factorA, capS  / max(tiny, fracS))
          if ((species .eq. m1_i_nue) .or. (species .eq. m1_i_nuebar)) then
            factorA = min(factorA, capYe / max(tiny, fracYe_rel))
            factorA = min(factorA, 1.0d0 / max(tiny, fracYe_abs))
          endif
          !if(m1_use_muons) then
          !  if ((species .eq. m1_i_mu) .or. (species .eq. m1_i_mu)) then
          !    factorA = min(factorA, capYmu / max(tiny, fracYmu_rel))
          !    factorA = min(factorA, 1.0d0 / max(tiny, fracYmu_abs))
          !  endif
          !end if 
  
          factorA = factorA * damping
          ! now apply the SAME factorA to all fluid channels
          wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * dE_source
          {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * dS_source^C \}
          {#IFDEF TABEOS
          if ((species .eq. m1_i_nue) .or. (species .eq. m1_i_nuebar)) then
            dye_tmp = wcons(ix^D,dye_) + factorA * dYeC_source
            ye_tmp = dye_tmp / wcons(ix^D,d_)
            Delta_Ye = dabs( ye_tmp - wcons(ix^D,dye_)/ wcons(ix^D,d_))
            error_Ye = Delta_Ye / rel_error_Ye
            if(ye_tmp < eos_yemin .or. ye_tmp > big_ye) then
               write(88,*) "m1-backreaction: ye is out of eos-bounds" 
               write(88,*) "ye yemin yemax",ye_tmp,eos_yemin,big_ye
            else if(error_Ye > 10.0d0) then
               write(88,*)"-----------------------------"
               write(88,*) "m1-backreaction: the relative ye-error is > 10"               
               write(88,*)"DYSP ixD",ix^D
               write(88,*)"x",x(ix^D,1),x(ix^D,2),x(ix^D,3)
               write(88,*)"in correct kappa: at t", qtC
               write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
               write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
               write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
               write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
               write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
               write(88,*)"wrad: N",N/metricM1%sqrtg(ix^D)
               write(88,*)"wrad: E",wrad(m1_energy_)/metricM1%sqrtg(ix^D)
               write(88,*) "dYe-source", dYeC_source
               {^C& write(88,*)"wrad: Fi",wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) \}
            else 
               wcons(ix^D,dye_) = wcons(ix^D,dye_) + factorA * dYeC_source
            end if 
            ! TEST !!!
            !wcons(ix^D,dye_) = wcons(ix^D,dye_) - factorA * dYeC_source
          endif
          !if(m1_use_muons) then
          !  if ((species .eq. m1_i_mu) .or. (species .eq. m1_i_mu)) then
          !    wcons(ix^D,dymu_) = wcons(ix^D,dymu_) + factorA * dYmuC_source
          !  endif
          !end if 
          }
          
       !-------------------------------------------
       else !not apply limits
          factorA = 1.0d0
          !factorA = 0.95d0
          !damping = 1.0d0
          factorA = factorA * damping
          ! now apply the SAME factorA to all fluid channels
          tau_tmp =  wcons(ix^D,tau_) + factorA * dE_source
          if(tau_tmp<0.0d0) then
             energy_good = .false.
             write(88,*) "m1-backreaction: energy tau is < 0"
          else
             energy_good = .true.
             wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * dE_source
          endif 
          ! TEST !!!
          !wcons(ix^D,tau_) = wcons(ix^D,tau_) - factorA * dE_source
          {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * dS_source^C \}
          ! TEST: !!
          !{C wcons(ix^D,s^C_) = wcons(ix^D,s^C_) - factorA * dS_source^C}
          {#IFDEF TABEOS
          if ((species .eq. m1_i_nue) .or. (species .eq. m1_i_nuebar)) then
            dye_tmp = wcons(ix^D,dye_) + factorA * dYeC_source
            ye_tmp = dye_tmp / wcons(ix^D,d_)
            Delta_Ye = dabs( ye_tmp - wcons(ix^D,dye_)/ wcons(ix^D,d_))
            error_Ye = Delta_Ye / rel_error_Ye
            if(ye_tmp < eos_yemin .or. ye_tmp > big_ye) then
              number_good = .false.
               write(88,*) "m1-backreaction: ye is out of eos-bounds" 
               write(88,*) "ye yemin yemax",ye_tmp,eos_yemin,big_ye
            else if(error_Ye > 10.0d0) then
              number_good = .false.
               write(88,*)"-----------------------------"
               write(88,*) "m1-backreaction: the relative ye-error is > 10"               
               write(88,*)"DYSP ixD",ix^D
               write(88,*)"x",x(ix^D,1),x(ix^D,2),x(ix^D,3)
               write(88,*)"in correct kappa: at t", qtC
               write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
               write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
               write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
               write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
               write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
               write(88,*)"wrad: N",N/metricM1%sqrtg(ix^D)
               write(88,*)"wrad: E",wrad(m1_energy_)/metricM1%sqrtg(ix^D)
               write(88,*) "dYe-source", dYeC_source
               {^C& write(88,*)"wrad: Fi",wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) \}
            else 
               number_good = .True.
               wcons(ix^D,dye_) = wcons(ix^D,dye_) + factorA * dYeC_source
            end if 
            ! TEST !!!
            !wcons(ix^D,dye_) = wcons(ix^D,dye_) - factorA * dYeC_source
          endif
          !if(m1_use_muons) then
          !  if ((species .eq. m1_i_mu) .or. (species .eq. m1_i_mu)) then
          !    wcons(ix^D,dymu_) = wcons(ix^D,dymu_) + factorA * dYmuC_source
          !  endif
          !end if 
          }
       end if 
  
  
    end subroutine m1_add_backreaction
    }
    }
 
end module mod_m1_backreaction
!

