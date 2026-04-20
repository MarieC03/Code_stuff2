module mod_m1_backreaction
    double precision :: t_back
   !************************

    contains
      
    {#IFNDEF UNIT_TESTS_EXPLICIT
    {#IFNDEF M1_EXPLICIT
    subroutine m1_add_backreaction(dtfactor,qdt,qtC,t_backreact,x,wcons,ix^D,ixI^L,sources,energy_good,number_good, proton_number_good)
  
      use mod_m1_internal
      use mod_m1_eas
      use mod_eos, only: small_rho, small_temp, big_ye, eos_yemin, eos_ymumin, eos_ymumax, &
           eos_uses_ye, eos_has_ymu
      use mod_Weakhub_reader, only: logtemp_min_IV, logtemp_max_IV,ye_min_IV, ye_max_IV,logrho_min_IV,logrho_max_IV
      use mod_m1_eas_param
      
     {#IFNDEF UNIT_TESTS
     use mod_m1_metric_interface
     include "amrvacdef.f"
     }
  
     integer, intent(in) :: ix^D, ixI^L
     double precision, intent(in) :: qdt,qtC
     double precision, intent(in) :: t_backreact
     double precision, intent(in)    :: x(ixI^S, 1:ndim)
     double precision, intent(in)    :: dtfactor   !< Timestep factor 
     double precision, intent(inout) :: wcons(ixI^S,1:nw)
     double precision, intent(inout) :: sources(1:m1_numvars_internal+1,^NS)   
     logical, intent(inout) :: energy_good
     logical, intent(inout) :: number_good
     logical, intent(inout) :: proton_number_good
     ! internal
     double precision :: sources_tmp, sources_dye_tmp, sources_dymu_tmp
     double precision :: tau_tmp, ye_tmp, ymu_tmp, dye_tmp, dymu_tmp, error_Ye, Delta_Ye
     double precision :: rel_error_Ye = 1.0d-3
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
  
     {^KSP&
     do i=1,m1_numvars_internal+1
     if(sources(i,^KSP) .ne. sources(i,^KSP)) then
       write(88,*) "SOURCES ARE NAN"
       write(88,*) "i sources", i, sources(i,^KSP)
       sources(i,^KSP) = 0.0d0
     end if 
     end do 
     \}
  
     !-------------------------------------------
     ! check if tau with sum of energy-sources is positive:
     sources_tmp = 0.0d0 
     {^KSP& sources_tmp  =  sources_tmp + sources(1,^KSP) \}
     tau_tmp =  wcons(ix^D,tau_) - qdt*dtfactor * sources_tmp 
     if(tau_tmp<0.0d0) then
        energy_good = .false.
     else
        energy_good = .true.
     endif 
     number_good = .true.
     proton_number_good = .true.
     {#IFDEF TABEOS
       if (eos_uses_ye()) then
         ! check if Ye with the sum of electron-flavor number sources is within range
         sources_dye_tmp = 0.0d0
         sources_dymu_tmp = 0.0d0
         {^KSP&
         if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
           sources_dye_tmp = sources_dye_tmp + signum(^KSP) * sources(1+^NC+1,^KSP)
         else if (m1_use_muons .and. ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mubar))) then
           sources_dymu_tmp = sources_dymu_tmp + signum(^KSP) * sources(1+^NC+1,^KSP)
         end if
         \}

         dye_tmp = wcons(ix^D,dye_) - qdt*dtfactor * sources_dye_tmp
         ye_tmp = dye_tmp / wcons(ix^D,d_)
         number_good = .true.
         if (ye_tmp < eos_yemin .or. ye_tmp > big_ye) then
           number_good = .false.
         end if

         proton_number_good = .true.
         if (eos_has_ymu()) then
           dymu_tmp = wcons(ix^D,dymu_) - qdt*dtfactor * sources_dymu_tmp
           ymu_tmp = dymu_tmp / wcons(ix^D,d_)
           if (ymu_tmp < eos_ymumin .or. ymu_tmp > eos_ymumax) then
             number_good = .false.
             proton_number_good = .false.
           end if
         end if
       else
         number_good = .true.
         proton_number_good = .true.
       end if
     }
     {^KSP&
       if(energy_good) then 
            wcons(ix^D,tau_) = wcons(ix^D,tau_) - qdt*dtfactor * sources(1,^KSP)
       end if 
       if(number_good) then
         if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
              wcons(ix^D,dye_) = wcons(ix^D,dye_) - qdt*dtfactor * signum(^KSP) * sources(1+^NC+1,^KSP)
         else if (eos_has_ymu() .and. ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mubar))) then
              wcons(ix^D,dymu_) = wcons(ix^D,dymu_) - qdt*dtfactor * signum(^KSP) * sources(1+^NC+1,^KSP)
         end if
       end if 
       
       {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) - qdt*dtfactor * sources(1+^C,^KSP) \}
     
     \}
  
  
    end subroutine m1_add_backreaction
    }
    }
 
end module mod_m1_backreaction
!

