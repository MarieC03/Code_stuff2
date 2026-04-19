module mod_m1_fermi
  use mod_m1_constants
  
  implicit none
  public
  
  ! Minimum eta threshold for switching between approximations
  double precision, parameter :: eta_min = 1.0d-3

  !************************************************************************
  
contains
 !===================================================================
 !==================== Fermi-Dirac Integrals using Takahashi et al. 1978
 !===================================================================
 
 ! Horner's method for polynomial evaluation
 ! Evaluates: c0 + c1*x + c2*x^2 + ... + cn*x^n
 pure function horner_2(x, c0, c1, c2) result(res)
   double precision, intent(in) :: x, c0, c1, c2
   double precision :: res
   res = c0 + x * (c1 + x * c2)
 end function horner_2
 
 pure function horner_3(x, c0, c1, c2, c3, c4) result(res)
   double precision, intent(in) :: x, c0, c1, c2, c3, c4
   double precision :: res
   res = c0 + x * (c1 + x * (c2 + x * (c3 + x * c4)))
 end function horner_3
 
 pure function horner_4(x, c0, c1, c2, c3, c4, c5, c6) result(res)
   double precision, intent(in) :: x, c0, c1, c2, c3, c4, c5, c6
   double precision :: res
   res = c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * c6)))))
 end function horner_4

 ! Fermi-Dirac integral of order N
 ! Uses approximate formulae from Takahashi et al. 1978
 pure function fermi_dirac(eta, N) result(res)
   double precision, intent(in) :: eta
   integer, intent(in) :: N
   double precision :: res
   
   select case(N)
     case(0)
       res = fermi_dirac_0(eta)
     case(1)
       res = fermi_dirac_1(eta)
     case(2)
       res = fermi_dirac_2(eta)
     case(3)
       res = fermi_dirac_3(eta)
     case(4)
       res = fermi_dirac_4(eta)
     case(5)
       res = fermi_dirac_5(eta)
     case default
       ! Fallback for unsupported orders
       res = 0.0d0
   end select
 end function fermi_dirac

 ! Fermi-Dirac integral F_0(eta)
 pure function fermi_dirac_0(eta) result(res)
   double precision, intent(in) :: eta
   double precision :: res
   
   if (eta > eta_min) then
     res = eta + log(1.0d0 + exp(-eta))
   else
     res = log(1.0d0 + exp(eta))
   end if
 end function fermi_dirac_0

 ! Fermi-Dirac integral F_1(eta)
 pure function fermi_dirac_1(eta) result(res)
   double precision, intent(in) :: eta
   double precision :: res
   
   if (eta > eta_min) then
     res = horner_2(eta, 1.6449d0, 0.0d0, 0.5d0) / (1.0d0 + exp(-1.6855d0 * eta))
   else
     res = exp(eta) / (1.0d0 + 0.2159d0 * exp(0.8857d0 * eta))
   end if
 end function fermi_dirac_1

 ! Fermi-Dirac integral F_2(eta)
 pure function fermi_dirac_2(eta) result(res)
   double precision, intent(in) :: eta
   double precision :: res
   
   if (eta > eta_min) then
     res = horner_2(eta, 0.0d0, 3.2899d0, 0.0d0 + 1.0d0/3.0d0 * eta) / &
           (1.0d0 - exp(-1.8246d0 * eta))
     ! Simpler form:
     res = (3.2899d0 * eta + eta**3 / 3.0d0) / (1.0d0 - exp(-1.8246d0 * eta))
   else
     res = 2.0d0 * exp(eta) / (1.0d0 + 0.1092d0 * exp(0.8908d0 * eta))
   end if
 end function fermi_dirac_2

 ! Fermi-Dirac integral F_3(eta)
 pure function fermi_dirac_3(eta) result(res)
   double precision, intent(in) :: eta
   double precision :: res
   
   if (eta > eta_min) then
     res = horner_3(eta, 11.3644d0, 0.0d0, 4.9348d0, 0.0d0, 0.25d0) / &
           (1.0d0 + exp(-1.9039d0 * eta))
   else
     res = 6.0d0 * exp(eta) / (1.0d0 + 0.0559d0 * exp(0.9069d0 * eta))
   end if
 end function fermi_dirac_3

 ! Fermi-Dirac integral F_4(eta)
 pure function fermi_dirac_4(eta) result(res)
   double precision, intent(in) :: eta
   double precision :: res
   
   if (eta > eta_min) then
     res = horner_3(eta, 0.0d0, 45.4576d0, 0.0d0, 6.5797d0, 0.2d0 * eta) / &
           (1.0d0 - exp(-1.9484d0 * eta))
     ! Simpler form:
     res = (45.4576d0 * eta + 6.5797d0 * eta**3 + 0.2d0 * eta**5) / &
           (1.0d0 - exp(-1.9484d0 * eta))
   else
     res = 24.0d0 * exp(eta) / (1.0d0 + 0.0287d0 * exp(0.9257d0 * eta))
   end if
 end function fermi_dirac_4

 ! Fermi-Dirac integral F_5(eta)
 pure function fermi_dirac_5(eta) result(res)
   double precision, intent(in) :: eta
   double precision :: res
   
   if (eta > eta_min) then
     res = horner_4(eta, 236.5323d0, 0.0d0, 113.6439d0, 0.0d0, 8.2247d0, 0.0d0, 1.0d0/6.0d0) / &
           (1.0d0 + exp(-1.9727d0 * eta))
   else
     res = 120.0d0 * exp(eta) / (1.0d0 + 0.0147d0 * exp(0.9431d0 * eta))
   end if
 end function fermi_dirac_5

 !===============================================================
 ! Your existing subroutines with minimal changes
 !===============================================================
 
 subroutine get_neutrino_temp_ixD(wrad,N,ix^D,Gamma,velU,Wlor,Jrad, eta, speciesKSP, fluid_Prim, av_energy, T_nu, T_fluid, timeCurrent1,x)    
   use mod_m1_internal
   use mod_m1_eas_param 
   use mod_Weakhub_reader, only: logtemp_min_IV
   use ieee_arithmetic
   {#IFDEF UNIT_TESTS
   use mod_m1_tests
   }
   {#IFNDEF UNIT_TESTS
   include "amrvacdef.f"
   }
   integer, intent(in) :: ix^D
   integer, intent(in) :: speciesKSP
   double precision, intent(inout) :: wrad(1:m1_numvars_internal)
   double precision, intent(inout) :: fluid_Prim(1:3)
   double precision, intent(in)    :: Gamma, Wlor, Jrad
   double precision, intent(in)    :: velU(1:^NC) 
   double precision, intent(in)    :: N
   double precision, intent(in)    :: eta
   double precision, intent(out)   :: av_energy
   double precision, intent(out)   :: T_nu
   double precision, intent(in)    :: T_fluid
   double precision, intent(in)    :: timeCurrent1
   double precision, intent(in)    :: x(1:ndim)
   ! internal:
   double precision :: E
   double precision :: n_dens
   double precision :: Fv
   double precision :: av_energy_MEV
   double precision :: av_energy1
   double precision :: F3

     Fv  = {^C&wrad(m1_flux^C_)*velU(^C)+}

     E = wrad(m1_energy_) 
     n_dens = N/Gamma

      ! neutrinos atmosphere?
      if((n_dens .lt. m1_E_atmo*100) .or. (E .lt. m1_E_atmo*100) .or. &
         (T_fluid .le. exp(logtemp_min_IV)*1.05)) then
           av_energy = m1_E_atmo
           T_nu = 1.d-2
           return
      endif

        av_energy1 = Wlor*(E-Fv)/N
        av_energy_MEV = av_energy1
        av_energy = av_energy_MEV

        F3 = fermi_dirac(eta, 3)
        if (F3 < 1.0d-15) then
          T_nu = 1.0d0/3.0d0 * av_energy_MEV
          return
        endif

        !> neutrino temperature: now in MeV
        T_nu = fermi_dirac(eta, 2) / F3 * av_energy_MEV

        if (.not. ieee_is_finite(T_nu)) then
          write(*,*) "T_nu is NaN/Inf!"
          write(*,*) "eta=", eta
          write(*,*) "fermi_dirac(eta,2)=", fermi_dirac(eta, 2)
          write(*,*) "fermi_dirac(eta,3)=", fermi_dirac(eta, 3)
          write(*,*) "av_energy_MEV=", av_energy_MEV
          write(*,*) "E=", E
          write(*,*) "Fv=", Fv
          write(*,*) "N=", N
          write(*,*) "n_dens=", n_dens
        end if
 end subroutine get_neutrino_temp_ixD 

 subroutine blackbody_ixD(wrad,ix^D,speciesKSP,fluid_Prim,T_fluid,eta,Black_er,Black_nr)
      use mod_m1_internal 
      {#IFDEF UNIT_TESTS
      use mod_m1_tests
      }
      use mod_m1_eas_param
      {#IFNDEF UNIT_TESTS
      include "amrvacdef.f"
      }
      integer, intent(in) :: ix^D, speciesKSP
      double precision, intent(in) :: eta
      double precision, intent(in) :: T_fluid 
      double precision, intent(inout) :: wrad(1:m1_numvars_internal)
      double precision, intent(in) :: fluid_Prim(1:3)
      double precision, intent(out) :: Black_er,Black_nr
      ! internal: 
      double precision, dimension(1:6) :: g_fact
      double precision :: Black_er_MEV, Black_nr_MEV
      double precision :: Units, Units2, UnitsCGS_er, UnitsCGS_nr
   
      Units2 = 4*PI/C2_CGS/HC_MEVCM

      Units = 4*PI*CLITE_CM/(HC_MEVCM)**3
      
      if(m1_use_muons) then
        if(m1_use_neutrinos) then
            g_fact(m1_i_nue) = 1.0d0
            g_fact(m1_i_nuebar) = 1.0d0
            g_fact(m1_i_nux) = 2.0d0
            g_fact(m1_i_mu) = 1.0d0
            g_fact(m1_i_mubar) = 1.0d0
        else
            g_fact(m1_i_mu) = 1.0d0
            g_fact(m1_i_mubar) = 1.0d0
        end if 
      else
        if(m1_use_neutrinos) then
            g_fact(m1_i_nue) = 1.0d0
            g_fact(m1_i_nuebar) = 1.0d0
            g_fact(m1_i_nux) = 4.0d0
            g_fact(4) = 1.0d0
            g_fact(5) = 1.0d0
        end if 
      end if 

       !> blackbody radiation in MeV:
       Black_nr_MEV = g_fact(speciesKSP) * Units * T_fluid**3 * fermi_dirac(eta, 2)
       Black_er_MEV = g_fact(speciesKSP) * Units * T_fluid**4 * fermi_dirac(eta, 3)

       !> blackbody radiation in code units:
       Black_nr = Black_nr_MEV * MNUC_CGS * RHOGF/TIMEGF * LENGTHGF 
       Black_er = Black_er_MEV * MEV_TO_ERG * EPSGF * RHOGF/TIMEGF * LENGTHGF 

   end subroutine blackbody_ixD

end module mod_m1_fermi



