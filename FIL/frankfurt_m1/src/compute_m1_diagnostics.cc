//
//  Copyright (C) 2021, Harry Ho-Yin Ng
//  You reach npe beta eqm does not mean you reach npmu beta eqm at the same time => Separated into two cases.
         // mu_delta_npe              = \mu_n - \mu_p - \mu_e
         // mu_delta_npmu             = \mu_n - \mu_p - \mu_mu
         // mu_delta_npenue            = \mu_n - \mu_p - \mu_e - \mu_{\bar{\nu}_e} 
         // mu_delta_npmunumu           = \mu_n - \mu_p - \mu_mu - \mu_{\bar{\nu}_\mu}
         // tau_beta_nu[index]        = 1./sqrt(kappa_nu[index]*( kappa_nu[index] + kappa_nu[index] ) + 1.e-45);
         // tau_h = -\nabla_\mu u^\mu = \sqrtgamm D / \dot{\sqrtgamma D}  FIXME: how to get time derivative qunatites?
//Nnue = Nnue_star / sqrtgamma / Gamma =  neutrino number density = n_b * Y_nu

#include "m1_diagnostics.hh"
#include "cctk.h"
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <array>
#include <algorithm>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "frankfurt_m1.h" /* Generic #define's and function prototypes */
#include "Margherita_M1.h"
#include "utils/bb.hh"
#include "Margherita_EOS.h"

// From supplemental material of Pedro Espino+ 2024, most accurate in regions with Ynu > 1e-4 for rho > 1e12 gcm^-3
// Calculating \mu_\nu by a rootfinding with all cgs units:
double M1_Diagnostics::munu__rho_temp_ynu(const double &rho, const double &T, const double &ynu){
   using namespace M1_Constants;
   using namespace Margherita_constants;

   double F2, etanu;
   // Should I neutron mass, not amu? since I am using mass_n for rho --> here finding number density of baryon
   F2 = rho / Margherita_constants::RHOGF * ::int_pow<3>(hc_mevcm) * ynu / (4.0 * pi * massn_cgs * ::int_pow<3>(T));

   // Root finding interface closure
   auto const func = [&](double &eta_nu) {
     double F2_etanu = Fermi_Dirac<2>::get(eta_nu);
     return F2 - F2_etanu;
   };

   if (ynu < 1e-20) { 
      return 0.0;
   } else {
      etanu = zero_brent(-5.0, 5.0, 1.e-10, func);
      return etanu * T;
   }
};


extern "C" void compute_m1_diagnostics(CCTK_ARGUMENTS){
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace Margherita_M1_EAS ;
  using namespace Margherita_constants;
  if (use_m1_diagnostics) {

  #pragma omp parallel for collapse(3)
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
  
       int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
  
           // mu_delta_npe = \mu_n - \mu_p - \mu_e
           // mu_delta_npmu = \mu_n - \mu_p - \mu_mu
           double Xa, Xh, Xn, Xp, Abar, Zbar, T, rho, rho_cgs, y_e, y_mu;
           double prs, eps, dummy;
  
           rho = rho_b[index];
           rho_cgs = rho/RHOGF;
           T = temp[index];
           y_e = ye[index];
           y_mu = ymu[index];
           yp[index] = ye[index] + ymu[index];
  
           if (EOS_Leptonic::use_muonic_eos) {
              // bound yp
              yp[index] = std::min(EOS_Leptonic::eos_yemax, std::max(EOS_Leptonic::eos_yemin, yp[index]) );

              typename EOS_Leptonic::error_type error;
              mu_e[index] = EOS_Leptonic::mue_mumu_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_yle_ymu(
                                                      mu_mu[index], mu_p[index], mu_n[index], Xa, Xh, Xn, Xp,
                                                      Abar, Zbar, T, rho, y_e, y_mu, error);
  
              mu_delta_npe[index] = mu_n[index] - mu_p[index] - mu_e[index] + Qnp;
              mu_delta_npmu[index] = mu_n[index] - mu_p[index] - mu_mu[index] + Qnp;
  
              //if (rho/Margherita_constants::RHOGF > 1.0e11) {
              //  prs = EOS_Leptonic::press_eps_yle_ymu__nu_beta_eq__rho_temp(eps, y_e,
              //                                                              y_mu, rho, T, mu_nue_bar, mu_numu_bar,
              //                                                              error);
              //} else { 
              //  prs = 0.0;
              //  eps = 0.0;
              //}
  
              //P_nubetaeq[index] = prs;
              //eps_nubetaeq[index] = eps;
              //ye_nubetaeq[index] = y_e;
              //ymu_nubetaeq[index] = y_mu;


  	      double ynue_l, ynue_bar_l, ynumu_l, ynumu_bar_l;
 	      if (rho_cgs > 1.0e10) {
                ynue_l      = Nnue[index] / (rho);
                ynue_bar_l  = Nnue_bar[index] / (rho);
                ynumu_l     = Nnumu[index] / (rho);
                ynumu_bar_l = Nnumu_bar[index] / (rho);
	          
	        ynue[index]      = ynue_l;
	        ynue_bar[index]  = ynue_bar_l;
	        ynumu[index]     = ynumu_l;
	        ynumu_bar[index] = ynumu_bar_l;

              	mu_nue[index]      = M1_Diagnostics::munu__rho_temp_ynu(rho, T,
              	                           ynue_l);
              	mu_nue_bar[index]  = M1_Diagnostics::munu__rho_temp_ynu(rho, T,
              	                           ynue_bar_l);
              	mu_numu[index]     = M1_Diagnostics::munu__rho_temp_ynu(rho, T,
              	                           ynumu_l);
              	mu_numu_bar[index] = M1_Diagnostics::munu__rho_temp_ynu(rho, T,
              	                           ynumu_bar_l);
	      } else {
		ynue[index]        = 0.0;
		ynue_bar[index]    = 0.0;
		ynumu[index]       = 0.0;
		ynumu_bar[index]   = 0.0;
                mu_nue[index]      = 0.0;
                mu_nue_bar[index]  = 0.0;
                mu_numu[index]     = 0.0;
                mu_numu_bar[index] = 0.0;
              }

              mu_delta_npenu[index] = mu_n[index] - mu_p[index] - mu_e[index] - mu_nue_bar[index] + Qnp;
              mu_delta_npmunu[index] = mu_n[index] - mu_p[index] - mu_mu[index] - mu_numu_bar[index] + Qnp;

	     
  	      double yle_minus_l, yle_plus_l, ymu_minus_l, ymu_plus_l;
              double press_e_minus_l, press_e_plus_l, press_mu_minus_l, press_mu_plus_l;
              double press_b_l;
              double eps_e_minus_l, eps_e_plus_l, eps_mu_minus_l, eps_mu_plus_l, eps_b_l;

              yle_minus_l = EOS_Leptonic::yl_plus_minus_press_l_press_b__temp_rho_yle_ymu(
						yle_plus_l, ymu_minus_l, ymu_plus_l,
      						press_e_minus_l, press_e_plus_l, press_mu_minus_l, press_mu_plus_l, 
						press_b_l,
      						eps_e_minus_l, eps_e_plus_l, eps_mu_minus_l, eps_mu_plus_l, eps_b_l,
      						T, rho, y_e, y_mu,
      						error);
	      yle_minus[index] = yle_minus_l;
	      yle_plus[index] = yle_plus_l;
	      ymu_minus[index] = ymu_minus_l;
	      ymu_plus[index] = ymu_plus_l;
	      press_e_minus[index] = press_e_minus_l;
	      press_e_plus[index] = press_e_plus_l;
	      press_mu_minus[index] = press_mu_minus_l;
	      press_mu_plus[index] = press_mu_plus_l;
	      press_b[index] = press_b_l;
              eps_e_minus[index] = eps_e_minus_l;
              eps_e_plus[index] = eps_e_plus_l;
              eps_mu_minus[index] = eps_mu_minus_l;
              eps_mu_plus[index] = eps_mu_plus_l;
	      eps_b[index] = eps_b_l;


  
           } else {
              // bound yp
              yp[index] = std::min(EOS_Tabulated::eos_yemax, std::max(EOS_Tabulated::eos_yemin, yp[index]) );

              typename EOS_Tabulated::error_type error;
              mu_mu[index] = 0.0;
              mu_e[index] = EOS_Tabulated::mue_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_yle_ymu(
                                                  mu_p[index], mu_n[index], Xa, Xh, Xn, Xp,
                                                  Abar, Zbar, T, rho, y_e, y_mu, error);

              mu_delta_npe[index] = mu_n[index] - mu_p[index] - mu_e[index] + Qnp;

              double ynue_l, ynue_bar_l;
 	      if (rho_cgs > 1.0e10) {
                //ynue_l      = Nnue[index] * ::int_pow<3>(LENGTHGF)  / (rho_cgs);
                //ynue_bar_l  = Nnue_bar[index] * ::int_pow<3>(LENGTHGF)  / (rho_cgs);
                ynue_l      = Nnue[index]   / (rho);
                ynue_bar_l  = Nnue_bar[index]  / (rho);
                //ynue_l      = Nnue[index] * mnuc_Msun / (rho);
                //ynue_bar_l  = Nnue_bar[index] * mnuc_Msun / (rho);

                ynue[index]      = ynue_l;
                ynue_bar[index]  = ynue_bar_l;

              	mu_nue[index]      = M1_Diagnostics::munu__rho_temp_ynu(rho, T,
              	                           ynue_l);
              	mu_nue_bar[index]  = M1_Diagnostics::munu__rho_temp_ynu(rho, T,
              	                           ynue_bar_l);
              } else {
	        ynue[index]       = 0.0;
	        ynue_bar[index]   = 0.0;
                mu_nue[index]     = 0.0;
                mu_nue_bar[index] = 0.0;
              }
              mu_delta_npenu[index] = mu_n[index] - mu_p[index] - mu_e[index] - mu_nue_bar[index] + Qnp;


              //press_b[index] = EOS_Tabulated::press__temp_rho_yle_ymu(T, rho, y_e, y_mu,
              //                                               error);
              //eps_b[index] = EOS_Tabulated::press__temp_rho_yle_ymu(T, rho, y_e, y_mu,
              //                                               error);

              yle_minus[index] = EOS_Tabulated::ye_plus_minus_press_e_press_b__temp_rho_yle_ymu(
                				 		 yle_plus[index],
                				 		 press_e_minus[index], press_e_plus[index], press_b[index],
                				 		 eps_e_minus[index], eps_e_plus[index], eps_b[index],
                				 		 T, rho, y_e, y_mu,
                				 		 error);


           }
    } //for loop
    //std::cout <<"Testing for M1_diagnostics with muons is finished"<<"\n";
    //CCTK_ERROR("End of m1_diagnostics tests") ;
  } //endif of m1_diagnostic

}

