/*
 * =====================================================================================
 *
 *       Filename:  fugacities.h
 *
 *    Description:  Class to store neutrino fugacities
 *
 *        Version:  1.0
 *        Created:  13/05/2017 21:14:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#include <cmath>
#include <limits>
#include "M1_constants.hh"
#include "../margherita.hh"
#include "../Margherita_EOS.h"
#include "inv_fermi.hxx"
#include "helpers.hh"
#include "constexpr_utils.hh"
#ifndef __H_FUGACITY
#define __H_FUGACITY

template<typename T>
class Fugacities {
 public:
//Harry: if I change it back to nue, nua, nux --> only nux has problem
   //                                        proton = 4 , neutron = 5
  //enum { NUE = 0, NUE_BAR, NUX, ELECTRON, PROTON, NEUTRON, HAT, NUM_INDICES };
  //                                                      proton = 7, neutron = 8 
  enum { NUE = 0, NUE_BAR, NUMU, NUMU_BAR, NUX, ELECTRON, MUON, PROTON, NEUTRON, HAT, NUM_INDICES };
  enum {
    ALPHA = 0,
    HEAVY,
    ABAR,
    ZBAR,
    MNUM_INDICES = ZBAR + 1 + 5 // PROTON is 7, NEUTRON = 8
    //MNUM_INDICES = ZBAR + 1 + 2
  };  // ALSO account for PROTON and NEUTRON
  std::array<T, NUM_INDICES> eta;
  std::array<T, MNUM_INDICES> Ym;
//Harry: why here is MNUM_INDICES???? WTF
  //std::array<T, MNUM_INDICES> ftau;
  std::array<T, NUM_INDICES> ftau;

  T eta_np;
  T eta_pn;
  T nb;

  T rho;
  T temp;
  T ye;
  T ymu;

  std::array<T, 5> tau_n{{std::numeric_limits<T>::quiet_NaN(),
                          std::numeric_limits<T>::quiet_NaN(),
                          std::numeric_limits<T>::quiet_NaN(),
                          std::numeric_limits<T>::quiet_NaN(),
                          std::numeric_limits<T>::quiet_NaN()}};
  std::array<T, 5> tau_e;

  inline void compute_fugacities() {
    using namespace M1_Constants;
    using namespace Margherita_constants;
    using namespace Margherita_helpers;

    // Making this call will automatically enforce table bounds
    double dummy;


    if (EOS_Leptonic::use_muonic_eos) {
         typename EOS_Leptonic::error_type error;
         eta[ELECTRON] =
             EOS_Leptonic::mue_mumu_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_yle_ymu(eta[MUON], 
                 eta[PROTON], eta[NEUTRON], Ym[ALPHA], Ym[HEAVY], Ym[NEUTRON],
                 Ym[PROTON], Ym[ABAR], Ym[ZBAR], temp, rho, ye, ymu, error);

         // We need to convert the density to cgs now
         rho *= INVRHOGF * (mnuc_cgs / EOS_Leptonic::baryon_mass);
    } else {
   	 typename EOS_Tabulated::error_type error;
   	 eta[ELECTRON] =
   	     EOS_Tabulated::mue_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_yle_ymu(
   	         eta[PROTON], eta[NEUTRON], Ym[ALPHA], Ym[HEAVY], Ym[NEUTRON],
   	         Ym[PROTON], Ym[ABAR], Ym[ZBAR], temp, rho, ye, dummy, error);
   	 /*eta[ELECTRON] = 1.; eta[PROTON] = 1.; eta[NEUTRON] = 1.; 
   	   Ym[ABAR] = 1.; Ym[ZBAR] = 1.;*/

   	 // We need to convert the density to cgs now
   	 rho *= INVRHOGF * (mnuc_cgs / EOS_Tabulated::baryon_mass);
    } // use_muonic_eos
    /*
    auto const nN = rho / mnuc_cgs * Ym[NEUTRON] ;
    auto const nP = rho / mnuc_cgs * Ym[PROTON] ;
    auto const eta_const = ::int_pow<3>(hc_mevcm)/(4.*M_PI*std::pow(2.*mnuc_cgs*clite*clite*temp,1.5));
    // compute equilibrium chemical potentials for degenerate
    // non-relativistic Fermi gas -- to be used for the computation
    // of the nucleon phase space blocking factors
    auto const etaN = inv_Fermi_Dirac_1h(eta_const * nN );
    auto const etaP = inv_Fermi_Dirac_1h(eta_const*nP) ;
    */
    // Eta still stores chemical potentials, we now fill the remaining gaps
    // and then normalize by T

//Harry: using leptonic eos, will be very easily having unphysical mu_e > 1000 mev
// if eta < - 1e2,  any problem? let me fix here for eta_nue_bar and eta_numu_bar

    eta[HAT] = eta[NEUTRON] - eta[PROTON] - Qnp;
    eta[NUE] = eta[ELECTRON] + eta[PROTON] - eta[NEUTRON] - Qnp;
    eta[NUE_BAR] = -eta[NUE];

    if (EOS_Leptonic::use_muonic_eos) {
    	eta[NUMU] = eta[MUON] + eta[PROTON] - eta[NEUTRON] - Qnp;
    	eta[NUMU_BAR] = -eta[NUMU];

        // WARNING: EXTENSION OF Bollig's trick for grey scheme to suppress muonic processes
       // this still leads to very wrong epsnumu
        //if (rho < 1.e11) {
        //   eta[NUMU] = 0.0;
        //   eta[NUMU_BAR] = 0.0;
        //}

    } else {
    	eta[NUMU] = 0;
    	eta[NUMU_BAR] = -eta[NUMU];
    } // use_muonic_eos



    eta[NUX] = 0;

    eta[PROTON] += Qnp;

    // Interpolate between optically thin and thick limit

    // Foucart 2015 (B7)
    for(int i=NUE; i <= NUX; i++) {
      if(tau_n[i] > 2) {
	      ftau[i] = 1;
      } else if(tau_n[i] < 2./3.) {
	      ftau[i] = 0.;
      } else {
	      ftau[i] = 3. * (tau_n[i] - 2./3. ) / 4.;
      }
    }

    // Nan check for tau_nu, usually having NaN and then everywhere = 0...
    for (int i = NUE; i <= NUX; i++) {
        if (tau_n[i] != tau_n[i]) {tau_n[i] = 1.e-10;}
        if (tau_n[i] < 1e-10) {tau_n[i] = 1.e-10;}
    }

    // Now we obtain the fugacities
    for (auto &etaL : eta) etaL /= temp;


    // code test, without using optical depth:
    //const auto eta_e_fac = 1.0;    
    const auto eta_e_fac = (1. - exp(-tau_n[NUE]));
    if (std::isfinite(eta_e_fac)) {
      eta[NUE] *= eta_e_fac;
    }
    //const auto eta_a_fac = 1.0;    
    const auto eta_a_fac = (1. - exp(-tau_n[NUE_BAR]));
    if (std::isfinite(eta_a_fac)) {
      eta[NUE_BAR] *= eta_a_fac;
    }

    // code test
   const double rescale_factor_tau  = 1./1.;
 
    if (EOS_Leptonic::use_muonic_eos) {
    	//const auto eta_numu_fac = 1.0;
    	const auto eta_numu_fac = (1. - exp(-tau_n[NUMU] * rescale_factor_tau ));
    	if (std::isfinite(eta_numu_fac)) {
    	  eta[NUMU] *= eta_numu_fac;
    	}

    	//const auto eta_numu_bar_fac = 1.0;
    	const auto eta_numu_bar_fac = (1. - exp(-tau_n[NUMU_BAR] * rescale_factor_tau));
    	if (std::isfinite(eta_numu_bar_fac)) {
    	  eta[NUMU_BAR] *= eta_numu_bar_fac;
    	}
    } // use_muonic_table


    if (EOS_Leptonic::use_muonic_eos) {
    // code test, when using tau_nu to kill eta_numu, but still have a small layer of large values when ymu > 0.005, that region made 
    // a lot of things wrong, I need to use ymu to limit eta_numu instead of using rho < 1e11. 
    // it should fix: 1.  pair processes BB problems, 2. numu_fact in corrections of kappa problems, 3.: 
       //if (ymu < 0.005 || tau_n[NUMU] < 10.0  ) { 
       //if (ymu < 0.005 || tau_n[NUMU] < 20.0 || tau_n[NUMU_BAR] < 20.0 ) {   // advopac needs taunu 30
       //if (ymu < 0.005 || tau_n[NUMU] < 1. || tau_n[NUMU_BAR] < 1. ) { 
       //if (tau_n[NUMU] < 10. || tau_n[NUMU_BAR] < 10. ) { 
       ////if (ymu < 0.005 || tau_n[NUMU] < 100.0 || tau_n[NUMU_BAR] < 50.0)  {
       ////if (ymu < 0.005 && (tau_n[NUMU] < 10.0 || tau_n[NUMU_BAR] < 10.0) ) {
       ////if (ymu < 0.005 && (tau_n[NUMU] < 10.0 || tau_n[NUMU_BAR] < 10.0) ) {
       ////if (ymu < 0.005) {
       ////if (tau_n[NUMU] < 10.0 || tau_n[NUMU_BAR] < 10.0) {
       //   eta[NUMU] = 0.0;
       //   eta[NUMU_BAR] = 0.0;
       //}
       const double tau_nu_thr = 1.;
       if (tau_n[NUMU] < tau_nu_thr) {
          eta[NUMU] = 0.0;
       }
       if (tau_n[NUMU_BAR] < tau_nu_thr) {
          eta[NUMU_BAR] = 0.0;
       }
       // density thr to prevent NaN tau_n to keep high values of eta_numu numubar outside
       // Never use temp to be a threshold
       //if (rho < 1.e11 ) {
       if (rho < 1.e11 || ymu < 0.005) {
          eta[NUMU] = 0.0;
          eta[NUMU_BAR] = 0.0;
       }
    // code test bounding eta_nu, no need bound eta_nux
    // if limit -5,5 in BNS case, a lot of problems
       for (int i = NUMU; i < NUX; i++) {
            if (eta[i] > 0.) {
                eta[i] = std::min(eta[i], 5.e0);
            } else if (eta[i] < 0.) {
                eta[i] = std::max(eta[i], -5.e0);
            }
       }
    }

    // Enforce positivity of mass fractions
    for (auto &YL : Ym) YL = max(0, YL);

    nb = rho * avogadro;  // Neutrino number density

    // nucleon final state blocking, see Rosswog
    // Ruffert et al. (A13)
    eta_np = nb * (Ym[PROTON] - Ym[NEUTRON]) / (exp(-eta[HAT]) - 1.);
    eta_pn = nb * (Ym[NEUTRON] - Ym[PROTON]) / (exp(eta[HAT]) - 1.);
    // eta_np = exp(eta[HAT])*eta_pn;

    // Bruenn (ApJSS 58 1985) (3.1): non degenerate matter limit.
    if (rho < 2.e11) {
      eta_pn = nb * Ym[PROTON];
      eta_np = nb * Ym[NEUTRON];
    }

    if (!std::isnormal(eta_np)) eta_np = 0.;
    if (!std::isnormal(eta_pn)) eta_pn = 0.;

    // Rosswog (A9)
    eta_pn = std::max(eta_pn, 0.);
    eta_np = std::max(eta_np, 0.);

  }

  Fugacities(const T &rho, const T &temp, const T &ye, const T &ymu)
      : rho(rho), temp(temp), ye(ye), ymu(ymu) {
    compute_fugacities();
  }

  template <typename F>
  Fugacities(const double &rho, const double &temp, const double &ye, const double &ymu, 
             F &&__tau_n)
      : rho(rho), temp(temp), ye(ye), ymu(ymu), tau_n(std::forward<F>(__tau_n)) {
    compute_fugacities();
  }


};
#endif
