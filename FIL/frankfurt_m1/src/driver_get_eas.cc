//Harry: all eas stuff should be <5>
/**
 * @file driver_get_eas.cc
 *
 * This file is part of frankfurt_m1
 *
 * The routines contained in this file
 * fill the GFs containing all the 
 * required microphysics for M1
 * these include emissivities (Q)
 * absorption and scattering opacities 
 * (k_a, k_s) number emission rates (R)
 * and transport opacities (k_n).
 *
 * The computations are done by Margherita,
 * here we just store the result.
 * If the beta-equilibrium fix a la Radice 
 * is enabled we additionally find (Y_e, T)_{\beta-eq}
 * and use these for the transport rates
 *
 * This file is heavily based on the older
 * version by Lukas Weih and Elias Most
 *
 * @author Carlo Musolino Elias R. Most Lukas Weih 
 */ 
#include "cctk.h"
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <array>
#include <algorithm>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "frankfurt_m1.h" /* Generic #define's and function prototypes */
#include "m1_closure.hh"
#include "Margherita_M1.h"
#include "Margherita_Weakhub.h"
#include "m1_closure.hh"
#include "utils/bb.hh"
//#define DEBUG_FILM1

void frankfurt_m1_compute_eas_microphysical(CCTK_ARGUMENTS);
void frankfurt_m1_compute_eas_corrected(CCTK_ARGUMENTS);
void frankfurt_m1_compute_eas_corrected_Weakhub(CCTK_ARGUMENTS);

extern "C" void frankfurt_m1_compute_eas(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if ( verbosity > 1 )
    CCTK_VInfo(CCTK_THORNSTRING,"***** Computing microphysics *****");

  if( CCTK_Equals( m1_which_eas, "microphysical" ) ) {
    frankfurt_m1_compute_eas_microphysical(CCTK_PASS_CTOC);
  } else if ( CCTK_Equals( m1_which_eas, "microphysical-corrected" ) ) {
    frankfurt_m1_compute_eas_corrected(CCTK_PASS_CTOC) ;
  } else if ( CCTK_Equals( m1_which_eas, "Weakhub" ) ) {
    frankfurt_m1_compute_eas_corrected_Weakhub(CCTK_PASS_CTOC) ;
  } else if (  CCTK_Equals( m1_which_eas, "test" ) ) {
    // this means we're running some kind of
    // test and the eas were set with the initial
    // data
    return ;
  }
  
}

void frankfurt_m1_compute_eas_microphysical(CCTK_ARGUMENTS) {
  	DECLARE_CCTK_ARGUMENTS;
  	DECLARE_CCTK_PARAMETERS;

	using namespace Margherita_M1_EAS ;

#pragma omp parallel for collapse(3)
	for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {

		int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

           if (use_5_spec_m1) {

		if( rho_b[index] < M1_rho_floor * ( 1.0 + M1_atmo_tol ) ) {


			Qnue[index] = 0. ;
			Qnue_bar[index] = 0. ;
			Qnumu[index] = 0. ;
			Qnumu_bar[index] = 0. ;
			Qnux[index] = 0. ;

			kappa_nue_a[index] = 0.;
			kappa_nue_bar_a[index] = 0.;
			kappa_numu_a[index] = 0.;
			kappa_numu_bar_a[index] = 0.;
			kappa_nux_a[index] = 0.;

			kappa_nue_s[index] = 0.;
			kappa_nue_bar_s[index] = 0.;
			kappa_numu_s[index] = 0.;
			kappa_numu_bar_s[index] = 0.;
			kappa_nux_s[index] = 0.;

			kappa_nue_n[index] = 0.;
			kappa_nue_bar_n[index] = 0.;
			kappa_numu_n[index] = 0.;
			kappa_numu_bar_n[index] = 0.;
			kappa_nux_n[index] = 0.;

			Rnue[index] = 0. ;
			Rnue_bar[index] = 0. ;
			Rnumu[index] = 0. ;
			Rnumu_bar[index] = 0. ;
			Rnux[index] = 0. ;


			continue ;
			
		}
		
		std::array<CCTK_REAL,REQ_FLUID_VARS> U;
		U[RHO]  = rho_b[index];
		U[YE]   = ye[index];
		U[YMU]  = ymu[index];
		U[TEMP] = temp[index];

		//The EAS class is provided by margherita with all the
		//functions and members called below. The EAS class is
		//basically the same as the Leakage class, but with some
		//modifications so that we get the correct opacities/emissivities
		//for M1.
		const auto tau_init = EAS<CCTK_REAL>::compute_analytic_opacity(U[RHO]);
		std::array<double,NUMSPECIES> tau_n {{tau_init, tau_init, tau_init, tau_init, tau_init}};
		if ( ! CCTK_Equals(m1_which_tau, "simple") ) {
		tau_n[NUE]      = tau_nue[index];
		tau_n[NUE_BAR]  = tau_nue_bar[index];
		tau_n[NUMU]     = tau_numu[index];
		tau_n[NUMU_BAR] = tau_numu_bar[index];
		tau_n[NUX]      = tau_nux[index];
		}
  //Harry: change it to read YMU and tau_n size 3 --> 5
		auto F= Fugacities<CCTK_REAL>(U[RHO], U[TEMP], U[YE], U[YMU], std::move(tau_n) );

		std::array<bool,M1_interactions::NINTERACTIONS> which_interactions { { beta_decay,
											plasmon_decay,
											bremstrahlung,
											pair_annihil }} ;
		
		EAS<CCTK_REAL> eas{F, which_interactions};
		eas.calc_eas(true);

		temp_nue[index]      = EAS<CCTK_REAL>::neutrino_temperature(eps_nue[index], F.eta[NUE]);
		temp_nue_bar[index]  = EAS<CCTK_REAL>::neutrino_temperature(eps_nue_bar[index], F.eta[NUE_BAR]);
		temp_numu[index]     = EAS<CCTK_REAL>::neutrino_temperature(eps_numu[index], F.eta[NUMU]);
		temp_numu_bar[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_numu_bar[index], F.eta[NUMU_BAR]);
		temp_nux[index]      = EAS<CCTK_REAL>::neutrino_temperature(eps_nux[index], F.eta[NUX]);

		CCTK_REAL const nue_fact = ::max( 1., ::int_pow<2>(temp_nue[index]/temp[index]));
		CCTK_REAL const nue_bar_fact = ::max( 1., ::int_pow<2>(temp_nue_bar[index]/temp[index]));
		CCTK_REAL const numu_fact = ::max( 1., ::int_pow<2>(temp_numu[index]/temp[index]));
		CCTK_REAL const numu_bar_fact = ::max( 1., ::int_pow<2>(temp_numu_bar[index]/temp[index]));
		CCTK_REAL const nux_fact = ::max( 1., ::int_pow<2>(temp_nux[index]/temp[index]));
		
		Qnue[index] = eas.Q_cactus<NUE>(); 
		Qnue_bar[index] = eas.Q_cactus<NUE_BAR>(); 
		Qnumu[index] = eas.Q_cactus<NUMU>(); 
		Qnumu_bar[index] = eas.Q_cactus<NUMU_BAR>(); 
		Qnux[index] = eas.Q_cactus<NUX>();
		
		kappa_nue_a[index] = eas.kappa_a_cactus<NUE>() * nue_fact ; 
		kappa_nue_bar_a[index] = eas.kappa_a_cactus<NUE_BAR>() * nue_bar_fact ; 
		kappa_numu_a[index] = eas.kappa_a_cactus<NUMU>() * numu_fact ; 
		kappa_numu_bar_a[index] = eas.kappa_a_cactus<NUMU_BAR>() * numu_bar_fact ; 
		kappa_nux_a[index] = eas.kappa_a_cactus<NUX>() ;
		
		kappa_nue_s[index] = eas.kappa_s_cactus<NUE>() * nue_fact ; 
		kappa_nue_bar_s[index] = eas.kappa_s_cactus<NUE_BAR>() * nue_bar_fact ; 
		kappa_numu_s[index] = eas.kappa_s_cactus<NUMU>() * numu_fact ; 
		kappa_numu_bar_s[index] = eas.kappa_s_cactus<NUMU_BAR>() * numu_bar_fact ; 
		kappa_nux_s[index] = eas.kappa_s_cactus<NUX>() * nux_fact ;

		Rnue[index] = eas.R_cactus<NUE>();
		Rnue_bar[index] = eas.R_cactus<NUE_BAR>();
		Rnumu[index] = eas.R_cactus<NUMU>();
		Rnumu_bar[index] = eas.R_cactus<NUMU_BAR>();
		Rnux[index] = eas.R_cactus<NUX>();
		
		kappa_nue_n[index] = eas.kappa_n_cactus<NUE>() * nue_fact ; 
		kappa_nue_bar_n[index] = eas.kappa_n_cactus<NUE_BAR>() * nue_bar_fact ; 
		kappa_numu_n[index] = eas.kappa_n_cactus<NUMU>() * numu_fact ; 
		kappa_numu_bar_n[index] = eas.kappa_n_cactus<NUMU_BAR>() * numu_bar_fact ; 
		kappa_nux_n[index] = eas.kappa_n_cactus<NUX>() ;

           } else {

		if( rho_b[index] < M1_rho_floor * ( 1.0 + M1_atmo_tol ) ) {


			Qnue[index] = 0. ;
			Qnue_bar[index] = 0. ;
			Qnux[index] = 0. ;

			kappa_nue_a[index] = 0.;
			kappa_nue_bar_a[index] = 0.;
			kappa_nux_a[index] = 0.;

			kappa_nue_s[index] = 0.;
			kappa_nue_bar_s[index] = 0.;
			kappa_nux_s[index] = 0.;

			kappa_nue_n[index] = 0.;
			kappa_nue_bar_n[index] = 0.;
			kappa_nux_n[index] = 0.;

			Rnue[index] = 0. ;
			Rnue_bar[index] = 0. ;
			Rnux[index] = 0. ;


			continue ;
			
		}
		
		std::array<CCTK_REAL,REQ_FLUID_VARS> U;
		U[RHO]  = rho_b[index];
		U[YE]   = ye[index];
		U[YMU]  = ymu[index];
		U[TEMP] = temp[index];

		//The EAS class is provided by margherita with all the
		//functions and members called below. The EAS class is
		//basically the same as the Leakage class, but with some
		//modifications so that we get the correct opacities/emissivities
		//for M1.
		const auto tau_init = EAS<CCTK_REAL>::compute_analytic_opacity(U[RHO]);
		std::array<double,NUMSPECIES> tau_n {{tau_init, tau_init, tau_init, tau_init, tau_init}};
		//std::array<double,3> tau_n {{tau_init, tau_init, tau_init}};
		if ( ! CCTK_Equals(m1_which_tau, "simple") ) {
		tau_n[NUE] = tau_nue[index];
		tau_n[NUE_BAR] = tau_nue_bar[index];
		tau_n[NUMU] = 0.0;
		tau_n[NUMU_BAR] = 0.0;
		tau_n[NUX] = tau_nux[index];
		}
		auto F= Fugacities<CCTK_REAL>(U[RHO], U[TEMP], U[YE], U[YMU], std::move(tau_n) );

		std::array<bool,M1_interactions::NINTERACTIONS> which_interactions { { beta_decay,
											plasmon_decay,
											bremstrahlung,
											pair_annihil }} ;
		
		EAS<CCTK_REAL> eas{F, which_interactions};
		eas.calc_eas(true);

		temp_nue[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue[index], F.eta[NUE]);
		temp_nue_bar[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue_bar[index], F.eta[NUE_BAR]);
		temp_numu[index] = 0.0 ;
		temp_numu_bar[index] = 0.0 ;
		temp_nux[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nux[index], F.eta[NUX]);

		CCTK_REAL const nue_fact = ::max( 1., ::int_pow<2>(temp_nue[index]/temp[index]));
		CCTK_REAL const nue_bar_fact = ::max( 1., ::int_pow<2>(temp_nue_bar[index]/temp[index]));
		CCTK_REAL const nux_fact = ::max( 1., ::int_pow<2>(temp_nux[index]/temp[index]));
		
		Qnue[index] = eas.Q_cactus<NUE>(); 
		Qnue_bar[index] = eas.Q_cactus<NUE_BAR>(); 
		Qnumu[index] = 0.0 ;
		Qnumu_bar[index] = 0.0 ;
		Qnux[index] = eas.Q_cactus<NUX>();
		
		kappa_nue_a[index] = eas.kappa_a_cactus<NUE>() * nue_fact ; 
		kappa_nue_bar_a[index] = eas.kappa_a_cactus<NUE_BAR>() * nue_bar_fact ; 
		kappa_numu_a[index] = 0.0 ;
		kappa_numu_bar_a[index] = 0.0 ;
		kappa_nux_a[index] = eas.kappa_a_cactus<NUX>() ;
		
		kappa_nue_s[index] = eas.kappa_s_cactus<NUE>() * nue_fact ; 
		kappa_nue_bar_s[index] = eas.kappa_s_cactus<NUE_BAR>() * nue_bar_fact ; 
		kappa_numu_s[index] = 0.0 ;
		kappa_numu_bar_s[index] = 0.0 ; // Harry kappa_s_numu numu bar should be non-zero, but store inside nux for 3 spec
		kappa_nux_s[index] = eas.kappa_s_cactus<NUX>() * nux_fact ;

		Rnue[index] = eas.R_cactus<NUE>();
		Rnue_bar[index] = eas.R_cactus<NUE_BAR>();
		Rnumu[index] = 0.0 ;
		Rnumu_bar[index] = 0.0 ;
		Rnux[index] = eas.R_cactus<NUX>();
		
		kappa_nue_n[index] = eas.kappa_n_cactus<NUE>() * nue_fact ; 
		kappa_nue_bar_n[index] = eas.kappa_n_cactus<NUE_BAR>() * nue_bar_fact ; 
		kappa_numu_n[index] = 0.0 ;
		kappa_numu_bar_n[index] = 0.0 ;
		kappa_nux_n[index] = eas.kappa_n_cactus<NUX>() ;

                // using 3 spec, numu numubar zeros:

                Qnumu[index] = 0. ;
                Qnumu_bar[index] = 0. ;

                kappa_numu_a[index] = 0.;
                kappa_numu_bar_a[index] = 0.;

                kappa_numu_s[index] = 0.;
                kappa_numu_bar_s[index] = 0.;

                kappa_numu_n[index] = 0.;
                kappa_numu_bar_n[index] = 0.;

                Rnumu[index] = 0. ;
                Rnumu_bar[index] = 0. ;
	   } // use_5_spec_m1
	} //for loop

} //compute_eas_microphysical

	// Compute transport rate but use T,Y_e at beta equilibrium
	// if the equilibration timescale is shorter than 0.5 * \Delta t
	// This requires a 2D rootfinding that is handled within Margherita
	// (see M1_find_betaeq.hh)
void frankfurt_m1_compute_eas_corrected(CCTK_ARGUMENTS) {
	DECLARE_CCTK_ARGUMENTS;
	DECLARE_CCTK_PARAMETERS;

	using namespace M1_Constants;
	using namespace Margherita_constants ;
	using namespace Margherita_M1_EAS ;
	using namespace frankfurt_m1 ; 
	//Global quantities for computing avergae over entire grid
	
	CCTK_INT nfix=0, nsucc=0;
	CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2] ;
	CCTK_REAL avgiter=0., error_int=0., error_int_denom=0.  ;

	
	#pragma omp parallel for collapse(3) reduction(+:nfix,nsucc,avgiter,error_int,error_int_denom) schedule(static)
	for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {

		int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

           if (use_5_spec_m1) {

		if( rho_b[index] < M1_rho_floor * ( 1.0 + M1_atmo_tol ) ) {

			Qnue[index] = 0. ;
			Qnue_bar[index] = 0. ;
			Qnumu[index] = 0. ;
			Qnumu_bar[index] = 0. ;
			Qnux[index] = 0. ;

			kappa_nue_a[index] = 0.;
			kappa_nue_bar_a[index] = 0.;
			kappa_numu_a[index] = 0.;
			kappa_numu_bar_a[index] = 0.;
			kappa_nux_a[index] = 0.;

			kappa_nue_s[index] = 0.;
			kappa_nue_bar_s[index] = 0.;
			kappa_numu_s[index] = 0.;
			kappa_numu_bar_s[index] = 0.;
			kappa_nux_s[index] = 0.;

			kappa_nue_n[index] = 0.;
			kappa_nue_bar_n[index] = 0.;
			kappa_numu_n[index] = 0.;
			kappa_numu_bar_n[index] = 0.;
			kappa_nux_n[index] = 0.;

			Rnue[index] = 0. ;
			Rnue_bar[index] = 0. ;
			Rnumu[index] = 0. ;
			Rnumu_bar[index] = 0. ;
			Rnux[index] = 0. ;

	                eta_nue[index] = 0. ;
	                eta_nue_bar[index] = 0. ;
	                eta_numu[index] = 0. ;
	                eta_numu_bar[index] = 0. ;
	                eta_nux[index] = 0. ;
			
			continue ;
		
		}
		
		CCTK_REAL _T,_ye;
		CCTK_REAL _ymu;
		
		std::array<CCTK_REAL,REQ_FLUID_VARS> U;
		U[RHO]  = rho_b[index];
		U[YE]   = ye[index];
		U[YMU]  = ymu[index];
		U[TEMP] = temp[index];

		_T = U[TEMP] ; _ye = U[YE] ; _ymu = U[YMU] ; 
		
		//The EAS class is provided by margherita with all the
		//functions and members called below. The EAS class is
		//basically the same as the Leakage class, but with some
		//modifications so that we get the correct opacities/emissivities
		//for M1.
		const auto tau_init = EAS<double>::compute_analytic_opacity(U[RHO]);
		std::array<double,5> tau_n {{tau_init, tau_init, tau_init, tau_init, tau_init}};
		if ( ! CCTK_Equals(m1_which_tau, "simple") ) {
		tau_n[NUE]       = tau_nue[index];
		tau_n[NUE_BAR]   = tau_nue_bar[index];
		tau_n[NUMU]      = tau_numu[index];
		tau_n[NUMU_BAR]  = tau_numu_bar[index];
		tau_n[NUX]       = tau_nux[index];
		}

		auto F= Fugacities<CCTK_REAL>(U[RHO], U[TEMP], U[YE], U[YMU], tau_n);
		std::array<bool,M1_interactions::NINTERACTIONS> which_interactions { { beta_decay,
											plasmon_decay,
											bremstrahlung,
											pair_annihil }} ;
		
		EAS<CCTK_REAL> eas{F, which_interactions};
		eas.calc_eas();

        if (use_grrt_beta_eq) {

			tau_betaeq_nue[index]       = 1./sqrt(eas.kappa_a_cactus<NUE>()*(eas.kappa_a_cactus<NUE>() + eas.kappa_s_cactus<NUE>() ) + 1.e-45) ;
			tau_betaeq_nue_bar[index]   = 1./sqrt(eas.kappa_a_cactus<NUE_BAR>()*(eas.kappa_a_cactus<NUE_BAR>() + eas.kappa_s_cactus<NUE_BAR>() ) + 1.e-45) ;
			tau_betaeq_numu[index]      = 1./sqrt(eas.kappa_a_cactus<NUMU>()*(eas.kappa_a_cactus<NUMU>() + eas.kappa_s_cactus<NUMU>() ) + 1.e-45) ;
			tau_betaeq_numu_bar[index]  = 1./sqrt(eas.kappa_a_cactus<NUMU_BAR>()*(eas.kappa_a_cactus<NUMU_BAR>() + eas.kappa_s_cactus<NUMU_BAR>() ) + 1.e-45) ;
			//tau_betaeq_nux[index]       = 1./sqrt(eas.kappa_a_cactus<NUX>()*(eas.kappa_a_cactus<NUX>() + eas.kappa_s_cactus<NUX>() ) + 1.e-45) ;
	
			// get the minimum beta-equil timescale across the three species 
			CCTK_REAL tau_betaeq = std::min(tau_betaeq_nue[index], std::min(tau_betaeq_nue_bar[index], std::min(tau_betaeq_numu[index], tau_betaeq_numu_bar[index] ))) ; 
			CCTK_REAL dt = CCTK_DELTA_TIME ;
			CCTK_REAL beta_equil_tscale_ = std::isnormal(tau_betaeq) ? (tau_betaeq / dt) : 0.0 ;
	
			if(beta_equil_tscale_ < 1.) {
			nfix++ ;
			
			CCTK_REAL ye_eq, T_eq;
			CCTK_REAL ymu_eq;
			bool success ;
			unsigned int niter;
			
			metric_c gamma( {gxx[index], gxy[index], gxz[index], gyy[index], gyz[index], gzz[index]},
			    {betax[index],betay[index],betaz[index]}, alp[index] ) ;
		
	//Harry: becaseful here!	
			std::array<CCTK_REAL,2*NUMSPECIES> Urad;

			auto const sqrtgamma = gamma.sqrtdet ;
	
		    closure_t cl( {zvecx[index],zvecy[index],zvecz[index]}, &gamma) ;

		    cl.update_closure<m1_closure_f>(Enue[index], {Fnue_x[index], Fnue_y[index], Fnue_z[index]}, zeta_nue[index],false) ; 
		    Urad[0] = cl.get_J() ;
		    
		    cl.update_closure<m1_closure_f>(Enue_bar[index], {Fnue_bar_x[index], Fnue_bar_y[index], Fnue_bar_z[index]}, zeta_nue_bar[index],false) ; 
		    Urad[1] = cl.get_J() ;

		    cl.update_closure<m1_closure_f>(Enumu[index], {Fnumu_x[index], Fnumu_y[index], Fnumu_z[index]}, zeta_numu[index],false) ; 
		    Urad[2] = cl.get_J() ;

		    cl.update_closure<m1_closure_f>(Enumu_bar[index], {Fnumu_bar_x[index], Fnumu_bar_y[index], Fnumu_bar_z[index]}, zeta_numu_bar[index],false) ; 
		    Urad[3] = cl.get_J() ;

		    cl.update_closure<m1_closure_f>(Enux[index], {Fnux_x[index], Fnux_y[index], Fnux_z[index]}, zeta_nux[index],false) ; 
		    Urad[4] = cl.get_J() ;


			//Urad[0] = Enue[index]; Urad[1] = Enue_bar[index] ; Urad[2] = Enumu[index] ; Urad[3] = Enumu_bar[index] ; Urad[4] = Enux[index] ;
			Urad[5] = Nnue[index] ; Urad[6] = Nnue_bar[index] ; Urad[7] = Nnumu[index] ; Urad[8] = Nnumu_bar[index] ; Urad[9] = Nnux[index] ;
	//Harry: change compute_T_ye_betaeq reading Urad <10> and tau_n <5> and read ymu, change to Leptonic
			std::tie(ye_eq,ymu_eq,T_eq,success,niter) = compute_T_ye_ymu_betaeq<EOS_Leptonic>(Urad,U, tau_n, which_interactions) ;
			
				if ( success ) {
					nsucc++;
					avgiter += niter ;
					if( beta_equil_tscale_ < 0.5 ) {
						_T = T_eq;
						_ye = ye_eq;
						_ymu = ymu_eq;
						error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) + fabs(_ymu-U[YMU]);
						error_int_denom += _ye + _T + _ymu;
					} else {
						CCTK_REAL fac = 2 * ( 1.0 - beta_equil_tscale_ ) ;
						_T = fac * T_eq  + ( 1 - fac ) * U[TEMP] ;
						_ye = fac * ye_eq + ( 1 - fac ) * U[YE] ;
						_ymu = fac * ymu_eq + ( 1 - fac ) * U[YMU] ;
						error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) + fabs(_ymu-U[YMU]);
						error_int_denom += _ye + _T + _ymu;
					}
	//Harry:	: add ymu
					F = Fugacities<CCTK_REAL>(U[RHO],_T, _ye, _ymu, tau_n) ;
					eas.set_F(F) ;
					eas.calc_eas() ;
					
				} // if success
			}
		} // end if of grrt_beta
			
		temp_nue[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue[index], F.eta[NUE]);
		temp_nue_bar[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue_bar[index], F.eta[NUE_BAR]);
		temp_numu[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_numu[index], F.eta[NUMU]);
		temp_numu_bar[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_numu_bar[index], F.eta[NUMU_BAR]);
		temp_nux[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nux[index], F.eta[NUX]);

		CCTK_REAL const nue_fact = ::max( 1., ::int_pow<2>(temp_nue[index]/temp[index]));
		CCTK_REAL const nue_bar_fact = ::max( 1., ::int_pow<2>(temp_nue_bar[index]/temp[index]));
		CCTK_REAL const numu_fact = ::max( 1., ::int_pow<2>(temp_numu[index]/temp[index]));
		CCTK_REAL const numu_bar_fact = ::max( 1., ::int_pow<2>(temp_numu_bar[index]/temp[index]));
		CCTK_REAL const nux_fact = ::max( 1., ::int_pow<2>(temp_nux[index]/temp[index]));

		Qnue[index] = eas.Q_cactus<NUE>() * nue_fact ;
		Qnue_bar[index] = eas.Q_cactus<NUE_BAR>() * nue_bar_fact ; 
		Qnumu[index] = eas.Q_cactus<NUMU>() * numu_fact ;
		Qnumu_bar[index] = eas.Q_cactus<NUMU_BAR>() * numu_bar_fact ; 
		Qnux[index] = eas.Q_cactus<NUX>() ;

		Rnue[index] = eas.R_cactus<NUE>() * nue_fact ; 
		Rnue_bar[index] = eas.R_cactus<NUE_BAR>() * nue_bar_fact ;
		Rnumu[index] = eas.R_cactus<NUMU>() * numu_fact ; 
		Rnumu_bar[index] = eas.R_cactus<NUMU_BAR>() * numu_bar_fact ;
		Rnux[index] = eas.R_cactus<NUX>() ;

		kappa_nue_a[index] = eas.kappa_a_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_a[index] = eas.kappa_a_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_a[index] = eas.kappa_a_cactus<NUMU>() * numu_fact ;
		kappa_numu_bar_a[index] = eas.kappa_a_cactus<NUMU_BAR>() * numu_bar_fact ;
		kappa_nux_a[index] = eas.kappa_a_cactus<NUX>() ;

		kappa_nue_s[index] = eas.kappa_s_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_s[index] = eas.kappa_s_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_s[index] = eas.kappa_s_cactus<NUMU>() * numu_fact ;
		kappa_numu_bar_s[index] = eas.kappa_s_cactus<NUMU_BAR>() * numu_bar_fact ;
		kappa_nux_s[index] = eas.kappa_s_cactus<NUX>() * nux_fact;

		kappa_nue_n[index] = eas.kappa_n_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_n[index] = eas.kappa_n_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_n[index] = eas.kappa_n_cactus<NUMU>() * numu_fact ;
		kappa_numu_bar_n[index] = eas.kappa_n_cactus<NUMU_BAR>() * numu_bar_fact ;
		kappa_nux_n[index] = eas.kappa_n_cactus<NUX>() ;

		
		eta_nue[index] = F.eta[NUE] ;
		eta_nue_bar[index] = F.eta[NUE_BAR] ;
		eta_numu[index] = F.eta[NUMU] ;
		eta_numu_bar[index] = F.eta[NUMU_BAR] ;
		eta_nux[index] = F.eta[NUX] ;

           } else {

		if( rho_b[index] < M1_rho_floor * ( 1.0 + M1_atmo_tol ) ) {

			Qnue[index] = 0. ;
			Qnue_bar[index] = 0. ;
			Qnux[index] = 0. ;

			kappa_nue_a[index] = 0.;
			kappa_nue_bar_a[index] = 0.;
			kappa_nux_a[index] = 0.;

			kappa_nue_s[index] = 0.;
			kappa_nue_bar_s[index] = 0.;
			kappa_nux_s[index] = 0.;

			kappa_nue_n[index] = 0.;
			kappa_nue_bar_n[index] = 0.;
			kappa_nux_n[index] = 0.;

			Rnue[index] = 0. ;
			Rnue_bar[index] = 0. ;
			Rnux[index] = 0. ;
			
	                eta_nue[index] = 0. ;
	                eta_nue_bar[index] = 0. ;
	                eta_nux[index] = 0. ;

			continue ;
		
		}
		
		CCTK_REAL _T,_ye,_ymu;
		
		std::array<CCTK_REAL,REQ_FLUID_VARS> U;
		U[RHO]  = rho_b[index];
		U[YE]   = ye[index];
		U[YMU]  = ymu[index];
		U[TEMP] = temp[index];

		_T = U[TEMP] ; _ye = U[YE] ; _ymu = U[YMU] ;
		
		//The EAS class is provided by margherita with all the
		//functions and members called below. The EAS class is
		//basically the same as the Leakage class, but with some
		//modifications so that we get the correct opacities/emissivities
		//for M1.
		const auto tau_init = EAS<double>::compute_analytic_opacity(U[RHO]);
		std::array<double,5> tau_n {{tau_init, tau_init, tau_init, tau_init, tau_init}};
		if ( ! CCTK_Equals(m1_which_tau, "simple") ) {
		tau_n[NUE]      = tau_nue[index];
		tau_n[NUE_BAR]  = tau_nue_bar[index];
		tau_n[NUMU]     = 0.0 ; 
		tau_n[NUMU_BAR] = 0.0 ;
		tau_n[NUX]      = tau_nux[index];
		}

		auto F= Fugacities<CCTK_REAL>(U[RHO], U[TEMP], U[YE], U[YMU], tau_n);
		std::array<bool,M1_interactions::NINTERACTIONS> which_interactions { { beta_decay,
											plasmon_decay,
											bremstrahlung,
											pair_annihil }} ;
		
		EAS<CCTK_REAL> eas{F, which_interactions};
		eas.calc_eas();

                if (use_grrt_beta_eq) {

			tau_betaeq_nue[index]      = 1./sqrt(eas.kappa_a_cactus<NUE>()*(eas.kappa_a_cactus<NUE>() + eas.kappa_s_cactus<NUE>() ) + 1.e-45) ;
			tau_betaeq_nue_bar[index]  = 1./sqrt(eas.kappa_a_cactus<NUE_BAR>()*(eas.kappa_a_cactus<NUE_BAR>() + eas.kappa_s_cactus<NUE_BAR>() ) + 1.e-45) ;
			//tau_betaeq_nux      = 1./sqrt(eas.kappa_a_cactus<NUX>()*(eas.kappa_a_cactus<NUX>() + eas.kappa_s_cactus<NUX>() ) + 1.e-45) ;
	
			// get the minimum beta-equil timescale across the three species 
			CCTK_REAL tau_betaeq = std::min(tau_betaeq_nue[index], tau_betaeq_nue_bar[index]) ; 
			CCTK_REAL dt = CCTK_DELTA_TIME ;
			CCTK_REAL beta_equil_tscale_ = std::isnormal(tau_betaeq) ? (tau_betaeq / dt) : 0.0 ;
	
			
			if(beta_equil_tscale_ < 1.) {
			nfix++ ;
			
			CCTK_REAL ye_eq, T_eq, ymu_eq;
			bool success ;
			unsigned int niter;
			
	                std::array<CCTK_REAL,2*NUMSPECIES> Urad;
	                Urad[0] = Enue[index]; Urad[1] = Enue_bar[index] ; Urad[2] = 0.0 ; Urad[3] = 0.0 ; Urad[4] = Enux[index] ;
	                Urad[5] = Nnue[index] ; Urad[6] = Nnue_bar[index] ; Urad[7] = 0.0 ; Urad[8] = 0.0 ; Urad[9] = Nnux[index] ;
	
			std::tie(ye_eq,ymu_eq,T_eq,success,niter) = compute_T_ye_ymu_betaeq<EOS_Tabulated>(Urad, U, tau_n, which_interactions) ;
			
			if ( success ) {
				nsucc++;
				avgiter += niter ;
				if( beta_equil_tscale_ < 0.5 ) {
					_T = T_eq;
					_ye = ye_eq;
					_ymu = ymu_eq;
					error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) + fabs(_ymu-U[YMU]);
					error_int_denom += _ye + _T + _ymu;
				} else {
					CCTK_REAL fac = 2 * ( 1.0 - beta_equil_tscale_ ) ;
					_T = fac * T_eq  + ( 1 - fac ) * U[TEMP] ;
					_ye = fac * ye_eq + ( 1 - fac ) * U[YE] ;
					_ymu = fac * ymu_eq + ( 1 - fac ) * U[YMU] ;
					error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) + fabs(_ymu-U[YMU]);
					error_int_denom += _ye + _T + _ymu;
				}
				F = Fugacities<CCTK_REAL>(U[RHO],_T, _ye, _ymu, tau_n) ;
				eas.set_F(F) ;
				eas.calc_eas() ;
				
			} // if success
			}
		} //endif of grrt_beta
			
		temp_nue[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue[index], F.eta[NUE]);
		temp_nue_bar[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue_bar[index], F.eta[NUE_BAR]);
		temp_numu[index] = 0.0 ; 
		temp_numu_bar[index] = 0.0 ;
		temp_nux[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nux[index], F.eta[NUX]);

		CCTK_REAL const nue_fact = ::max( 1., ::int_pow<2>(temp_nue[index]/temp[index]));
		CCTK_REAL const nue_bar_fact = ::max( 1., ::int_pow<2>(temp_nue_bar[index]/temp[index]));
		CCTK_REAL const nux_fact = ::max( 1., ::int_pow<2>(temp_nux[index]/temp[index]));

		Qnue[index] = eas.Q_cactus<NUE>() * nue_fact ;
		Qnue_bar[index] = eas.Q_cactus<NUE_BAR>() * nue_bar_fact ; 
		Qnumu[index] = 0.0 ;
		Qnumu_bar[index] = 0.0 ;
		Qnux[index] = eas.Q_cactus<NUX>() ;

		Rnue[index] = eas.R_cactus<NUE>() * nue_fact ; 
		Rnue_bar[index] = eas.R_cactus<NUE_BAR>() * nue_bar_fact ;
		Rnumu[index] = 0.0 ;
		Rnumu_bar[index] = 0.0 ;
		Rnux[index] = eas.R_cactus<NUX>() ;

		kappa_nue_a[index] = eas.kappa_a_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_a[index] = eas.kappa_a_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_a[index] = 0.0 ;
		kappa_numu_bar_a[index] = 0.0 ;
		kappa_nux_a[index] = eas.kappa_a_cactus<NUX>() ;

		kappa_nue_s[index] = eas.kappa_s_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_s[index] = eas.kappa_s_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_s[index] = 0.0 ;
		kappa_numu_bar_s[index] = 0.0 ;
		kappa_nux_s[index] = eas.kappa_s_cactus<NUX>() * nux_fact ;

		kappa_nue_n[index] = eas.kappa_n_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_n[index] = eas.kappa_n_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_n[index] = 0.0 ;
		kappa_numu_bar_n[index] = 0.0 ;
		kappa_nux_n[index] = eas.kappa_n_cactus<NUX>() ;
		
		eta_nue[index] = F.eta[NUE] ;
		eta_nue_bar[index] = F.eta[NUE_BAR] ;
		eta_numu[index] = 0.0 ;
		eta_numu_bar[index] = 0.0 ;
		eta_nux[index] = F.eta[NUX] ;

                // using 3 spec, numu numubar zeros:

                Qnumu[index] = 0. ;
                Qnumu_bar[index] = 0. ;

                kappa_numu_a[index] = 0.;
                kappa_numu_bar_a[index] = 0.;

                kappa_numu_s[index] = 0.;
                kappa_numu_bar_s[index] = 0.;

                kappa_numu_n[index] = 0.;
                kappa_numu_bar_n[index] = 0.;

                Rnumu[index] = 0. ;
                Rnumu_bar[index] = 0. ;

		eta_numu[index] = 0.0 ;
		eta_numu_bar[index] = 0.0 ;

	   }  // use_5_spec_m1
	} //for loop


	if(verbosity>1)
		CCTK_VInfo(CCTK_THORNSTRING,
			"***** EAS: Lev: %d |  Beta equilibrium fixes: applied:success:total %d:%d:%d \t | avg iter/point %1.4f | relative error %1.4g denominator %1.4g *****\n",
			(int)GetRefinementLevel(cctkGH),
			nfix,nsucc,npoints,avgiter/npoints,
			error_int/(error_int_denom+1.e-40),error_int_denom);
  
} // FIL_M1_compute_eas_corrected


void frankfurt_m1_compute_eas_corrected_Weakhub(CCTK_ARGUMENTS) {
	DECLARE_CCTK_ARGUMENTS;
	DECLARE_CCTK_PARAMETERS;

	using namespace M1_Constants;
	using namespace Margherita_constants ;
	using namespace Margherita_M1_EAS ;
//code test
  using namespace frankfurt_m1 ;


	//Global quantities for computing avergae over entire grid
	
	CCTK_INT nfix=0, nsucc=0;
	CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2] ;
	CCTK_REAL avgiter=0., error_int=0., error_int_denom=0.  ;

	
	#pragma omp parallel for collapse(3) reduction(+:nfix,nsucc,avgiter,error_int,error_int_denom) schedule(static)
	for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {

		int index = CCTK_GFINDEX3D(cctkGH,i,j,k);


// code test
//                        Qnue[index] = 1e-45 ;
//                        Qnue_bar[index] = 1e-45 ;
//                        Qnumu[index] = 1e-45 ;
//                        Qnumu_bar[index] = 1e-45 ;
//                        Qnux[index] = 1e-45 ;
//
//                        kappa_nue_a[index] = 1e-45;
//                        kappa_nue_bar_a[index] = 1e-45;
//                        kappa_numu_a[index] = 1e-45;
//                        kappa_numu_bar_a[index] = 1e-45;
//                        kappa_nux_a[index] = 1e-45;
//
//                        kappa_nue_s[index] = 1e-45;
//                        kappa_nue_bar_s[index] = 1e-45;
//                        kappa_numu_s[index] = 1e-45;
//                        kappa_numu_bar_s[index] = 1e-45;
//                        kappa_nux_s[index] = 1e-45;
//
//                        kappa_nue_n[index] = 1e-45;
//                        kappa_nue_bar_n[index] = 1e-45;
//                        kappa_numu_n[index] = 1e-45;
//                        kappa_numu_bar_n[index] = 1e-45;
//                        kappa_nux_n[index] = 1e-45;
//
//                        Rnue[index] = 1e-45 ;
//                        Rnue_bar[index] = 1e-45 ;
//                        Rnumu[index] = 1e-45 ;
//                        Rnumu_bar[index] = 1e-45 ;
//                        Rnux[index] = 1e-45 ;
//			continue ;



           if (use_5_spec_m1) {

		if( rho_b[index] < M1_rho_floor * ( 1.0 + M1_atmo_tol ) ) {

			Qnue[index] = 0. ;
			Qnue_bar[index] = 0. ;
			Qnumu[index] = 0. ;
			Qnumu_bar[index] = 0. ;
			Qnux[index] = 0. ;

			kappa_nue_a[index] = 0.;
			kappa_nue_bar_a[index] = 0.;
			kappa_numu_a[index] = 0.;
			kappa_numu_bar_a[index] = 0.;
			kappa_nux_a[index] = 0.;

			kappa_nue_s[index] = 0.;
			kappa_nue_bar_s[index] = 0.;
			kappa_numu_s[index] = 0.;
			kappa_numu_bar_s[index] = 0.;
			kappa_nux_s[index] = 0.;

			kappa_nue_n[index] = 0.;
			kappa_nue_bar_n[index] = 0.;
			kappa_numu_n[index] = 0.;
			kappa_numu_bar_n[index] = 0.;
			kappa_nux_n[index] = 0.;

			Rnue[index] = 0. ;
			Rnue_bar[index] = 0. ;
			Rnumu[index] = 0. ;
			Rnumu_bar[index] = 0. ;
			Rnux[index] = 0. ;
			
	                eta_nue[index] = 0. ;
	                eta_nue_bar[index] = 0. ;
	                eta_numu[index] = 0. ;
	                eta_numu_bar[index] = 0. ;
	                eta_nux[index] = 0. ;

			continue ;
		
		}
		
		CCTK_REAL _T,_ye;
		CCTK_REAL _ymu;
		
		std::array<CCTK_REAL,REQ_FLUID_VARS> U;
		U[RHO]  = rho_b[index];
		U[YE]   = ye[index];
		U[YMU]  = ymu[index];
		U[TEMP] = temp[index];

		_T = U[TEMP] ; _ye = U[YE] ; _ymu = U[YMU] ; 
		
		//The EAS class is provided by margherita with all the
		//functions and members called below. The EAS class is
		//basically the same as the Leakage class, but with some
		//modifications so that we get the correct opacities/emissivities
		//for M1.
		const auto tau_init = EAS<double>::compute_analytic_opacity(U[RHO]);
		std::array<double,5> tau_n {{tau_init, tau_init, tau_init, tau_init, tau_init}};
		if ( ! CCTK_Equals(m1_which_tau, "simple") ) {
		tau_n[NUE]       = tau_nue[index];
		tau_n[NUE_BAR]   = tau_nue_bar[index];
		tau_n[NUMU]      = tau_numu[index];
		tau_n[NUMU_BAR]  = tau_numu_bar[index];
		tau_n[NUX]       = tau_nux[index];
		}

		auto F= Fugacities<CCTK_REAL>(U[RHO], U[TEMP], U[YE], U[YMU], tau_n);
		std::array<bool,M1_interactions::NINTERACTIONS> which_interactions { { beta_decay,
											plasmon_decay,
											bremstrahlung,
											pair_annihil }} ;
		
		EAS<CCTK_REAL> eas{F, which_interactions};
		eas.calc_eas();

                if (use_grrt_beta_eq) {
	
	
			tau_betaeq_nue[index]       = 1./sqrt(eas.kappa_a_cactus<NUE>()*(eas.kappa_a_cactus<NUE>() + eas.kappa_s_cactus<NUE>() ) + 1.e-45) ;
			tau_betaeq_nue_bar[index]   = 1./sqrt(eas.kappa_a_cactus<NUE_BAR>()*(eas.kappa_a_cactus<NUE_BAR>() + eas.kappa_s_cactus<NUE_BAR>() ) + 1.e-45) ;
			tau_betaeq_numu[index]      = 1./sqrt(eas.kappa_a_cactus<NUMU>()*(eas.kappa_a_cactus<NUMU>() + eas.kappa_s_cactus<NUMU>() ) + 1.e-45) ;
			tau_betaeq_numu_bar[index]  = 1./sqrt(eas.kappa_a_cactus<NUMU_BAR>()*(eas.kappa_a_cactus<NUMU_BAR>() + eas.kappa_s_cactus<NUMU_BAR>() ) + 1.e-45) ;
			//tau_betaeq_nux       = 1./sqrt(eas.kappa_a_cactus<NUX>()*(eas.kappa_a_cactus<NUX>() + eas.kappa_s_cactus<NUX>() ) + 1.e-45) ;
	
			// get the minimum beta-equil timescale across the three species 
			CCTK_REAL tau_betaeq = std::min(tau_betaeq_nue[index], std::min(tau_betaeq_nue_bar[index], std::min(tau_betaeq_numu[index], tau_betaeq_numu_bar[index]))) ; 
			CCTK_REAL dt = CCTK_DELTA_TIME ;
			CCTK_REAL beta_equil_tscale_ = std::isnormal(tau_betaeq) ? (tau_betaeq / dt) : 0.0 ;
	
			
			if(beta_equil_tscale_ < 1.) {
			nfix++ ;
			
			CCTK_REAL ye_eq, T_eq;
			CCTK_REAL ymu_eq;
			bool success ;
			unsigned int niter;
		
	//Harry: becaseful here!	
			std::array<CCTK_REAL,2*NUMSPECIES> Urad;
			Urad[0] = Enue[index]; Urad[1] = Enue_bar[index] ; Urad[2] = Enumu[index] ; Urad[3] = Enumu_bar[index] ; Urad[4] = Enux[index] ;
			Urad[5] = Nnue[index] ; Urad[6] = Nnue_bar[index] ; Urad[7] = Nnumu[index] ; Urad[8] = Nnumu_bar[index] ; Urad[9] = Nnux[index] ;
	//Harry: change compute_T_ye_betaeq reading Urad <10> and tau_n <5> and read ymu, change to Leptonic
			std::tie(ye_eq,ymu_eq,T_eq,success,niter) = compute_T_ye_ymu_betaeq<EOS_Leptonic>(Urad,U, tau_n, which_interactions) ;
			
				if ( success ) {
					nsucc++;
					avgiter += niter ;
					if( beta_equil_tscale_ < 0.5 ) {
						_T = T_eq;
						_ye = ye_eq;
						_ymu = ymu_eq;
						error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) + fabs(_ymu-U[YMU]);
						error_int_denom += _ye + _T + _ymu;
					} else {
						CCTK_REAL fac = 2 * ( 1.0 - beta_equil_tscale_ ) ;
						_T = fac * T_eq  + ( 1 - fac ) * U[TEMP] ;
						_ye = fac * ye_eq + ( 1 - fac ) * U[YE] ;
						_ymu = fac * ymu_eq + ( 1 - fac ) * U[YMU] ;
						error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) + fabs(_ymu-U[YMU]);
						error_int_denom += _ye + _T + _ymu;
					}
	//Harry:	: add ymu
					F = Fugacities<CCTK_REAL>(U[RHO],_T, _ye, _ymu, tau_n) ;
					eas.set_F(F) ;
					eas.calc_eas() ;
					
				} // if success
			}
		} // use_grrt_beta_eq
//code test
//        metric_c gamma( {gxx[index], gxy[index], gxz[index], gyy[index], gyz[index], gzz[index]},
//                        {betax[index],betay[index],betaz[index]},alp[index] ) ;
//
//        closure_t cl( {zvecx[index],zvecy[index],zvecz[index]}, &gamma ) ;
//
//                // unit correct???
//	 	eps_nue[index] = std::max(cl.get_average_energy(Nnue_star[index]) * Margherita_constants::mnuc_MeV,0.0) ; 
//	 	eps_nue_bar[index] = std::max(cl.get_average_energy(Nnue_bar_star[index]) * Margherita_constants::mnuc_MeV,0.0) ; 
//	 	eps_numu[index] = std::max(cl.get_average_energy(Nnumu_star[index]) * Margherita_constants::mnuc_MeV,0.0) ; 
//	 	eps_numu_bar[index] = std::max(cl.get_average_energy(Nnumu_bar_star[index]) * Margherita_constants::mnuc_MeV,0.0) ; 
//	 	eps_nux[index] = std::max(cl.get_average_energy(Nnux_star[index]) * Margherita_constants::mnuc_MeV,0.0) ; 

		temp_nue[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue[index], F.eta[NUE]);
		temp_nue_bar[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue_bar[index], F.eta[NUE_BAR]);
		temp_numu[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_numu[index], F.eta[NUMU]);
		temp_numu_bar[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_numu_bar[index], F.eta[NUMU_BAR]);
		temp_nux[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nux[index], F.eta[NUX]);

		CCTK_REAL const nue_fact = ::max( 1., ::int_pow<2>(temp_nue[index]/temp[index]));
		CCTK_REAL const nue_bar_fact = ::max( 1., ::int_pow<2>(temp_nue_bar[index]/temp[index]));
//code test
		CCTK_REAL const numu_fact = 1.0;
		CCTK_REAL const numu_bar_fact = 1.0;
		//CCTK_REAL const numu_fact = ::max( 1., ::int_pow<2>(temp_numu[index]/temp[index]));
		//CCTK_REAL const numu_bar_fact = ::max( 1., ::int_pow<2>(temp_numu_bar[index]/temp[index]));
		CCTK_REAL const nux_fact = ::max( 1., ::int_pow<2>(temp_nux[index]/temp[index]));


		Qnue[index] = eas.Q_cactus<NUE>() * nue_fact ;
		Qnue_bar[index] = eas.Q_cactus<NUE_BAR>() * nue_bar_fact ; 
		Qnumu[index] = eas.Q_cactus<NUMU>() * numu_fact ;
		Qnumu_bar[index] = eas.Q_cactus<NUMU_BAR>() * numu_bar_fact ; 
		Qnux[index] = eas.Q_cactus<NUX>() ;

		Rnue[index] = eas.R_cactus<NUE>() * nue_fact ; 
		Rnue_bar[index] = eas.R_cactus<NUE_BAR>() * nue_bar_fact ;
		Rnumu[index] = eas.R_cactus<NUMU>() * numu_fact ; 
		Rnumu_bar[index] = eas.R_cactus<NUMU_BAR>() * numu_bar_fact ;
		Rnux[index] = eas.R_cactus<NUX>() ;

		kappa_nue_a[index] = eas.kappa_a_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_a[index] = eas.kappa_a_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_a[index] = eas.kappa_a_cactus<NUMU>() * numu_fact ;
		kappa_numu_bar_a[index] = eas.kappa_a_cactus<NUMU_BAR>() * numu_bar_fact ;
		kappa_nux_a[index] = eas.kappa_a_cactus<NUX>() ;

		kappa_nue_s[index] = eas.kappa_s_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_s[index] = eas.kappa_s_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_s[index] = eas.kappa_s_cactus<NUMU>() * numu_fact ;
		kappa_numu_bar_s[index] = eas.kappa_s_cactus<NUMU_BAR>() * numu_bar_fact ;
		kappa_nux_s[index] = eas.kappa_s_cactus<NUX>() * nux_fact ;

		kappa_nue_n[index] = eas.kappa_n_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_n[index] = eas.kappa_n_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_n[index] = eas.kappa_n_cactus<NUMU>() * numu_fact ;
		kappa_numu_bar_n[index] = eas.kappa_n_cactus<NUMU_BAR>() * numu_bar_fact ;
		kappa_nux_n[index] = eas.kappa_n_cactus<NUX>() ;

		
		eta_nue[index] = F.eta[NUE] ;
		eta_nue_bar[index] = F.eta[NUE_BAR] ;
		eta_numu[index] = F.eta[NUMU] ;
		eta_numu_bar[index] = F.eta[NUMU_BAR] ;
		eta_nux[index] = F.eta[NUX] ;

           } else {

		if( rho_b[index] < M1_rho_floor * ( 1.0 + M1_atmo_tol ) ) {

			Qnue[index] = 0. ;
			Qnue_bar[index] = 0. ;
			Qnux[index] = 0. ;

			kappa_nue_a[index] = 0.;
			kappa_nue_bar_a[index] = 0.;
			kappa_nux_a[index] = 0.;

			kappa_nue_s[index] = 0.;
			kappa_nue_bar_s[index] = 0.;
			kappa_nux_s[index] = 0.;

			kappa_nue_n[index] = 0.;
			kappa_nue_bar_n[index] = 0.;
			kappa_nux_n[index] = 0.;

			Rnue[index] = 0. ;
			Rnue_bar[index] = 0. ;
			Rnux[index] = 0. ;
			
	                eta_nue[index] = 0. ;
	                eta_nue_bar[index] = 0. ;
	                eta_nux[index] = 0. ;
			continue ;
		
		}
		
		CCTK_REAL _T,_ye,_ymu;
		
		std::array<CCTK_REAL,REQ_FLUID_VARS> U;
		U[RHO]  = rho_b[index];
		U[YE]   = ye[index];
		U[YMU]  = ymu[index];
		U[TEMP] = temp[index];

		_T = U[TEMP] ; _ye = U[YE] ; _ymu = U[YMU] ;
		
		//The EAS class is provided by margherita with all the
		//functions and members called below. The EAS class is
		//basically the same as the Leakage class, but with some
		//modifications so that we get the correct opacities/emissivities
		//for M1.
		const auto tau_init = EAS<double>::compute_analytic_opacity(U[RHO]);
		std::array<double,5> tau_n {{tau_init, tau_init, tau_init, tau_init, tau_init}};
		if ( ! CCTK_Equals(m1_which_tau, "simple") ) {
		tau_n[NUE]      = tau_nue[index];
		tau_n[NUE_BAR]  = tau_nue_bar[index];
		tau_n[NUMU]     = 0.0 ; 
		tau_n[NUMU_BAR] = 0.0 ;
		tau_n[NUX]      = tau_nux[index];
		}

		auto F= Fugacities<CCTK_REAL>(U[RHO], U[TEMP], U[YE], U[YMU], tau_n);
		std::array<bool,M1_interactions::NINTERACTIONS> which_interactions { { beta_decay,
											plasmon_decay,
											bremstrahlung,
											pair_annihil }} ;
		
		EAS<CCTK_REAL> eas{F, which_interactions};
		eas.calc_eas();

                if (use_grrt_beta_eq) {
			tau_betaeq_nue[index]      = 1./sqrt(eas.kappa_a_cactus<NUE>()*(eas.kappa_a_cactus<NUE>() + eas.kappa_s_cactus<NUE>() ) + 1.e-45) ;
			tau_betaeq_nue_bar[index]  = 1./sqrt(eas.kappa_a_cactus<NUE_BAR>()*(eas.kappa_a_cactus<NUE_BAR>() + eas.kappa_s_cactus<NUE_BAR>() ) + 1.e-45) ;
			//tau_betaeq_nux      = 1./sqrt(eas.kappa_a_cactus<NUX>()*(eas.kappa_a_cactus<NUX>() + eas.kappa_s_cactus<NUX>() ) + 1.e-45) ;
	
			// get the minimum beta-equil timescale across the three species 
			CCTK_REAL tau_betaeq = std::min(tau_betaeq_nue[index], tau_betaeq_nue_bar[index]) ; 
			CCTK_REAL dt = CCTK_DELTA_TIME ;
			CCTK_REAL beta_equil_tscale_ = std::isnormal(tau_betaeq) ? (tau_betaeq / dt) : 0.0 ;
	
			
			if(beta_equil_tscale_ < 1.) {
			nfix++ ;
			
			CCTK_REAL ye_eq, T_eq, ymu_eq;
			bool success ;
			unsigned int niter;
			
	                std::array<CCTK_REAL,2*NUMSPECIES> Urad;
	                Urad[0] = Enue[index]; Urad[1] = Enue_bar[index] ; Urad[2] = 0.0 ; Urad[3] = 0.0 ; Urad[4] = Enux[index] ;
	                Urad[5] = Nnue[index] ; Urad[6] = Nnue_bar[index] ; Urad[7] = 0.0 ; Urad[8] = 0.0 ; Urad[9] = Nnux[index] ;
	
			std::tie(ye_eq,ymu_eq,T_eq,success,niter) = compute_T_ye_ymu_betaeq<EOS_Tabulated>(Urad, U, tau_n, which_interactions) ;
			
			if ( success ) {
				nsucc++;
				avgiter += niter ;
				if( beta_equil_tscale_ < 0.5 ) {
					_T = T_eq;
					_ye = ye_eq;
					_ymu = ymu_eq;
					error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) + fabs(_ymu-U[YMU]);
					error_int_denom += _ye + _T + _ymu;
				} else {
					CCTK_REAL fac = 2 * ( 1.0 - beta_equil_tscale_ ) ;
					_T = fac * T_eq  + ( 1 - fac ) * U[TEMP] ;
					_ye = fac * ye_eq + ( 1 - fac ) * U[YE] ;
					_ymu = fac * ymu_eq + ( 1 - fac ) * U[YMU] ;
					error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) + fabs(_ymu-U[YMU]);
					error_int_denom += _ye + _T + _ymu;
				}
				F = Fugacities<CCTK_REAL>(U[RHO],_T, _ye, _ymu, tau_n) ;
				eas.set_F(F) ;
				eas.calc_eas() ;
				
			} // if success
			} // beta_equil_tscale_ < 1.
		} // use_grrt_beta_eq

		temp_nue[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue[index], F.eta[NUE]);
		temp_nue_bar[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue_bar[index], F.eta[NUE_BAR]);
		temp_numu[index] = 0.0 ; 
		temp_numu_bar[index] = 0.0 ;
		temp_nux[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nux[index], F.eta[NUX]);

		CCTK_REAL const nue_fact = ::max( 1., ::int_pow<2>(temp_nue[index]/temp[index]));
		CCTK_REAL const nue_bar_fact = ::max( 1., ::int_pow<2>(temp_nue_bar[index]/temp[index]));
		CCTK_REAL const nux_fact = ::max( 1., ::int_pow<2>(temp_nux[index]/temp[index]));

		Qnue[index] = eas.Q_cactus<NUE>() * nue_fact ;
		Qnue_bar[index] = eas.Q_cactus<NUE_BAR>() * nue_bar_fact ; 
		Qnumu[index] = 0.0 ;
		Qnumu_bar[index] = 0.0 ;
		Qnux[index] = eas.Q_cactus<NUX>() ;

		Rnue[index] = eas.R_cactus<NUE>() * nue_fact ; 
		Rnue_bar[index] = eas.R_cactus<NUE_BAR>() * nue_bar_fact ;
		Rnumu[index] = 0.0 ;
		Rnumu_bar[index] = 0.0 ;
		Rnux[index] = eas.R_cactus<NUX>() ;

		kappa_nue_a[index] = eas.kappa_a_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_a[index] = eas.kappa_a_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_a[index] = 0.0 ;
		kappa_numu_bar_a[index] = 0.0 ;
		kappa_nux_a[index] = eas.kappa_a_cactus<NUX>() ;

		kappa_nue_s[index] = eas.kappa_s_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_s[index] = eas.kappa_s_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_s[index] = 0.0 ;
		kappa_numu_bar_s[index] = 0.0 ;
		kappa_nux_s[index] = eas.kappa_s_cactus<NUX>() * nux_fact ;

		kappa_nue_n[index] = eas.kappa_n_cactus<NUE>() * nue_fact ;
		kappa_nue_bar_n[index] = eas.kappa_n_cactus<NUE_BAR>() * nue_bar_fact ;
		kappa_numu_n[index] = 0.0 ;
		kappa_numu_bar_n[index] = 0.0 ;
		kappa_nux_n[index] = eas.kappa_n_cactus<NUX>() ;
		
		eta_nue[index] = F.eta[NUE] ;
		eta_nue_bar[index] = F.eta[NUE_BAR] ;
		eta_numu[index] = 0.0 ;
		eta_numu_bar[index] = 0.0 ;
		eta_nux[index] = F.eta[NUX] ;

                // using 3 spec, numu numubar zeros:

                Qnumu[index] = 0. ;
                Qnumu_bar[index] = 0. ;

                kappa_numu_a[index] = 0.;
                kappa_numu_bar_a[index] = 0.;

                kappa_numu_s[index] = 0.;
                kappa_numu_bar_s[index] = 0.;

                kappa_numu_n[index] = 0.;
                kappa_numu_bar_n[index] = 0.;

                Rnumu[index] = 0. ;
                Rnumu_bar[index] = 0. ;

		eta_numu[index] = 0.0 ;
		eta_numu_bar[index] = 0.0 ;

	   }  // use_5_spec_m1
	} //for loop


	if(verbosity>1)
		CCTK_VInfo(CCTK_THORNSTRING,
			"***** EAS: Lev: %d |  Beta equilibrium fixes: applied:success:total %d:%d:%d \t | avg iter/point %1.4f | relative error %1.4g denominator %1.4g *****\n",
			(int)GetRefinementLevel(cctkGH),
			nfix,nsucc,npoints,avgiter/npoints,
			error_int/(error_int_denom+1.e-40),error_int_denom);
  
} // FIL_M1_compute_eas_Weakhub



