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
#include "Margherita_M1.h"
#include "utils/bb.hh"
//#define DEBUG_FILM1

void frankfurt_m1_compute_eas_microphysical(CCTK_ARGUMENTS);
void frankfurt_m1_compute_eas_corrected(CCTK_ARGUMENTS);

extern "C" void frankfurt_m1_compute_eas(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if ( verbosity > 1 )
    CCTK_VInfo(CCTK_THORNSTRING,"***** Computing microphysics *****");

  if( CCTK_Equals( m1_which_eas, "microphysical" ) ) {
    frankfurt_m1_compute_eas_microphysical(CCTK_PASS_CTOC);
  } else if ( CCTK_Equals( m1_which_eas, "microphysical-corrected" ) ) {
    frankfurt_m1_compute_eas_corrected(CCTK_PASS_CTOC) ;
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

		if( rho_b[index] < M1_rho_floor * ( 1.0 + M1_atmo_tol ) ) {

			Qnue[index] = 0. ;
			Qnua[index] = 0. ;
			Qnux[index] = 0. ;

			kappa_nue_a[index] = 0.;
			kappa_nua_a[index] = 0.;
			kappa_nux_a[index] = 0.;

			kappa_nue_s[index] = 0.;
			kappa_nua_s[index] = 0.;
			kappa_nux_s[index] = 0.;

			kappa_nue_n[index] = 0.;
			kappa_nua_n[index] = 0.;
			kappa_nux_n[index] = 0.;

			Rnue[index] = 0. ;
			Rnua[index] = 0. ;
			Rnux[index] = 0. ;

			continue ;
			
		}
		
		std::array<CCTK_REAL,REQ_FLUID_VARS> U;
		U[RHO]  = rho_b[index];
		U[YE]   = ye[index];
		U[TEMP] = temp[index];

		//The EAS class is provided by margherita with all the
		//functions and members called below. The EAS class is
		//basically the same as the Leakage class, but with some
		//modifications so that we get the correct opacities/emissivities
		//for M1.
		const auto tau_init = EAS<CCTK_REAL>::compute_analytic_opacity(U[RHO]);
		std::array<double,3> tau_n {{tau_init, tau_init, tau_init}};
		if ( ! CCTK_Equals(m1_which_tau, "simple") ) {
		tau_n[NUE] = tau_nue[index];
		tau_n[NUA] = tau_nua[index];
		tau_n[NUX] = tau_nux[index];
		}
		auto F= Fugacities<CCTK_REAL>(U[RHO], U[TEMP], U[YE], std::move(tau_n) );

		std::array<bool,M1_interactions::NINTERACTIONS> which_interactions { { beta_decay,
											plasmon_decay,
											bremstrahlung,
											pair_annihil }} ;
		
		EAS<CCTK_REAL> eas{F, which_interactions};
		eas.calc_eas(true);

		temp_nue[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue[index], F.eta[NUE]);
		temp_nua[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nua[index], F.eta[NUA]);
		temp_nux[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nux[index], F.eta[NUX]);

		CCTK_REAL const nue_fact = ::max( 1., ::int_pow<2>(temp_nue[index]/temp[index]));
		CCTK_REAL const nua_fact = ::max( 1., ::int_pow<2>(temp_nua[index]/temp[index]));
		CCTK_REAL const nux_fact = ::max( 1., ::int_pow<2>(temp_nux[index]/temp[index]));
		
		Qnue[index] = eas.Q_cactus<NUE>(); 
		Qnua[index] = eas.Q_cactus<NUA>(); 
		Qnux[index] = eas.Q_cactus<NUX>();
		
		kappa_nue_a[index] = eas.kappa_a_cactus<NUE>() * nue_fact ; 
		kappa_nua_a[index] = eas.kappa_a_cactus<NUA>() * nua_fact ; 
		kappa_nux_a[index] = eas.kappa_a_cactus<NUX>() ;
		
		kappa_nue_s[index] = eas.kappa_s_cactus<NUE>() * nue_fact ; 
		kappa_nua_s[index] = eas.kappa_s_cactus<NUA>() * nua_fact ; 
		kappa_nux_s[index] = eas.kappa_s_cactus<NUX>() ;

		Rnue[index] = eas.R_cactus<NUE>();
		Rnua[index] = eas.R_cactus<NUA>();
		Rnux[index] = eas.R_cactus<NUX>();
		
		kappa_nue_n[index] = eas.kappa_n_cactus<NUE>() * nue_fact ; 
		kappa_nua_n[index] = eas.kappa_n_cactus<NUA>() * nua_fact ; 
		kappa_nux_n[index] = eas.kappa_n_cactus<NUX>() ;
		
	} //for loop

	}

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

	//Global quantities for computing avergae over entire grid
	
	CCTK_INT nfix=0, nsucc=0;
	CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2] ;
	CCTK_REAL avgiter=0., error_int=0., error_int_denom=0.  ;

	
	#pragma omp parallel for collapse(3) reduction(+:nfix,nsucc,avgiter,error_int,error_int_denom) schedule(static)
	for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {

		int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

		if( rho_b[index] < M1_rho_floor * ( 1.0 + M1_atmo_tol ) ) {

			Qnue[index] = 0. ;
			Qnua[index] = 0. ;
			Qnux[index] = 0. ;

			kappa_nue_a[index] = 0.;
			kappa_nua_a[index] = 0.;
			kappa_nux_a[index] = 0.;

			kappa_nue_s[index] = 0.;
			kappa_nua_s[index] = 0.;
			kappa_nux_s[index] = 0.;

			kappa_nue_n[index] = 0.;
			kappa_nua_n[index] = 0.;
			kappa_nux_n[index] = 0.;

			Rnue[index] = 0. ;
			Rnua[index] = 0. ;
			Rnux[index] = 0. ;
			
			continue ;
		
		}
		
		CCTK_REAL _T,_ye ;
		
		std::array<CCTK_REAL,REQ_FLUID_VARS> U;
		U[RHO]  = rho_b[index];
		U[YE]   = ye[index];
		U[TEMP] = temp[index];

		_T = U[TEMP] ; _ye = U[YE] ;
		
		//The EAS class is provided by margherita with all the
		//functions and members called below. The EAS class is
		//basically the same as the Leakage class, but with some
		//modifications so that we get the correct opacities/emissivities
		//for M1.
		const auto tau_init = EAS<double>::compute_analytic_opacity(U[RHO]);
		std::array<double,3> tau_n {{tau_init, tau_init, tau_init}};
		if ( ! CCTK_Equals(m1_which_tau, "simple") ) {
		tau_n[NUE] = tau_nue[index];
		tau_n[NUA] = tau_nua[index];
		tau_n[NUX] = tau_nux[index];
		}

		auto F= Fugacities<CCTK_REAL>(U[RHO], U[TEMP], U[YE], tau_n);
		std::array<bool,M1_interactions::NINTERACTIONS> which_interactions { { beta_decay,
											plasmon_decay,
											bremstrahlung,
											pair_annihil }} ;
		
		EAS<CCTK_REAL> eas{F, which_interactions};
		eas.calc_eas();


		CCTK_REAL tau_betaeq_e  = 1./sqrt(eas.kappa_a_cactus<NUE>()*(eas.kappa_a_cactus<NUE>() + eas.kappa_s_cactus<NUE>() ) + 1.e-45) ;
		CCTK_REAL tau_betaeq_a  = 1./sqrt(eas.kappa_a_cactus<NUA>()*(eas.kappa_a_cactus<NUA>() + eas.kappa_s_cactus<NUA>() ) + 1.e-45) ;
		CCTK_REAL tau_betaeq_x  = 1./sqrt(eas.kappa_a_cactus<NUX>()*(eas.kappa_a_cactus<NUX>() + eas.kappa_s_cactus<NUX>() ) + 1.e-45) ;

		// get the minimum beta-equil timescale across the three species 
		CCTK_REAL tau_betaeq = std::min(tau_betaeq_e, std::min(tau_betaeq_a,tau_betaeq_x)) ; 
		CCTK_REAL dt = CCTK_DELTA_TIME ;
		CCTK_REAL beta_equil_tscale_ = std::isnormal(tau_betaeq) ? (tau_betaeq / dt) : 0.0 ;

		
		if(beta_equil_tscale_ < 1.) {
		nfix++ ;
		
		CCTK_REAL ye_eq, T_eq ;
		bool success ;
		unsigned int niter;
		
		std::array<CCTK_REAL,6> Urad;
		Urad[0] = Enue[index]; Urad[1] = Enua[index] ; Urad[2] = Enux[index] ;
		Urad[3] = Nnue[index] ; Urad[4] = Nnua[index] ; Urad[5] = Nnux[index] ;
		std::tie(ye_eq,T_eq,success,niter) = compute_T_ye_betaeq<EOS_Tabulated>(Urad,U, tau_n, which_interactions) ;
		
		if ( success ) {
			nsucc++;
			avgiter += niter ;
			if( beta_equil_tscale_ < 0.5 ) {
				_T = T_eq;
				_ye = ye_eq;
				error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) ;
				error_int_denom += _ye + _T ;
			} else {
				CCTK_REAL fac = 2 * ( 1.0 - beta_equil_tscale_ ) ;
				_T = fac * T_eq  + ( 1 - fac ) * U[TEMP] ;
				_ye = fac * ye_eq + ( 1 - fac ) * U[YE] ;
				error_int += fabs(_T-U[TEMP]) + fabs(_ye-U[YE]) ;
				error_int_denom += _ye + _T ;
			}
			F = Fugacities<CCTK_REAL>(U[RHO],_T, _ye, tau_n) ;
			eas.set_F(F) ;
			eas.calc_eas() ;
			
		} // if success
		}
			
		temp_nue[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nue[index], F.eta[NUE]);
		temp_nua[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nua[index], F.eta[NUA]);
		temp_nux[index] = EAS<CCTK_REAL>::neutrino_temperature(eps_nux[index], F.eta[NUX]);

		CCTK_REAL const nue_fact = ::max( 1., ::int_pow<2>(temp_nue[index]/temp[index]));
		CCTK_REAL const nua_fact = ::max( 1., ::int_pow<2>(temp_nua[index]/temp[index]));

		Qnue[index] = eas.Q_cactus<NUE>() * nue_fact ;
		Qnua[index] = eas.Q_cactus<NUA>() * nua_fact ; 
		Qnux[index] = eas.Q_cactus<NUX>() ;

		Rnue[index] = eas.R_cactus<NUE>() * nue_fact ; 
		Rnua[index] = eas.R_cactus<NUA>() * nua_fact ;
		Rnux[index] = eas.R_cactus<NUX>() ;

		kappa_nue_a[index] = eas.kappa_a_cactus<NUE>() * nue_fact ;
		kappa_nua_a[index] = eas.kappa_a_cactus<NUA>() * nua_fact ;
		kappa_nux_a[index] = eas.kappa_a_cactus<NUX>() ;

		kappa_nue_s[index] = eas.kappa_s_cactus<NUE>() * nue_fact ;
		kappa_nua_s[index] = eas.kappa_s_cactus<NUA>() * nua_fact ;
		kappa_nux_s[index] = eas.kappa_s_cactus<NUX>() ;

		kappa_nue_n[index] = eas.kappa_n_cactus<NUE>() * nue_fact ;
		kappa_nua_n[index] = eas.kappa_n_cactus<NUA>() * nua_fact ;
		kappa_nux_n[index] = eas.kappa_n_cactus<NUX>() ;

		
		eta_nue[index] = F.eta[NUE] ;
		eta_nua[index] = F.eta[NUA] ;
		eta_nux[index] = F.eta[NUX] ;
		
	} //for loop


	if(verbosity>1)
		CCTK_VInfo(CCTK_THORNSTRING,
			"***** EAS: Lev: %d |  Beta equilibrium fixes: applied:success:total %d:%d:%d \t | avg iter/point %1.4f | relative error %1.4g denominator %1.4g *****\n",
			(int)GetRefinementLevel(cctkGH),
			nfix,nsucc,npoints,avgiter/npoints,
			error_int/(error_int_denom+1.e-40),error_int_denom);
  
} // FIL_M1_compute_eas



