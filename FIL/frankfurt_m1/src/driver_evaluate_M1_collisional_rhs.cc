// Harry:  nux in 3 species already include 4 heavy lepton contributions
//  nux in 5 species include 2 heavy lepton contributions
/**
 * @file driver_evaluate_M1_collisional_rhs.cc
 *
 * This file is part of frankfurt_m1
 *
 * Evaluate collisional source terms for M1 
 * equations implicitly via NR rootfinding.
 * The rootfinding is handled by the closure_t
 * class (see m1_collisional.cc)
 *
 * This routine is called after MoL_Add 
 * but BEFORE AddToTmunu ( when C2P happens).
 *
 * Since we treat the sources implicitly 
 * we also need to add them manually to 
 * both radiation and fluid GFs, this
 * is also done in this routine.
 *
 * This routine is based on the initial
 * version of FIL-M1 by Lukas Weih
 * and Elias Most. 
 *
 * @author Carlo Musolino
 */

#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <tuple>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "frankfurt_m1.h" 
#include "m1_implicit_conventions.hh"
#include "Margherita_EOS.h"
#include "m1_closure.hh" 


extern "C" void frankfurt_m1_evaluate_collisional_rhs(CCTK_ARGUMENTS) {
  using namespace frankfurt_m1 ;
  using namespace m1_implicit_conventions ;
  using error_t = m1_implicit_error_bits ;
    
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if ( ! evolve_radiation )
	return ;
  // For now MoL handles the pseudo evolution
  // not sure whether this is actually good
  // from a design perspective
  int const substep = 2 - (*RK_counter);
  
  bool use_initguess = true ;
  
  double dt = CCTK_DELTA_TIME ;
  
  // output stats
  CCTK_INT imin = cctk_nghostzones[0]; CCTK_INT jmin = cctk_nghostzones[1]; CCTK_INT kmin = cctk_nghostzones[2];
  CCTK_INT imax = cctk_lsh[0] - cctk_nghostzones[0] ; CCTK_INT jmax = cctk_lsh[1] - cctk_nghostzones[1] ; CCTK_INT kmax = cctk_lsh[2] -	cctk_nghostzones[2] ;
  CCTK_INT npoints = (imax-imin) * (jmax-jmin) * (kmax-kmin) ;
  CCTK_REAL avgiter=0.;
  CCTK_INT nfix = 0, nfixe_low = 0, nfixe_high = 0, nfail_inv = 0 , nfix_flux = 0, nyefix=0, nymufix=0;
  
#pragma omp parallel for collapse(3) reduction(+:nfix,nfixe_low,nfixe_high,nfail_inv,nfix_flux,avgiter,nyefix,nymufix) 
  for(int k=cctk_nghostzones[2] ;k<cctk_lsh[2]-(cctk_nghostzones[2]);k++) for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-(cctk_nghostzones[1]);j++) for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-(cctk_nghostzones[0]);i++) {

  
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	
	// ATMO check: if we are in Hydro atmosphere no need to make any calculation
	if (rho_b[index] < M1_rho_floor)
	{
	    continue;
	}
	
	taufix[index] = 0.; yefix[index] = 0.; ymufix[index] = 0.;
	m1_e_lepton_source[index] = 0. ;
	m1_mu_lepton_source[index] = 0. ;
	m1_heatcool[index]  = 0. ;
	m1_entropy_source[index] = 0. ;
	
	metric_c gamma( {gxx[index], gxy[index], gxz[index], gyy[index], gyz[index], gzz[index]},
			{betax[index],betay[index],betaz[index]}, alp[index] ) ;

	auto const sqrtgamma = gamma.sqrtdet ;
	
	closure_t closure( {zvecx[index],zvecy[index],zvecz[index]}, &gamma) ;

	
	error_t error_bits_nue, error_bits_nue_bar, error_bits_numu, error_bits_numu_bar, error_bits_nux ;

        if (use_5_spec_m1) {
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		// Step 1: compute collisional sources
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		CCTK_REAL dtau_nue, dtau_nue_bar, dtau_numu, dtau_numu_bar, dtau_nux;
		CCTK_REAL dSx_nue, dSx_nue_bar, dSx_numu, dSx_numu_bar, dSx_nux;
		CCTK_REAL dSy_nue, dSy_nue_bar, dSy_numu, dSy_numu_bar, dSy_nux;
		CCTK_REAL dSz_nue, dSz_nue_bar, dSz_numu, dSz_numu_bar, dSz_nux;
		CCTK_REAL dYE_nue, dYE_nue_bar, dYE_numu, dYE_numu_bar, dYE_nux;
		CCTK_REAL dYMU_nue, dYMU_nue_bar, dYMU_numu, dYMU_numu_bar, dYMU_nux;

		// -------------------------------------------------------------------------------------------------------------------
		//  -- NUE
		// -------------------------------------------------------------------------------------------------------------------
		CCTK_REAL Nnu_e = Nnue_star[index];
		// The solvers works on primitives
		CCTK_REAL Enu_e = Enue_star[index] / gamma.sqrtdet;
		CCTK_REAL Fnux_e = Ft_nue_x[index] / gamma.sqrtdet;
		CCTK_REAL Fnuy_e = Ft_nue_y[index] / gamma.sqrtdet;
		CCTK_REAL Fnuz_e = Ft_nue_z[index] / gamma.sqrtdet;

		// the order here must be consistent with the
		// indices defined in frankfurt_m1.h
		std::array<CCTK_REAL,5> eas_nue { Qnue[index], kappa_nue_a[index], kappa_nue_s[index], Rnue[index], kappa_nue_n[index]};
		
		double zold = zeta_nue[index] ;
		
		// Set internal state before solving
		zeta_nue[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue, Nnu_e, Enu_e, Fnux_e, Fnuy_e, Fnuz_e,
										 dtau_nue, dSx_nue, dSy_nue, dSz_nue, dYE_nue, dYMU_nue, 
										 zold, dt, true,  error_bits_nue ) ;

		if( error_bits_nue[ERRNOCONV] ) {
		  error_bits_nue[ERRNOCONV] = false ;
		  zeta_nue[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue, Nnu_e, Enu_e, Fnux_e, Fnuy_e, Fnuz_e,
										   dtau_nue, dSx_nue, dSy_nue, dSz_nue, dYE_nue, dYMU_nue, 
										   0., dt, false,  error_bits_nue ) ;
		}
		
		// -------------------------------------------------------------------------------------------------------------------
		//  -- NUE_BAR
		// -------------------------------------------------------------------------------------------------------------------
		CCTK_REAL Nnu_e_bar = Nnue_bar_star[index];
		CCTK_REAL Enu_e_bar = Enue_bar_star[index] / gamma.sqrtdet;
		CCTK_REAL Fnux_e_bar = Ft_nue_bar_x[index] / gamma.sqrtdet;
		CCTK_REAL Fnuy_e_bar = Ft_nue_bar_y[index] / gamma.sqrtdet;
		CCTK_REAL Fnuz_e_bar = Ft_nue_bar_z[index] / gamma.sqrtdet;

		zold = zeta_nue_bar[index];
		
		std::array<CCTK_REAL,5> eas_nue_bar { Qnue_bar[index], kappa_nue_bar_a[index], kappa_nue_bar_s[index], Rnue_bar[index], kappa_nue_bar_n[index]};
		
		zeta_nue_bar[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue_bar, Nnu_e_bar, Enu_e_bar, Fnux_e_bar, Fnuy_e_bar, Fnuz_e_bar,
										 dtau_nue_bar, dSx_nue_bar, dSy_nue_bar, dSz_nue_bar, dYE_nue_bar, dYMU_nue_bar, 
										 zold, dt, true,  error_bits_nue_bar );
		if ( error_bits_nue_bar[ERRNOCONV] ) {
		  error_bits_nue_bar[ERRNOCONV] = false ;
		  zeta_nue_bar[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue_bar, Nnu_e_bar, Enu_e_bar, Fnux_e_bar, Fnuy_e_bar, Fnuz_e_bar,
										   dtau_nue_bar, dSx_nue_bar, dSy_nue_bar, dSz_nue_bar, dYE_nue_bar, dYMU_nue_bar, 
										   0., dt, false,  error_bits_nue_bar );
		}

                // -------------------------------------------------------------------------------------------------------------------
                //  -- NUMU
                // -------------------------------------------------------------------------------------------------------------------
                CCTK_REAL Nnu_mu = Nnumu_star[index];
                CCTK_REAL Enu_mu = Enumu_star[index] / gamma.sqrtdet;
                CCTK_REAL Fnux_mu = Ft_numu_x[index] / gamma.sqrtdet;
                CCTK_REAL Fnuy_mu = Ft_numu_y[index] / gamma.sqrtdet;
                CCTK_REAL Fnuz_mu = Ft_numu_z[index] / gamma.sqrtdet;

                zold = zeta_numu[index];

                std::array<CCTK_REAL,5> eas_numu { Qnumu[index], kappa_numu_a[index], kappa_numu_s[index], Rnumu[index], kappa_numu_n[index]};

                zeta_numu[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_numu, Nnu_mu, Enu_mu, Fnux_mu, Fnuy_mu, Fnuz_mu,
                                                                                 dtau_numu, dSx_numu, dSy_numu, dSz_numu, dYE_numu, dYMU_numu, 
                                                                                 zold, dt, true,  error_bits_numu );
                if ( error_bits_numu[ERRNOCONV] ) {
                  error_bits_numu[ERRNOCONV] = false ;
                  zeta_numu[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_numu, Nnu_mu, Enu_mu, Fnux_mu, Fnuy_mu, Fnuz_mu,
                                                                                   dtau_numu, dSx_numu, dSy_numu, dSz_numu, dYE_numu, dYMU_numu, 
                                                                                   0., dt, false,  error_bits_numu );
                }

                // -------------------------------------------------------------------------------------------------------------------
                //  -- NUMU_BAR
                // -------------------------------------------------------------------------------------------------------------------
                CCTK_REAL Nnu_mu_bar = Nnumu_bar_star[index];
                CCTK_REAL Enu_mu_bar = Enumu_bar_star[index] / gamma.sqrtdet;
                CCTK_REAL Fnux_mu_bar = Ft_numu_bar_x[index] / gamma.sqrtdet;
                CCTK_REAL Fnuy_mu_bar = Ft_numu_bar_y[index] / gamma.sqrtdet;
                CCTK_REAL Fnuz_mu_bar = Ft_numu_bar_z[index] / gamma.sqrtdet;

                zold = zeta_numu_bar[index];

                std::array<CCTK_REAL,5> eas_numu_bar { Qnumu_bar[index], kappa_numu_bar_a[index], kappa_numu_bar_s[index], Rnumu_bar[index], kappa_numu_bar_n[index]};

                zeta_numu_bar[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_numu_bar, Nnu_mu_bar, Enu_mu_bar, Fnux_mu_bar, Fnuy_mu_bar, Fnuz_mu_bar,
                                                                                 dtau_numu_bar, dSx_numu_bar, dSy_numu_bar, dSz_numu_bar, dYE_numu_bar, dYMU_numu_bar, 
                                                                                 zold, dt, true,  error_bits_numu_bar );
                if ( error_bits_numu_bar[ERRNOCONV] ) {
                  error_bits_numu_bar[ERRNOCONV] = false ;
                  zeta_numu_bar[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_numu_bar, Nnu_mu_bar, Enu_mu_bar, Fnux_mu_bar, Fnuy_mu_bar, Fnuz_mu_bar,
                                                                                   dtau_numu_bar, dSx_numu_bar, dSy_numu_bar, dSz_numu_bar, dYE_numu_bar, dYMU_numu_bar, 
                                                                                   0., dt, false,  error_bits_numu_bar );
                }


		// -------------------------------------------------------------------------------------------------------------------
		//  -- NUX
		// -------------------------------------------------------------------------------------------------------------------
		CCTK_REAL Nnu_x = Nnux_star[index];
		CCTK_REAL Enu_x = Enux_star[index] / gamma.sqrtdet;
		CCTK_REAL Fnux_x = Ft_nux_x[index] / gamma.sqrtdet;
		CCTK_REAL Fnuy_x = Ft_nux_y[index] / gamma.sqrtdet;
		CCTK_REAL Fnuz_x = Ft_nux_z[index] / gamma.sqrtdet;

		zold = zeta_nux[index] ;
		
		std::array<CCTK_REAL,5> eas_nux { Qnux[index], kappa_nux_a[index], kappa_nux_s[index], Rnux[index], kappa_nux_n[index]};
			
		zeta_nux[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nux, Nnu_x, Enu_x, Fnux_x, Fnuy_x, Fnuz_x,
										 dtau_nux, dSx_nux, dSy_nux, dSz_nux, dYE_nux, dYMU_nux, 
										 zold, dt, true,  error_bits_nux ) ;

		if( error_bits_nux[ERRNOCONV] ) {
		  error_bits_nux[ERRNOCONV] = false ;
		  zeta_nux[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nux, Nnu_x, Enu_x, Fnux_x, Fnuy_x, Fnuz_x,
										   dtau_nux, dSx_nux, dSy_nux, dSz_nux, dYE_nux, dYMU_nux, 
										   0., dt, false,  error_bits_nux ) ;
		}
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		
		// -------------------------------------------------------------------------------------------------------------------
		// Step 2: Finally perform the update. We check that tau>0 and ye, ymu is within table bounds after the update 
		// -------------------------------------------------------------------------------------------------------------------
		bool  number_bad{false}, energy_bad{false};
		// Energy
		double const tau_rhs2    = -couple_fluid * couple_energy * (dtau_nue + dtau_nue_bar + dtau_numu + dtau_numu_bar + dtau_nux);
		// Momenta
		double const st_x_rhs2  = -couple_momenta * couple_fluid * (dSx_nue + dSx_nue_bar + dSx_numu + dSx_numu_bar + dSx_nux);
		double const st_y_rhs2  = -couple_momenta * couple_fluid * (dSy_nue + dSy_nue_bar + dSy_numu + dSy_numu_bar + dSy_nux);
		double const st_z_rhs2  = -couple_momenta * couple_fluid * (dSz_nue + dSz_nue_bar + dSz_numu + dSz_numu_bar + dSz_nux);
        //Harry: does this ye_star_rhs2 have correct with 1/m_b ?
		// Electron fraction
		double ye_star_rhs2 = - couple_fluid * couple_ye   * ( dYE_nue - dYE_nue_bar ) ; 
		// Muon fraction
        //Harry: check here, dYMU_numu - dYMU_numu_bar  good?
		double ymu_star_rhs2 = - couple_fluid * couple_ymu   * ( dYMU_numu - dYMU_numu_bar ) ; 
		// Entropy - Check is this correct? -- CM: I believe so I checked with the paper 
		double s_star_rhs2{0};
		// Check that the update does not make tau < 0
		double const tmptau = tau[index] + dt * tau_rhs2 ;
		if ( tmptau < 0 ) {
		  energy_bad = true ;
		} else {
		  tau[index] += dt * tau_rhs2 ;
		}
		
		mhd_st_x[index]  += dt * st_x_rhs2;
		mhd_st_y[index]  += dt * st_y_rhs2;
		mhd_st_z[index]  += dt * st_z_rhs2;
		// Check that the update doesn't send ye out of table bounds 
		CCTK_REAL const tmpye = ( ye_star[index] + dt * ye_star_rhs2 ) / rho_star[index]; // this is just ye after the implicit step
		// Check that the update doesn't send ymu out of table bounds 
		CCTK_REAL const tmpymu = ( ymu_star[index] + dt * ymu_star_rhs2 ) / rho_star[index]; // this is just ymu after the implicit step
		
		if( (tmpye <  EOS_Leptonic::eos_yemin) ||  (tmpye > EOS_Leptonic::eos_yemax )  ){
		  nyefix ++ ; 
		  number_bad = true ;
		} else {
		  ye_star[index] += dt * ye_star_rhs2;
		}	
// Harry, source term too large --> a lot of nymufix, what should I do?
                if( (tmpymu <  EOS_Leptonic::eos_ymumin) ||  (tmpymu > EOS_Leptonic::eos_ymumax )  ){
                  nymufix ++ ;
                  number_bad = true ;
                } else {
                  ymu_star[index] += dt * ymu_star_rhs2;
                }

	    s_star[index] += dt * s_star_rhs2 ;

	    // just for output 
		m1_e_lepton_source[index] = ye_star_rhs2 * dt / rho_star[index] ;
		m1_mu_lepton_source[index] = ymu_star_rhs2 * dt / rho_star[index] ;
		m1_heatcool[index] = tau_rhs2 * dt ;
		m1_entropy_source[index] = s_star_rhs2 * dt / rho_star[index];
		  
     		
		// Only update Enu if tau was updated 
		if ( ! energy_bad ) {
		  Enue_star[index]       = Enu_e * sqrtgamma;
		  Enue_bar_star[index]   = Enu_e_bar * sqrtgamma;
		  Enumu_star[index]      = Enu_mu * sqrtgamma;
		  Enumu_bar_star[index]  = Enu_mu_bar * sqrtgamma;
		  Enux_star[index]       = Enu_x * sqrtgamma;
		}
		
		Ft_nue_x[index] = Fnux_e* sqrtgamma; Ft_nue_bar_x[index] = Fnux_e_bar* sqrtgamma; Ft_nux_x[index] = Fnux_x* sqrtgamma; 
		Ft_nue_y[index] = Fnuy_e* sqrtgamma; Ft_nue_bar_y[index] = Fnuy_e_bar* sqrtgamma; Ft_nux_y[index] = Fnuy_x* sqrtgamma; 
		Ft_nue_z[index] = Fnuz_e* sqrtgamma; Ft_nue_bar_z[index] = Fnuz_e_bar* sqrtgamma; Ft_nux_z[index] = Fnuz_x* sqrtgamma; 
		Ft_numu_x[index] = Fnux_mu* sqrtgamma; Ft_numu_bar_x[index] = Fnux_mu_bar* sqrtgamma;
		Ft_numu_y[index] = Fnuy_mu* sqrtgamma; Ft_numu_bar_y[index] = Fnuy_mu_bar* sqrtgamma;
		Ft_numu_z[index] = Fnuz_mu* sqrtgamma; Ft_numu_bar_z[index] = Fnuz_mu_bar* sqrtgamma;

		// N is conservative in the solver
		if( !number_bad ){ 
		  Nnue_star[index]      = Nnu_e ;
		  Nnue_bar_star[index]  = Nnu_e_bar ;
		  Nnumu_star[index]     = Nnu_mu ;
		  Nnumu_bar_star[index] = Nnu_mu_bar ;
		  Nnux_star[index]      = Nnu_x ;
		}
		
		// Finally we collect some information about the failures.. 
		nfail_inv += error_bits_nue[ERRNOCONV] + error_bits_nue_bar[ERRNOCONV] + error_bits_numu[ERRNOCONV] + error_bits_numu_bar[ERRNOCONV]+ error_bits_nux[ERRNOCONV] ;

	} else {
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		// Step 1: compute collisional sources
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		CCTK_REAL dtau_nue, dtau_nue_bar, dtau_nux;
		CCTK_REAL dSx_nue, dSx_nue_bar, dSx_nux;
		CCTK_REAL dSy_nue, dSy_nue_bar, dSy_nux;
		CCTK_REAL dSz_nue, dSz_nue_bar, dSz_nux;
		CCTK_REAL dYE_nue, dYE_nue_bar, dYE_nux;
		CCTK_REAL dummy;

		// -------------------------------------------------------------------------------------------------------------------
		//  -- NUE
		// -------------------------------------------------------------------------------------------------------------------
		CCTK_REAL Nnu_e = Nnue_star[index];
		// The solvers works on primitives
		CCTK_REAL Enu_e = Enue_star[index] / gamma.sqrtdet;
		CCTK_REAL Fnux_e = Ft_nue_x[index] / gamma.sqrtdet;
		CCTK_REAL Fnuy_e = Ft_nue_y[index] / gamma.sqrtdet;
		CCTK_REAL Fnuz_e = Ft_nue_z[index] / gamma.sqrtdet;

		// the order here must be consistent with the
		// indices defined in frankfurt_m1.h
		std::array<CCTK_REAL,5> eas_nue { Qnue[index], kappa_nue_a[index], kappa_nue_s[index], Rnue[index], kappa_nue_n[index]};
		
		double zold = zeta_nue[index] ;
		
		// Set internal state before solving
		zeta_nue[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue, Nnu_e, Enu_e, Fnux_e, Fnuy_e, Fnuz_e,
										 dtau_nue, dSx_nue, dSy_nue, dSz_nue, dYE_nue, dummy, 
										 zold, dt, true,  error_bits_nue ) ;

		if( error_bits_nue[ERRNOCONV] ) {
		  error_bits_nue[ERRNOCONV] = false ;
		  zeta_nue[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue, Nnu_e, Enu_e, Fnux_e, Fnuy_e, Fnuz_e,
										   dtau_nue, dSx_nue, dSy_nue, dSz_nue, dYE_nue, dummy, 
										   0., dt, false,  error_bits_nue ) ;
		}
		
		// -------------------------------------------------------------------------------------------------------------------
		//  -- NUE_BAR
		// -------------------------------------------------------------------------------------------------------------------
		CCTK_REAL Nnu_e_bar = Nnue_bar_star[index];
		CCTK_REAL Enu_e_bar = Enue_bar_star[index] / gamma.sqrtdet;
		CCTK_REAL Fnux_e_bar = Ft_nue_bar_x[index] / gamma.sqrtdet;
		CCTK_REAL Fnuy_e_bar = Ft_nue_bar_y[index] / gamma.sqrtdet;
		CCTK_REAL Fnuz_e_bar = Ft_nue_bar_z[index] / gamma.sqrtdet;

		zold = zeta_nue_bar[index];
		
		std::array<CCTK_REAL,5> eas_nue_bar { Qnue_bar[index], kappa_nue_bar_a[index], kappa_nue_bar_s[index], Rnue_bar[index], kappa_nue_bar_n[index]};
		
		zeta_nue_bar[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue_bar, Nnu_e_bar, Enu_e_bar, Fnux_e_bar, Fnuy_e_bar, Fnuz_e_bar,
										 dtau_nue_bar, dSx_nue_bar, dSy_nue_bar, dSz_nue_bar, dYE_nue_bar, dummy, 
										 zold, dt, true,  error_bits_nue_bar );
		if ( error_bits_nue_bar[ERRNOCONV] ) {
		  error_bits_nue_bar[ERRNOCONV] = false ;
		  zeta_nue_bar[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue_bar, Nnu_e_bar, Enu_e_bar, Fnux_e_bar, Fnuy_e_bar, Fnuz_e_bar,
										   dtau_nue_bar, dSx_nue_bar, dSy_nue_bar, dSz_nue_bar, dYE_nue_bar, dummy, 
										   0., dt, false,  error_bits_nue_bar );
		}
		// -------------------------------------------------------------------------------------------------------------------
		//  -- NUX
		// -------------------------------------------------------------------------------------------------------------------
		CCTK_REAL Nnu_x = Nnux_star[index];
		CCTK_REAL Enu_x = Enux_star[index] / gamma.sqrtdet;
		CCTK_REAL Fnux_x = Ft_nux_x[index] / gamma.sqrtdet;
		CCTK_REAL Fnuy_x = Ft_nux_y[index] / gamma.sqrtdet;
		CCTK_REAL Fnuz_x = Ft_nux_z[index] / gamma.sqrtdet;

		zold = zeta_nux[index] ;
		
		std::array<CCTK_REAL,5> eas_nux { Qnux[index], kappa_nux_a[index], kappa_nux_s[index], Rnux[index], kappa_nux_n[index]};
			
		zeta_nux[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nux, Nnu_x, Enu_x, Fnux_x, Fnuy_x, Fnuz_x,
										 dtau_nux, dSx_nux, dSy_nux, dSz_nux, dYE_nux, dummy, 
										 zold, dt, true,  error_bits_nux ) ;

		if( error_bits_nux[ERRNOCONV] ) {
		  error_bits_nux[ERRNOCONV] = false ;
		  zeta_nux[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nux, Nnu_x, Enu_x, Fnux_x, Fnuy_x, Fnuz_x,
										   dtau_nux, dSx_nux, dSy_nux, dSz_nux, dYE_nux, dummy, 
										   0., dt, false,  error_bits_nux ) ;
		}
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		// -------------------------------------------------------------------------------------------------------------------
		
		// -------------------------------------------------------------------------------------------------------------------
		// Step 2: Finally perform the update. We check that tau>0 and ye is within table bounds after the update 
		// -------------------------------------------------------------------------------------------------------------------
		bool  number_bad{false}, energy_bad{false};
		// Energy
		double const tau_rhs2    = -couple_fluid * couple_energy * (dtau_nue + dtau_nue_bar + dtau_nux);
		// Momenta
		double const st_x_rhs2  = -couple_momenta * couple_fluid * (dSx_nue + dSx_nue_bar + dSx_nux);
		double const st_y_rhs2  = -couple_momenta * couple_fluid * (dSy_nue + dSy_nue_bar + dSy_nux);
		double const st_z_rhs2  = -couple_momenta * couple_fluid * (dSz_nue + dSz_nue_bar + dSz_nux);
		// Electron fraction
		double ye_star_rhs2 = - couple_fluid * couple_ye   * ( dYE_nue - dYE_nue_bar ) ; 
		// Muon fraction
		double ymu_star_rhs2 = 0.0 ;
		// Entropy - Check is this correct? -- CM: I believe so I checked with the paper 
		double s_star_rhs2{0};
		// Check that the update does not make tau < 0
		double const tmptau = tau[index] + dt * tau_rhs2 ;
		if ( tmptau < 0 ) {
		  energy_bad = true ;
		} else {
		  tau[index] += dt * tau_rhs2 ;
		}
		
		mhd_st_x[index]  += dt * st_x_rhs2;
		mhd_st_y[index]  += dt * st_y_rhs2;
		mhd_st_z[index]  += dt * st_z_rhs2;
		// Check that the update doesn't send ye out of table bounds 
		CCTK_REAL const tmpye = ( ye_star[index] + dt * ye_star_rhs2 ) / rho_star[index]; // this is just ye after the implicit step
		
		if( (tmpye <  EOS_Tabulated::eos_yemin) ||  (tmpye > EOS_Tabulated::eos_yemax )  ){
		  nyefix ++ ; 
		  number_bad = true ;
		} else {
		  ye_star[index] += dt * ye_star_rhs2;
		}	
		ymu_star[index] += dt * 0.0 ;
		s_star[index] += dt * s_star_rhs2 ;
		// just for output 
		m1_e_lepton_source[index] = ye_star_rhs2 * dt / rho_star[index] ;
		m1_mu_lepton_source[index] = 0.0 ;
		m1_heatcool[index] = tau_rhs2 * dt ;
		m1_entropy_source[index] = s_star_rhs2 * dt / rho_star[index];
	
		// Only update Enu if tau was updated 
		if ( ! energy_bad ) {
		  Enue_star[index]       = Enu_e * sqrtgamma;
		  Enue_bar_star[index]   = Enu_e_bar * sqrtgamma;
		  Enux_star[index]       = Enu_x * sqrtgamma;
		}
		
		Ft_nue_x[index] = Fnux_e* sqrtgamma; Ft_nue_bar_x[index] = Fnux_e_bar* sqrtgamma; Ft_nux_x[index] = Fnux_x* sqrtgamma; 
		Ft_nue_y[index] = Fnuy_e* sqrtgamma; Ft_nue_bar_y[index] = Fnuy_e_bar* sqrtgamma; Ft_nux_y[index] = Fnuy_x* sqrtgamma; 
		Ft_nue_z[index] = Fnuz_e* sqrtgamma; Ft_nue_bar_z[index] = Fnuz_e_bar* sqrtgamma; Ft_nux_z[index] = Fnuz_x* sqrtgamma; 

		// N is conservative in the solver
		if( !number_bad ){ 
		  Nnue_star[index] = Nnu_e ;
		  Nnue_bar_star[index] = Nnu_e_bar ;
		  Nnux_star[index] = Nnu_x ;
		}
		
		// Finally we collect some information about the failures.. 
		nfail_inv += error_bits_nue[ERRNOCONV] + error_bits_nue_bar[ERRNOCONV] + error_bits_nux[ERRNOCONV] ;

        } // use_5_spec_m1
    } // for loop
	
  avgiter /= npoints ;
  nfix = nyefix + nfail_inv ;
  
  if (verbosity)
    CCTK_VInfo(CCTK_THORNSTRING,
	       "****** Implicit dt at Lev. %d | avg, iter/point: %1.5f | nfix : npoints \t %d : %d | ye_outofbound Rootfinder failure \t %d : %d ",
	       (int)GetRefinementLevel(cctkGH),
	       avgiter,nfix, npoints,
	       nyefix, nfail_inv ) ;
  return;
}

