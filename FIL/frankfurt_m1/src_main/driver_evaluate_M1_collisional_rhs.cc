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
  CCTK_INT nfix = 0, nfixe_low = 0, nfixe_high = 0, nfail_inv = 0 , nfix_flux = 0, nyefix=0 ;
  
#pragma omp parallel for collapse(3) reduction(+:nfix,nfixe_low,nfixe_high,nfail_inv,nfix_flux,avgiter,nyefix) 
  for(int k=cctk_nghostzones[2] ;k<cctk_lsh[2]-(cctk_nghostzones[2]);k++) for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-(cctk_nghostzones[1]);j++) for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-(cctk_nghostzones[0]);i++) {

  
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	
	// ATMO check: if we are in Hydro atmosphere no need to make any calculation
	if (rho_b[index] < M1_rho_floor)
	{
	    continue;
	}
	
	taufix[index] = 0.; yefix[index] = 0.;
	m1_lepton_source[index] = 0. ;
	m1_heatcool[index]  = 0. ;
	m1_entropy_source[index] = 0. ;
	
	metric_c gamma( {gxx[index], gxy[index], gxz[index], gyy[index], gyz[index], gzz[index]},
			{betax[index],betay[index],betaz[index]}, alp[index] ) ;

	auto const sqrtgamma = gamma.sqrtdet ;
	
	closure_t closure( {zvecx[index],zvecy[index],zvecz[index]}, &gamma) ;

	
	error_t error_bits_e, error_bits_a, error_bits_x ;
	
	// -------------------------------------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------------
	// Step 1: compute collisional sources
	// -------------------------------------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------------
	CCTK_REAL dtau_nue, dtau_nua, dtau_nux;
	CCTK_REAL dSx_nue, dSx_nua, dSx_nux;
	CCTK_REAL dSy_nue, dSy_nua, dSy_nux;
	CCTK_REAL dSz_nue, dSz_nua, dSz_nux;
	CCTK_REAL dYE_nue, dYE_nua, dYE_nux;

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
	
	double zold = zeta_e[index] ;
	
	// Set internal state before solving
	zeta_e[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue, Nnu_e, Enu_e, Fnux_e, Fnuy_e, Fnuz_e,
									 dtau_nue, dSx_nue, dSy_nue, dSz_nue, dYE_nue,
									 zold, dt, true,  error_bits_e ) ;

	if( error_bits_e[ERRNOCONV] ) {
	  error_bits_e[ERRNOCONV] = false ;
	  zeta_e[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nue, Nnu_e, Enu_e, Fnux_e, Fnuy_e, Fnuz_e,
									   dtau_nue, dSx_nue, dSy_nue, dSz_nue, dYE_nue,
									   0., dt, false,  error_bits_e ) ;
	}
	
	// -------------------------------------------------------------------------------------------------------------------
	//  -- NUA
	// -------------------------------------------------------------------------------------------------------------------
	CCTK_REAL Nnu_a = Nnua_star[index];
	CCTK_REAL Enu_a = Enua_star[index] / gamma.sqrtdet;
	CCTK_REAL Fnux_a = Ft_nua_x[index] / gamma.sqrtdet;
	CCTK_REAL Fnuy_a = Ft_nua_y[index] / gamma.sqrtdet;
	CCTK_REAL Fnuz_a = Ft_nua_z[index] / gamma.sqrtdet;

	zold = zeta_a[index];
	
	std::array<CCTK_REAL,5> eas_nua { Qnua[index], kappa_nua_a[index], kappa_nua_s[index], Rnua[index], kappa_nua_n[index]};
	
	zeta_a[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nua, Nnu_a, Enu_a, Fnux_a, Fnuy_a, Fnuz_a,
									 dtau_nua, dSx_nua, dSy_nua, dSz_nua, dYE_nua,
									 zold, dt, true,  error_bits_a );
	if ( error_bits_a[ERRNOCONV] ) {
	  error_bits_a[ERRNOCONV] = false ;
	  zeta_a[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nua, Nnu_a, Enu_a, Fnux_a, Fnuy_a, Fnuz_a,
									   dtau_nua, dSx_nua, dSy_nua, dSz_nua, dYE_nua,
									   0., dt, false,  error_bits_a );
	}
	// -------------------------------------------------------------------------------------------------------------------
	//  -- NUX
	// -------------------------------------------------------------------------------------------------------------------
	CCTK_REAL Nnu_x = Nnux_star[index];
	CCTK_REAL Enu_x = Enux_star[index] / gamma.sqrtdet;
	CCTK_REAL Fnux_x = Ft_nux_x[index] / gamma.sqrtdet;
	CCTK_REAL Fnuy_x = Ft_nux_y[index] / gamma.sqrtdet;
	CCTK_REAL Fnuz_x = Ft_nux_z[index] / gamma.sqrtdet;

	zold = zeta_x[index] ;
	
	std::array<CCTK_REAL,5> eas_nux { Qnux[index], kappa_nux_a[index], kappa_nux_s[index], Rnux[index], kappa_nux_n[index]};
		
	zeta_x[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nux, Nnu_x, Enu_x, Fnux_x, Fnuy_x, Fnuz_x,
									 dtau_nux, dSx_nux, dSy_nux, dSz_nux, dYE_nux,
									 zold, dt, true,  error_bits_x ) ;

	if( error_bits_x[ERRNOCONV] ) {
	  error_bits_x[ERRNOCONV] = false ;
	  zeta_x[index] = closure.get_collisional_rootfinder<m1_closure_f>(eas_nux, Nnu_x, Enu_x, Fnux_x, Fnuy_x, Fnuz_x,
									   dtau_nux, dSx_nux, dSy_nux, dSz_nux, dYE_nux,
									   0., dt, false,  error_bits_x ) ;
	}
	// -------------------------------------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------------
	
	// -------------------------------------------------------------------------------------------------------------------
	// Step 2: Finally perform the update. We check that tau>0 and ye is within table bounds after the update 
	// -------------------------------------------------------------------------------------------------------------------
	bool  number_bad{false}, energy_bad{false};
	if ( substep==1 ) { // only update fluid on second substep of the method
	  // Energy
	  double const tau_rhs2    = -couple_fluid * couple_energy * (dtau_nue + dtau_nua + dtau_nux);
	  // Momenta
	  double const st_x_rhs2  = -couple_momenta * couple_fluid * (dSx_nue + dSx_nua + dSx_nux);
	  double const st_y_rhs2  = -couple_momenta * couple_fluid * (dSy_nue + dSy_nua + dSy_nux);
	  double const st_z_rhs2  = -couple_momenta * couple_fluid * (dSz_nue + dSz_nua + dSz_nux);
	  // Electron fraction
	  double ye_star_rhs2 = - couple_fluid * couple_ye   * ( dYE_nue - dYE_nua ) ; 
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
	  s_star[index] += dt * s_star_rhs2 ;

	  // just for output 
	  m1_lepton_source[index] = ye_star_rhs2 * dt / rho_star[index] ;
	  m1_heatcool[index] = tau_rhs2 * dt ;
	  m1_entropy_source[index] = s_star_rhs2 * dt / rho_star[index];
	  
	} // if substep == 1

       	
	// Only update Enu if tau was updated 
	if ( ! energy_bad ) {
	  Enue_star[index]   = Enu_e * sqrtgamma;
	  Enua_star[index]   = Enu_a * sqrtgamma;
	  Enux_star[index]   = Enu_x * sqrtgamma;
	}
	
	Ft_nue_x[index] = Fnux_e* sqrtgamma; Ft_nua_x[index] = Fnux_a* sqrtgamma; Ft_nux_x[index] = Fnux_x* sqrtgamma; 
	Ft_nue_y[index] = Fnuy_e* sqrtgamma; Ft_nua_y[index] = Fnuy_a* sqrtgamma; Ft_nux_y[index] = Fnuy_x* sqrtgamma; 
	Ft_nue_z[index] = Fnuz_e* sqrtgamma; Ft_nua_z[index] = Fnuz_a* sqrtgamma; Ft_nux_z[index] = Fnuz_x* sqrtgamma; 

	// N is conservative in the solver
	if( !number_bad ){ 
	  Nnue_star[index] = Nnu_e ;
	  Nnua_star[index] = Nnu_a ;
	  Nnux_star[index] = Nnu_x ;
	}
	
	// Finally we collect some information about the failures.. 
	nfail_inv += error_bits_e[ERRNOCONV] + error_bits_a[ERRNOCONV] + error_bits_x[ERRNOCONV] ;
    }
	
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

