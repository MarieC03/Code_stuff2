/**
 * @file m1_closure.cc
 * This file is part of frankfurt_m1
 *
 * This file contains:
 * Ctor for closure_t class
 * update_closure routine plus its helpers 
 *
 * @author Carlo Musolino
 */

#include "m1_closure.hh"
#include <iostream>
#include "utils/brent.h"
#include "utils/newton.h"

#include <tuple>
namespace frankfurt_m1 {

// ===================================
// Minerbo closure
// ===================================

  /** 
   * Ctor. Initialises velocity and metric pter
   * also computes v2 and all needed powers of W
   */ 
  closure_t::closure_t(vec_t && __zvec, metric_c * __gamma ) : vel(std::move(__zvec)), gamma(__gamma)  
    {
      // vel contains zvec right now
      W2 = 1. + gamma->square_3_upper<0,3>(vel); 
      W  = sqrt(W2) ;

      if ( W > Wmax ) {
	auto const factor = sqrt(zmax2 / (W2-1.)) ;
	for( auto& vv: vel) vv *= factor ;
	W = Wmax;
	W2 = W*W ;
      }

      for( auto& vv: vel) vv/=W ;
      
      vlow = gamma->lower_index<0,3>(vel) ;
      v2 = gamma->square_3_upper<0,3>(vel) ;
      W3 = W*W2 ;
    };
  
  // ===================================
  // Closure helpers 
  // ===================================
  /**
   * When computing the closure (zeta)
   * we solve for xi(z) = 0
   * so this function computes xi
   */
  template<class closure_f>
  CCTK_REAL closure_t::get_xi(CCTK_REAL const& zeta) {
    const CCTK_REAL chi = closure_f::closure_func(zeta) ;
    const CCTK_REAL d_thin  = 1.5*chi - 0.5 ;
    const CCTK_REAL d_thick = 1. - d_thin ;
    const CCTK_REAL J_closure = B0 + d_thin * Bthin + d_thick*Bthick ;
    const CCTK_REAL Hsq = HH0 + d_thick*HHT + d_thin * HHt +
      d_thick * d_thin * HHTt + d_thick * d_thick * HHTT
      + d_thin * d_thin * HHtt ;
    return ( (J_closure*J_closure*zeta*zeta) - Hsq ) / ( E_closure * E_closure ) ;
  }
  /**
   * Same as above with derivative 
   */
  template<class closure_f>
  std::pair<CCTK_REAL,CCTK_REAL> closure_t::get_xi_dxi(CCTK_REAL const & zeta ) {
    CCTK_REAL dchidz ;
    const CCTK_REAL chi = closure_f::closure_func(zeta, dchidz) ;
    const CCTK_REAL d_thin  = 1.5*chi - 0.5 ;
    const CCTK_REAL d_thick = 1. - d_thin ;
    const CCTK_REAL d_thin_dz = 1.5 * dchidz ;
    const CCTK_REAL d_thick_dz = -1.5*dchidz ;
    
    const CCTK_REAL J_closure = B0 + d_thin * Bthin + d_thick*Bthick ;
    const CCTK_REAL dJ_dz = Bthin * d_thin_dz + Bthick * d_thick_dz ;
    
    const CCTK_REAL Hsq = HH0 + d_thick*HHT + d_thin * HHt +
      d_thick * d_thin * HHTt + d_thick * d_thick * HHTT
      + d_thin * d_thin * HHtt ;
    const CCTK_REAL dHsq_dd_thin = HHt + HHTt * d_thick  + 2.* HHtt * d_thin ;
    const CCTK_REAL dHsq_dd_thick = HHT + HHTt * d_thin  + 2.* HHTT * d_thick ;
    
    const CCTK_REAL f  = ( ::int_pow<2>(J_closure*zeta) - Hsq ) / ( E_closure * E_closure ) ;
    CCTK_REAL df = 2.*J_closure*dJ_dz*::int_pow<2>(zeta) + 2.*::int_pow<2>(J_closure)*zeta
      - d_thin_dz*dHsq_dd_thin - d_thick_dz*dHsq_dd_thick;
    df /= (E_closure*E_closure);
    return std::make_pair(f, df) ;
    
  }

  // ========================================================
  // Compute closure and all quantities needed for Jacobian
  // =========================================================
  
  /**
   * First: fill the Eulerian radiation moments 
   * and limit F such that F^i F_i <= E^2 (causality)
   * Second: if requested update the closure via NR or
   * Brent if that fails
   * Third: Compute P^{\mu\nu} H_{\mu} and J (fluid frame moments)
   * based on the updated closure
   * @return Eddington factor
   */
  template<class closure_f>
  CCTK_REAL closure_t::update_closure(CCTK_REAL const& __E,
				      std::array<CCTK_REAL,3> && __F_lo,
				      CCTK_REAL const& zold, bool get_zeta=true) {
  
    
    F_lo = std::move(__F_lo) ;
    F2 = gamma->square_3_lower<0,3>(F_lo) + TINY;
    E = std::max(__E, 1.-15) ;
    // Enforce causality
    
    if( (F2 > E*E) ) {
      auto const factor = sqrt(E*E/F2) * 0.999999 ;
      for( auto &FF: F_lo) FF*=factor ;
      F2 = gamma->square_3_lower<0,3>(F_lo) + TINY;
    }
    F = sqrt(F2) ;
    F_hi = gamma->raise_index<0,3>(F_lo) ;
    for( int ii=0; ii<3; ii++) {
    f_lo[ii] = F_lo[ii]/F ;
    f_hi[ii] = F_hi[ii]/F ;
    }			  
    fv = gamma->scalar_product(vel,f_lo) + TINY;
    Fv = gamma->scalar_product(vel,F_lo) + TINY;
    Ff = gamma->scalar_product(F_lo,f_hi) + TINY;

    // Pieces of J
    B0 = W2 * ( E - 2.*Fv ) ;
    Bthin = W2 * E * ::int_pow<2>(fv) ;
    Bthick = (W2-1.) / (1. + 2. * W2) * (4. * W2 * Fv + (3.-2.*W2)*E) ;

    // coefficients for various parts
    // of H_\mu proportional to v
    // n F and f in optically thick
    // and thin regimes
    an = W*B0 + W*(Fv-E) ;
    av = W*B0 ;
    aF = -W ;

    an_thick = W*Bthick ;
    av_thick = W*Bthick + W/(2.*W2 + 1.)*((3.-2.*W2)*E
					+ (2.*W2-1.)*Fv ) ;

    aF_thick = W*v2 ;
    
    an_thin = W*Bthin ;
    av_thin = an_thin ;
    af_thin = W*E*fv  ;

    if ( get_zeta ) {
      
      if ( v2 < 1.e-15 ) {
	zeta = std::sqrt(F2 / (E*E+TINY)) ;
      } else {
	
	E_closure = E ;
	
	// Pieces of H_\mu H^\mu
	HH0 = ::int_pow<2>(av) * v2 + ::int_pow<2>(aF)*F2 + 2.*av*aF*Fv - ::int_pow<2>(an) ;
	
	
	HHTT = ::int_pow<2>(av_thick)*v2 + ::int_pow<2>(aF_thick)*F2
	  + 2. * aF_thick*av_thick * Fv - an_thick*an_thick ;
	HHtt = ::int_pow<2>(av_thin) * v2 + ::int_pow<2>(af_thin) + 2.* af_thin*av_thin*fv 
	  - ::int_pow<2>(an_thin) ;
	
	HHt = 2.*(av*av_thin*v2 + aF*af_thin*Ff
		  + af_thin*av*fv + av_thin * aF* Fv - an_thin*an ) ;
	HHT = 2.*(av * av_thick * v2 + aF*aF_thick*F2
		  + aF_thick*av*Fv + av_thick*aF*Fv - an_thick*an ) ;
	
	HHTt = 2.*(av_thin*av_thick*v2 + af_thin*aF_thick*Ff + af_thin*av_thick*fv
		   + av_thin*aF_thick*Fv - an_thin*an_thick ) ;
	
	constexpr const CCTK_REAL macheps = std::numeric_limits<CCTK_REAL>::epsilon() ;
	auto const stopif = [] (CCTK_REAL const& a, CCTK_REAL const& b) {
			      return std::fabs(b-a) < m1_implicit_conventions::CLOSURE_TOL + std::fabs(a)*macheps*2. ;
			    };

	size_t iter{ m1_implicit_conventions::MAXITER_NR_CLOSURE } ; 
	
	auto const f = [&] ( CCTK_REAL const& z) {
	  return get_xi<closure_f>(z) ;
	} ;
	
	auto const fdf = [&] ( CCTK_REAL const& z) {
	  return get_xi_dxi<closure_f>(z) ;
	} ;
	
	zeta = m1_rootfinders::rootfind_constrained_newton_raphson(0., 1., fdf, stopif, iter ) ;
	
#ifdef DBG_CLOSURE
	if ( iter>= m1_implicit_conventions::MAXITER_NR_CLOSURE ) { 
	  std::cout << "Newton Raphson failed..." << std::endl ;
	} else {
	  std::cout << "NR converged after " << iter << " iterations\n";
	}
#endif
	bool do_brent = false;
	
	if ( zeta  > 1. || zeta < 0. )
	  do_brent = true;
      
	if ( !do_brent && iter >= m1_implicit_conventions::MAXITER_NR_CLOSURE){
	  do_brent = true ;
	}// good old brent show us the way
	
	if  ( do_brent ) 
	  zeta = zero_brent(0.,1.,m1_implicit_conventions::CLOSURE_TOL , f ) ;
	
      } // if v2 > 1.e-15 
    }else { // if get_zeta 
      zeta = zold ;
    }
    /*
    if ( zeta > 1. - 1.e-10 )
      zeta = 1. ;
    else if ( zeta < 1.e-10 )
      zeta = 0. ;
    */
    chi = closure_f::closure_func(zeta) ;
    dthin  = 1.5*chi - 0.5 ;
    dthick = 1.5 - 1.5*chi ;

    // J and H_hi in the optically thick limit
    // Foucart 2015 A5 A6
    // Radice 2021 12 13 
    const CCTK_REAL JT = 3./(2.*W2 + 1. ) * ((2.*W2 - 1.)*E - 2.*W2*Fv) ;  
    vec_t tHt ;
    for( int ii=0; ii<3; ii++ )
      tHt[ii] = F_hi[ii]/W + vel[ii]*W/(2.*W2+1.)*((4.*W2+1.)*Fv-4.*W2*E);
    //tHt[ii] = F_hi[ii]/W + W*vel[ii]*(Fv - E + JT) - 4./3.*W*vel[ii]*JT ; < -- same thing 
    
    // compute pressure
    // Thick: Radice (8)
    // Thin: Weih Radice Foucart Shibata
    int idx =0 ;
    for( int ii=0; ii<3; ii++)
      for( int jj=ii; jj<3; jj++) {
	// dthick Pthick
	P[idx] = dthick * ( JT/3.* (4.*W2*vel[ii]*vel[jj] + gamma->invmetric_comp(idx)) + W*(tHt[ii]*vel[jj] + tHt[jj]*vel[ii]) ) ;
	// + dthin Pthin
	P[idx]+= dthin  * E*f_hi[ii]*f_hi[jj] ; // Zelmani has an extra W2 here 
	idx++ ;
      }
  
    // And fluid rest frame quantities
    // Radice (48)
    J = B0 + Bthin*dthin + Bthick*dthick;
    Hn = an + an_thick*dthick + an_thin * dthin ;
    // Radice (49)
    for( int ii=0; ii<3; ii++){
      H_lo[ii] = -(av + dthin*av_thin + dthick*av_thick) * vlow[ii]
	         -(aF + dthick*aF_thick) * F_lo[ii]
	         -dthin*af_thin*f_lo[ii] ; // corrected Zelmani
    }
    H_hi = gamma->raise_index<0,3>(H_lo) ;
    
    // normalisation for neutrino n density
    Gamma = W * ( E - Fv ) / (J+TINY) + TINY;
    
    return zeta ;
  }
  
  // Instantiate templates
  template CCTK_REAL closure_t::update_closure<m1_closure_f>(CCTK_REAL const&,
							     std::array<CCTK_REAL,3>&&,
							     CCTK_REAL const&,
							     bool) ;
  template std::pair<CCTK_REAL,CCTK_REAL> closure_t::get_xi_dxi<m1_closure_f>(CCTK_REAL const& ) ;
  template CCTK_REAL closure_t::get_xi<m1_closure_f>(CCTK_REAL const& ) ;
  
}
