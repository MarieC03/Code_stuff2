/**
 * @file m1_fluxes.cc
 *
 * This file is part of frankfurt_m1
 *
 * This file contains routines to compute
 * wavespeeds and fluxes
 * @author Carlo Musolino
 */
#include "m1_closure.hh"
#include <cmath>

namespace frankfurt_m1 {
  /** 
   * Get wavespeeds
   * The equations can be found in Shibata
   * and Foucart 
   * Note that the closure must always be up to date
   * before this routine is called 
   * @return: maximum left and right going eigenspeeds
   */ 
  template<int flux_dirn>
  void closure_t::get_wavespeeds(CCTK_REAL& cp, CCTK_REAL& cm) {
    auto const alp = gamma->lapse ; auto const sqrtg = gamma->sqrtdet ;
    auto const betad = gamma->get_shift(flux_dirn) ;
    int index = 0 * (flux_dirn==0) + 3 * ( flux_dirn==1 ) + 5 * (flux_dirn==2) ;
    auto const gddu  = gamma->invmetric_comp(index);
    
    // All of this is Foucart+ 2015 (C2-8)
    CCTK_REAL l_thin_m = -betad - alp*fabs(F_hi[flux_dirn])/F ;
    CCTK_REAL l_thin_p = -betad + alp*fabs(F_hi[flux_dirn])/F ;
    
    CCTK_REAL p = alp*vel[flux_dirn]/W ;
    CCTK_REAL r = sqrt(alp*alp*gddu*(2.*W2+1.)-2.*W2*p*p ) ;
    
    CCTK_REAL l_thick_m = std::min( -betad+(2.*W2*p-r)/(2.*W2+1.), -betad+p) ;
    CCTK_REAL l_thick_p = std::max( -betad+(2.*W2*p+r)/(2.*W2+1.), -betad+p) ;
    
    cp = dthin*l_thin_p + dthick*l_thick_p ;
    cm = dthin*l_thin_m + dthick*l_thick_m ;
  }
  // ==========================================================================
  /** 
   * Get numerical fluxes -- See Radice Foucart Shibata Weih for the eqs
   * \partial_t \sqrtgamma E = (...) - \partial_i \sqrtgamma  ( \alpha F^i - \beta^i E )
   * \partial_t \sqrtgamma F_k = (...) - \partial_i \sqrtgamma  ( \alpha P^i_k - \beta^i F_k )
   */
  template<int flux_dirn> 
  void closure_t::get_fluxes(CCTK_REAL const& n,
			     CCTK_REAL &Nflux, CCTK_REAL &Eflux,
			     CCTK_REAL &F_x_flux, CCTK_REAL &F_y_flux,
			     CCTK_REAL &F_z_flux ) {
    
  // vec_t == std::array<CCTK_REAL,3>
    vec_t const Pdmu = { (flux_dirn==0)*P[0] + (flux_dirn==1)*P[1] + (flux_dirn==2)*P[2],
			 (flux_dirn==0)*P[1] + (flux_dirn==1)*P[3] + (flux_dirn==2)*P[4],
			 (flux_dirn==0)*P[2] + (flux_dirn==1)*P[4] + (flux_dirn==2)*P[5] } ;
    
    // This is P^(flux_dirn)_i 
    vec_t const Pdi_lo = gamma->lower_index<0,3>(Pdmu) ;
    
    auto const alp = gamma->lapse ; auto const sqrtg = gamma->sqrtdet ;
    auto const betad = gamma->get_shift(flux_dirn) ;
    
    // DR's way, like it more -- Radice+ 2021 (23) and (28)
    Nflux = sqrtg*alp*n*( W*(vel[flux_dirn]-betad/alp) + H_hi[flux_dirn]/(TINY+J) ); 
    //Nflux = sqrtg * (alp*W*vel[flux_dirn]*J +
    //		   alp*H_hi[flux_dirn] - betad* N ); // -- this was Foucart - Zelmani
    
    // sqrtg alpha * ( F^i - \beta^i/alpha * E ) 
    Eflux = sqrtg * (alp*F_hi[flux_dirn] - betad* E );

    // sqrtg alpha * ( P^i_k - \beta^i/alpha * F_k ) 
    F_x_flux = sqrtg * ( alp * Pdi_lo[0] - betad*F_lo[0] ) ;
    F_y_flux = sqrtg * ( alp * Pdi_lo[1] - betad*F_lo[1] ) ;
    F_z_flux = sqrtg * ( alp * Pdi_lo[2] - betad*F_lo[2] ) ;
    return ;
  }
  // ==========================================================================
  // ==========================================================================

  // Instantiate the templates
  template void closure_t::get_fluxes<0>(CCTK_REAL const&,
					 CCTK_REAL &, CCTK_REAL &,
					 CCTK_REAL &, CCTK_REAL &,
					 CCTK_REAL &) ;
  template void closure_t::get_fluxes<1>(CCTK_REAL const&,
					 CCTK_REAL &, CCTK_REAL &,
					 CCTK_REAL &, CCTK_REAL &,
					 CCTK_REAL &) ;
  template void closure_t::get_fluxes<2>(CCTK_REAL const&,
					 CCTK_REAL &, CCTK_REAL &,
					 CCTK_REAL &, CCTK_REAL &,
					 CCTK_REAL &) ;

  template void closure_t::get_wavespeeds<0>(CCTK_REAL&, CCTK_REAL&) ;
  template void closure_t::get_wavespeeds<1>(CCTK_REAL&, CCTK_REAL&) ;
  template void closure_t::get_wavespeeds<2>(CCTK_REAL&, CCTK_REAL&) ;
  
}
