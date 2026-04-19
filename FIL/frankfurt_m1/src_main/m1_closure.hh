#ifndef FRANKFURT_M1_CLOSURE_HH
#define FRANKFURT_M1_CLOSURE_HH
/**
 * @file m1_closure.hh
 * @author Carlo Musolino
 *
 * M1 closure class
 * 
 * This class does all the hard work 
 * for the M1 scheme. It contains methods
 * to update the closure (obviously) and 
 * compute the pressure and the fluid frame 
 * moments, a method to compute the wavespeeds
 * and the numerical fluxes, and a routine to 
 * compute the collisional sources implicitly
 * via rootfinding. The design is heavily 
 * based on Ott et al.'s implementation in 
 * the publicly available ZelmaniM1 code.
 *
 */



#include "frankfurt_m1.h"
#include "m1_implicit_conventions.hh"
#include "utils/metric.hh"
#include "constexpr_utils.hh"

#ifndef NO_GSL__
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#endif

#ifdef STANDALONE
#define CCTK_REAL double
#endif 

namespace frankfurt_m1 {

  class closure_t {
    using vec_t = std::array<CCTK_REAL,3> ;
    using error_t = m1_implicit_conventions::m1_implicit_error_bits ;
    
    static constexpr CCTK_REAL TINY = 1.e-90 ;

    static constexpr double const Wmax = 10. ;
    static constexpr double const zmax2 = Wmax*Wmax - 1. ;
  protected : 
    metric_c * gamma ;
    // Eulerian and fluid frame
    // moments of the radiation field
    CCTK_REAL E, F, F2 ;
    CCTK_REAL E_closure ;
    vec_t F_lo, F_hi ;
    vec_t f_lo, f_hi ;
    CCTK_REAL J ;
    std::array<CCTK_REAL,3> H_lo, H_hi;
    CCTK_REAL Hn ;
    // Pressure tensor,
    // indices UP, symmetric
    std::array<CCTK_REAL,6> P;
    // Fluid velocity, Lorentz
    // factor and its powers 
    vec_t vel, vlow ;
    CCTK_REAL W, W2, W3, v2 ;
    // Contractions of F with v
    CCTK_REAL Fv, fv, Ff ;
    // Eddington factor (closure)
    CCTK_REAL zeta ;
    
    // Stuff needed for the Jacobian
    CCTK_REAL B0, Bthin, Bthick ;
    CCTK_REAL an, av, aF;
    CCTK_REAL an_thick, av_thick, aF_thick ;
    CCTK_REAL an_thin, av_thin, af_thin ;
    CCTK_REAL dthick, dthin, chi ;
    
    // H_sq helpers
    CCTK_REAL HH0, HHTT, HHtt, HHt, HHT, HHTt ;
    
  public:
    // for now this guy is public
    CCTK_REAL Gamma ;
    
    // Constructor, everything is initialised except for
    // P and the fluid frame variables.
    // These are computed when update_closure is called 
    closure_t(vec_t && __zvec, metric_c * __gamma );
    /**
     * Closure update, initialises all
     * radiation moments and computes
     * P and zeta. Uses NR or Brent
     * rootfinding for the Eddington factor
     */
    template<class closure_f>
    CCTK_REAL update_closure( CCTK_REAL const& E,
			   vec_t && Flow, CCTK_REAL const& zold, bool get_zeta) ;
    
    /*
     * Wavespeeds and fluxes for
     * the Rieman solver, called in
     * add_fluxes_and_source_terms
     */ 
    template< int flux_dirn> 
    void get_wavespeeds(CCTK_REAL &cp, CCTK_REAL &cm ) ;
    
    template< int flux_dirn> 
    void get_fluxes(CCTK_REAL const& N,
		    CCTK_REAL &Nflux, CCTK_REAL &Eflux,
		    CCTK_REAL &F_x_flux, CCTK_REAL &F_y_flux,
		    CCTK_REAL &F_z_flux ) ;
    /**
     * Solve the system of 4 equations to compute
     * implicit collisional update. Internally 
     * uses GSL for the rootfinding. 
     * @return Updates E and F_i
     * and returns sources for the fluid. Also updates
     * N and returns Y_e source
     */
#ifndef NO_GSL__
    template< class closure_f >
    CCTK_REAL get_collisional_rootfinder(std::array<CCTK_REAL, 5> const& eas,
				      CCTK_REAL &N, CCTK_REAL &E,
				      CCTK_REAL &Fx, CCTK_REAL &Fy, CCTK_REAL &Fz,
				      CCTK_REAL &dE, CCTK_REAL &dSx, CCTK_REAL &dSy, CCTK_REAL &dSz,
				      CCTK_REAL &dYe, 
				      CCTK_REAL zold,
				      CCTK_REAL const& dt, 
				      bool compute_zeta,
				      error_t& error_bits ) ;
#endif 
    // Getter for pressure tensor, needed to fill GF
    // and for source term computation
    std::array<CCTK_REAL,6> get_pressure() const { return P;  } 

    CCTK_REAL inline get_average_energy(const CCTK_REAL& N) {
      return W*(E-Fv)/N*gamma->sqrtdet ;
    }
    
    // gsl wrappers
    template< class closure_f>
    friend int get_gsl_sources( gsl_vector const* u, void * param, gsl_vector * f) ;
    template< class closure_f >
    friend int get_gsl_jacobian( gsl_vector const* u, void * param, gsl_matrix * j ) ;
    template<class closure_f >
    friend int get_gsl_sources_and_jacobian( gsl_vector const* u, void * param,
					     gsl_vector * f,
					     gsl_matrix * j ) ;
    
  private:
    // Internal functions for:
    // Closure
    template<class closure_f>
    CCTK_REAL get_xi(CCTK_REAL const& zeta) ;
    template<class closure_f>
    std::pair<CCTK_REAL,CCTK_REAL> get_xi_dxi(CCTK_REAL const & zeta );
    // Implicit sources rootfinding 
#ifndef NO_GSL__
    void get_initguess_iteration(std::array<CCTK_REAL,NUM_EAS>const& eas, CCTK_REAL const& dt);
    void get_sources( std::array<CCTK_REAL, 4>& coll_sources, std::array<CCTK_REAL, NUM_EAS>const& eas) const ;
    void get_jacobian( std::array<CCTK_REAL, 16>& jac, std::array<CCTK_REAL, NUM_EAS>const& eas) const ;
    void get_number_source(CCTK_REAL &N, CCTK_REAL const& dt, std::array<CCTK_REAL, NUM_EAS> const& eas, CCTK_REAL& dYE) const ;
    void limit_primitives(gsl_vector* x) const ;
#endif 
    
  } ;
  
  
  struct implicit_solver_param_t {
    std::array<CCTK_REAL, NUM_EAS> eas ;
    std::array<CCTK_REAL, 4> Ut_explicit ;
    CCTK_REAL dt;
    bool get_zeta ;
    closure_t * closure ;
  } ;
  
}

#endif 
