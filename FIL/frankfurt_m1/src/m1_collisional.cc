/**
 * @file m1_collisional.cc
 *
 * This file is part of frankfurt_m1
 *
 * This file contains routines to compute
 * the collisional source terms of the M1
 * equations implicitly 
 * @author Carlo Musolino
 */
#include "m1_closure.hh"
#include "frankfurt_m1.h"

#include <limits>

namespace frankfurt_m1 {
  
  /**
   * Given microphysics, computes
   * the collisional source terms 
   * accordint to the current state  
   * of the closure
   */
  void closure_t::get_sources(std::array<CCTK_REAL,4>& coll_sources, std::array<CCTK_REAL,NUM_EAS> const& eas ) const {
    // assumes the closure is up to date !
    auto const alp = gamma->lapse ;
    auto const kappa = eas[K_A] + eas[K_S] ;
    
    coll_sources[0] = alp*W * ( eas[ETA] + eas[K_S]*J - kappa*(E-Fv) );
    for( int ii=0; ii<3; ii++){
      coll_sources[ii+1] = alp*W*(eas[ETA]-eas[K_A]*J)*vlow[ii]
	                 - alp*kappa*H_lo[ii] ;
    }
    return ;
  }
  // ==========================================================================
  /**
   * Given microphysics, computes
   * the collisional source terms 
   * derivatives according to the 
   * current state of the closure
   */
  // ==========================================================================
  void closure_t::get_jacobian( std::array<CCTK_REAL,16>& jac, std::array<CCTK_REAL,NUM_EAS> const& eas) const {

    auto const alp = gamma->lapse ;
    auto const kappa = eas[K_A] + eas[K_S] ;

    // Auxiliaries for derivatives
    // Radice (69-77)
    CCTK_REAL const JvF = 2.*W2*(-1.+dthin*E*fv/F + 2.*dthick*(W2-1.)/(1.+2*W2) ) ;
    CCTK_REAL const JfF = -2.*dthin*W2*E*::int_pow<2>(fv)/F;
    
    CCTK_REAL const HvE = W3*( -1. -dthin*::int_pow<2>(fv) + dthick*(2.*W2-3.)/(1.+2*W2) );
    CCTK_REAL const HfE = -dthin*W*fv ;

    CCTK_REAL const HdF = W*(1.-dthick*v2 - dthin*E*fv/F) ;
    CCTK_REAL const HvvF = 2.*W3 * ( 1. -dthin*E*fv/F -dthick*(v2+1./(2*W2*(1.+2*W2))) ) ;
    CCTK_REAL const HffF = 2.*dthin*W*E*fv/F ;
    CCTK_REAL const HvfF = 2.*dthin*W3*E*::int_pow<2>(fv)/F ;
    CCTK_REAL const HfvF = -dthin*W*E/F ;
			    
    // derivatives
    // Radice (65-68)
    CCTK_REAL const dJdE = W2 + dthin*W2*::int_pow<2>(fv) + dthick*(3.-2.*W2)*(W2-1.)/(1.+2.*W2) ;
    // index UP
    std::array<CCTK_REAL, 3> dJdF;
    for( int ii=0; ii<3; ii++)
      dJdF[ii] = JvF*vel[ii] + JfF*f_hi[ii] ;
    // index DOWN
    std::array<CCTK_REAL, 3> dHdE;
    for( int ii=0; ii<3; ii++)
      dHdE[ii] = HvE*vlow[ii] + HfE*f_lo[ii] ;
    // indices DOWN-UP
    std::array<CCTK_REAL, 9> dHdF;
    for( int ii=0; ii<3; ii++)
      for( int jj=0; jj<3; jj++){
	int idx = jj + 3*ii ;
	dHdF[idx] = HdF * kronecker_delta(ii,jj)
	  + HvvF * vlow[ii] * vel[jj]
	  + HffF * f_lo[ii] * f_hi[jj]
	  + HvfF * vlow[ii] * f_hi[jj]
	  + HfvF * f_lo[ii] * vel[jj] ;
      }
    // Jacobian
    // Radice (61-64)
    // J00
    jac[0] = -alp*W*( kappa - eas[K_S] * dJdE) ;
    // J0i
    for( int ii=1; ii<4; ii++)
      jac[ii] = alp*W*eas[K_S]*dJdF[ii-1] + alp*W*kappa*vel[ii-1] ;
    // Jj0
    for( int jj=1; jj<4; jj++)
      jac[4*jj] = -alp*(kappa*dHdE[jj-1] + W*eas[K_A]*dJdE*vlow[jj-1]) ;
    // Jij
    for( int ii=1; ii<4; ii++)
      for( int jj=1; jj<4; jj++) {
	int const idx_3 = (jj-1) + 3*(ii-1) ;
	int const idx_4 = jj+4*ii ;
	jac[idx_4] = -alp* ( kappa*dHdF[idx_3] + W*eas[K_A]*vlow[ii-1]*dJdF[jj-1]);
      }
  
    // all done!
    return ;
  }
  // ==========================================================================
  // Friend wrappers for GSL calls
  // they are not part of the class
  // because (*)(double) != (closure_t::*)(double)
  // ==========================================================================
  template< class closure_f >
  int get_gsl_sources(gsl_vector const* u, void* param, gsl_vector* f) {
    // read parameters 
    std::array<CCTK_REAL,NUM_EAS>const& eas = ((struct implicit_solver_param_t *) param)->eas ;
    std::array<CCTK_REAL,4>const& Ut_expl   = ((struct implicit_solver_param_t *) param)->Ut_explicit ;
    CCTK_REAL dt                            = ((struct implicit_solver_param_t *) param)->dt ;
    bool compute_zeta                       = ((struct implicit_solver_param_t *) param)->get_zeta ; 
    closure_t* cl                           = ((struct implicit_solver_param_t *) param)->closure ;

    // update closure
    cl->update_closure<closure_f>( gsl_vector_get(u,0),
				   { gsl_vector_get(u,1),
				     gsl_vector_get(u,2),
				     gsl_vector_get(u,3)} , cl->zeta, compute_zeta ) ;
    // get sources 
//Harry: ask  <4> mean which 4?   should I modify it for ymu?
    std::array<CCTK_REAL,4> coll_sources ;
    cl->get_sources(coll_sources, eas) ;
    // fill gsl vector
    // We're solving
    // U_explicit + M(U) = U 
    for( int ii=0; ii<4; ii++)
      {
      CCTK_REAL const y = ( Ut_expl[ii] + dt * coll_sources[ii] - gsl_vector_get(u,ii) );
      gsl_vector_set(f,ii,y) ;
      }
    return GSL_SUCCESS ;
  }
  // ==========================================================================
  template< class closure_f >
  int get_gsl_jacobian( gsl_vector const* u, void* param, gsl_matrix* J) {
    
    std::array<CCTK_REAL,NUM_EAS>const& eas = ((struct implicit_solver_param_t *) param)->eas ;
    CCTK_REAL dt                            = ((struct implicit_solver_param_t *) param)->dt ;
    bool compute_zeta                       = ((struct implicit_solver_param_t *) param)->get_zeta ; 
    closure_t* cl                           = ((struct implicit_solver_param_t *) param)->closure ;
    
    // update closure
    cl->update_closure<closure_f>( gsl_vector_get(u,0),
				   { gsl_vector_get(u,1),
				     gsl_vector_get(u,2),
				     gsl_vector_get(u,3)} , cl->zeta, compute_zeta ) ;
    // get jacobian
    std::array<CCTK_REAL,16> jac;
    cl->get_jacobian(jac, eas) ;
    // fill gsl matrix
    for( int ii=0; ii<4; ii++)
      for( int jj=0; jj<4; jj++)
	{
	  CCTK_REAL const j =  ( dt * jac[jj+4*ii] - kronecker_delta(ii,jj) );
	  gsl_matrix_set(J,ii,jj,j) ;
	}
    return GSL_SUCCESS ;
  }
  // ==========================================================================
  template<class closure_f>
  int get_gsl_sources_and_jacobian( gsl_vector const* u, void* param,
				    gsl_vector* f,
				    gsl_matrix* J) {
    
    std::array<CCTK_REAL,NUM_EAS>const& eas = ((struct implicit_solver_param_t *) param)->eas ;
    std::array<CCTK_REAL,4>const& Ut_expl   = ((struct implicit_solver_param_t *) param)->Ut_explicit ;
    CCTK_REAL dt                            = ((struct implicit_solver_param_t *) param)->dt ;
    bool compute_zeta                       = ((struct implicit_solver_param_t *) param)->get_zeta ; 
    closure_t* cl                           = ((struct implicit_solver_param_t *) param)->closure ;
    
    // update closure 
    cl->update_closure<closure_f>( gsl_vector_get(u,0),
				   { gsl_vector_get(u,1),
				     gsl_vector_get(u,2),
				     gsl_vector_get(u,3)} , cl->zeta, compute_zeta ) ;
    
    
    // get sources
    std::array<CCTK_REAL,4> coll_sources ;
    cl->get_sources(coll_sources, eas) ;
    // fill gsl vector
    for( int ii=0; ii<4; ii++)
      {
	CCTK_REAL const y =  ( Ut_expl[ii] + dt * coll_sources[ii] - gsl_vector_get(u,ii) );
	gsl_vector_set(f,ii,y) ;
      }
    // get jacobian
    std::array<CCTK_REAL,16> jac;
    cl->get_jacobian(jac, eas) ;
    // fill gsl matrix
    for( int ii=0; ii<4; ii++)
      for( int jj=0; jj<4; jj++)
	{
	  CCTK_REAL const j =  ( dt * jac[jj+4*ii] - kronecker_delta(ii,jj) );
	  gsl_matrix_set(J,ii,jj,j) ;
	}
    return GSL_SUCCESS ;
  }
  // ==========================================================================
  // Here N is N_star -- i.e. sqrtgamma Gamma n 
  void closure_t::get_number_source(CCTK_REAL &N, CCTK_REAL const& dt, std::array<CCTK_REAL,NUM_EAS> const& eas, 
				    CCTK_REAL& dYE, CCTK_REAL& dYMU) const {
    auto const alp = gamma->lapse ;
    auto const sqrtg = gamma->sqrtdet ;
    // the source for N decouples and
    // can be inverted analytically 
    N = ( N + dt*alp*sqrtg*eas[ETA_N]) / ( 1.+dt*alp*eas[K_N]/Gamma) ;
    CCTK_REAL const n = N/Gamma/sqrtg ;
    // get source for YE, contains sqrtgamma already!
    dYE = alp * sqrtg * ( eas[ETA_N] - eas[K_N]*n ) ;

//Harry: wtf?? ETA_N and K_N are nonzeros for any species? is it correct?
    // get source for YMU, contains sqrtgamma already!
    dYMU = alp * sqrtg * ( eas[ETA_N] - eas[K_N]*n ) ;

    return ;
  }
  // ==========================================================================
  /**
   * Compute initial guess for rootfinder
   * by solving the equations to O(v/c)
   * in the fluid frame and boosting back
   * under the assumption zeta=0 (optically
   * thick medium)
   */
  void closure_t::get_initguess_iteration(std::array<CCTK_REAL,NUM_EAS>const& eas, CCTK_REAL const& dt) {
    // DR is missing a factor of alp
    // on each dt
    auto const alp = gamma->lapse ; 
    // This implements the O(v/c) solution in the
    // fluid frame by Radice (94-95)
    J = ( J + alp*dt/W*eas[ETA] ) / ( 1. + eas[K_A]*alp*dt/W ) ;
    for( int ii=0; ii<3; ii++)
      H_lo[ii] /= ( 1. + alp*dt/W * (eas[K_A]+eas[K_S]) ) ;
    H_hi = gamma->raise_index<0,3>(H_lo) ;
    // H^\mu u_\mu = 0 --> H \dot n = - H^i v_i 
    CCTK_REAL H_n = -gamma->scalar_product(vel,H_lo);
    
    // Boost to Eulerian frame assuming
    // optically thick closure ( zeta = 0 )
    // Radice (9) and (10)
    E = J/3.*(4.*W2-1.) - 2.*H_n*W ;
    
    for( int ii=0; ii<3; ii++)
      F_lo[ii] = W*H_lo[ii] + ( 4./3.*W2*J - W*H_n) * vlow[ii];
    
    return ;
  }
  // ==========================================================================
  //          -- Main interface for collisional sources --
  //          called in driver_evaluate_M1_collisional_rhs.C
  // ==========================================================================
  /**
   * Takes the microphysics ( Q k_a k_s R k_n )
   * and the radiation state and solves the 
   * equations for the implicit source terms 
   * via Newton Raphson rootfinding. We're 
   * currently using GSL's multiroot interface
   * @return updated eulerian frame moments and
   * source terms 
   */
  template<class closure_f>
  CCTK_REAL closure_t::get_collisional_rootfinder( std::array<CCTK_REAL, 5> const& eas,
						   CCTK_REAL &_N, CCTK_REAL &_E,
						   CCTK_REAL &_Fx, CCTK_REAL &_Fy, CCTK_REAL &_Fz,
						   CCTK_REAL &dE, CCTK_REAL &dSx, CCTK_REAL &dSy, CCTK_REAL &dSz,
						   CCTK_REAL &dYe,  CCTK_REAL &dYmu,
						   CCTK_REAL zold,
						   CCTK_REAL const& dt,
						   bool compute_zeta,
						   m1_implicit_conventions::m1_implicit_error_bits& error_bits ){
    using namespace m1_implicit_conventions ;
    
    static constexpr const double macheps = std::numeric_limits<double>::epsilon() ;
    
    int status ;
    size_t iter = 0;
    size_t const n = 4 ;
    
    zeta = zold ;
    
    
    update_closure<closure_f>( _E,
		                         { _Fx,
		                           _Fy,
		                           _Fz } , zeta, compute_zeta ) ;
    
    implicit_solver_param_t parameters ;
    parameters.eas         = eas ;
    parameters.Ut_explicit = {E, F_lo[0], F_lo[1], F_lo[2]} ;
    parameters.dt          = dt ;
    parameters.get_zeta    = compute_zeta ;
    parameters.closure     = this ;
    
    
    gsl_vector *x = gsl_vector_alloc(n) ;
    // set initial conditions for rootfinder  
    get_initguess_iteration(eas, dt)  ;
    
    // gsl iterates once when the solver is set...
    gsl_vector_set(x,0,E); gsl_vector_set(x,1,F_lo[0]); gsl_vector_set(x,2,F_lo[1]); gsl_vector_set(x,3,F_lo[2]); 
    
    gsl_multiroot_function_fdf f = { &(get_gsl_sources<closure_f>), &(get_gsl_jacobian<closure_f>),
				     &(get_gsl_sources_and_jacobian<closure_f>), n, reinterpret_cast<void*>(&parameters) } ;
    
    // hybridsj and hybridj both perform comparatively well
    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridsj;
    gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, n) ;
    gsl_multiroot_fdfsolver_set(s, &f, x) ;
    
    
    // prevent GSL from aborting the whole run if it does not converge
    gsl_error_handler_t* old_handler = gsl_set_error_handler_off() ;
    
    do {
      iter ++ ;
      // when the sources are computed the internal state of
      // the closure is set, and the primitives are limited...
      // Therefore to keep things consistent we feed the
      // limited radiation moments back into the solver
      gsl_vector_set(x,0,E); gsl_vector_set(x,1,F_lo[0]); gsl_vector_set(x,2,F_lo[1]); gsl_vector_set(x,3,F_lo[2]); 
      status = gsl_multiroot_fdfsolver_iterate(s) ;
      if(status)
        break ;
#ifdef DEBUG_ROOTFINDING
      printf ("iter = %3u | x = % .6g % .6g % .6g  % .6g | "
              "err_f(x) = % .8e % .8e % .8e % .8e \n"
	      "zeta = %1.5f\n\n",	      
              iter,
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              gsl_vector_get (s->x, 2),
              gsl_vector_get (s->x, 3),
              gsl_vector_get (s->f, 0),
              gsl_vector_get (s->f, 1),
              gsl_vector_get (s->f, 2),
              gsl_vector_get (s->f, 3),zeta);
#endif
      // the constants are defined in m1_implicit_conventions.hh
      //status = gsl_multiroot_test_residual(s->f, SOURCES_TOL) ;
      status = gsl_multiroot_test_delta(s->dx,s->x,SOURCES_TOL,2.*macheps);
    } while( status==GSL_CONTINUE && iter < 100 );
    
#ifdef DEBUG_ROOTFINDING
    std::cout << "Status = " << status << std::endl ;
    if( ! status ) {
      printf (" Success!  x = % .6g % .6g % .6g  % .6g "
	      "f(x) = % .8e % .8e % .8e % .8e \n\n",
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              gsl_vector_get (s->x, 2),
              gsl_vector_get (s->x, 3),
              gsl_vector_get (s->f, 0),
              gsl_vector_get (s->f, 1),
              gsl_vector_get (s->f, 2),
              gsl_vector_get (s->f, 3));
     }
#endif

    auto const sqrtg = gamma->sqrtdet ;
    /*
    if (status == GSL_ENOPROG || status == GSL_ENOPROGJ ) 
      {
	status = gsl_multiroot_test_residual(s->f, SOURCES_TOL_BACKUP) ;	
      }
    */
    if( status != GSL_SUCCESS ) {
      
      error_bits[ERRNOCONV] = true ;
      // just set back to initial state
      // and pray for the best
      // We may want to use the initguess (?)
      
      _E = parameters.Ut_explicit[0] ;
      _Fx = parameters.Ut_explicit[1];
      _Fy = parameters.Ut_explicit[2];
      _Fz = parameters.Ut_explicit[3] ;
      dE = 0.; dSx = 0.; dSy=0.; dSz=0.; dYe=0.; dYmu=0.;
      
      gsl_multiroot_fdfsolver_free(s) ;
      gsl_vector_free(x) ;
      
      return 0;
    }
    
    // update closure one more time, make sure everything is consistent 
    update_closure<closure_f>( gsl_vector_get(s->x,0),
                             { gsl_vector_get(s->x,1),
                               gsl_vector_get(s->x,2),
                               gsl_vector_get(s->x,3)} , zeta, compute_zeta ) ;
    
    
    std::array<CCTK_REAL,4> coll_sources ;
    // this call does NOT update the closure 
    get_sources(coll_sources, eas) ;
    // The source terms must be densitized before we return them
    dE = sqrtg*coll_sources[0]; dSx=sqrtg*coll_sources[1]; dSy=sqrtg*coll_sources[2]; dSz=sqrtg*coll_sources[3];
    
    // return the updated primitives too
    _E  = E;
    _Fx = F_lo[0];
    _Fy = F_lo[1];
    _Fz = F_lo[2];
    
    gsl_multiroot_fdfsolver_free(s) ;
    gsl_vector_free(x) ;
    
    // finally compute ye and ymu source term 
    get_number_source(_N, dt, eas, dYe, dYmu) ;
    
    return zeta ;
    
  }

  // Instantiate templates
  template CCTK_REAL closure_t::get_collisional_rootfinder<m1_closure_f>( std::array<CCTK_REAL,5> const&,
									  CCTK_REAL&, CCTK_REAL&,
									  CCTK_REAL&, CCTK_REAL&, CCTK_REAL&,
									  CCTK_REAL&, CCTK_REAL&,CCTK_REAL&, CCTK_REAL&,
									  CCTK_REAL&, CCTK_REAL&, 
									  CCTK_REAL,
									  CCTK_REAL const&,
									  bool, 
									  m1_implicit_conventions::m1_implicit_error_bits& );

  template int get_gsl_sources<m1_closure_f>( gsl_vector const*, void*, gsl_vector*);
  template int get_gsl_jacobian<m1_closure_f>( gsl_vector const*, void*, gsl_matrix*);
  template int get_gsl_sources_and_jacobian<m1_closure_f>( gsl_vector const*, \
	                                                         void*,
	                                                         gsl_vector*,
	                                                         gsl_matrix*);
							       
  
}
// ========================================================================
//                        == END OF FILE == 
// ==========================================================================
