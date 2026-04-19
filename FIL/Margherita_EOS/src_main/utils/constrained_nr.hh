#ifndef HH_CONSTRAINED_NEWTON_RAPHSON
#define HH_CONSTRAINED_NEWTON_RAPHSON

#include <cmath>
#include <limits>

static constexpr double POISON_VAL = 1e100 ;
static constexpr double TINY       = 1e-10 ;

template<typename F_t, typename dF_t>
static double zero_constrained_newton_raphson(double const& A, double const& B, double const tol, F_t& f, dF_t& f_prime, int max_iter, int& err) {

  err = 0 ;

  double xroot, xnew, xa, xb ,ff, df, fa, fb ;
  int it ;
  xa = A    ; xb = B    ;
  fa = f(A) ; fb = f(B) ;

  constexpr double macheps = std::numeric_limits<double>::epsilon();
  
  if( fabs(fa) <= tol ) {
    return A ;
  } else if ( fabs(fb) <= tol ) {
    return B ;
  } else if ( fa*fb > 0 ) {
    err = -1 ;
    return A-1 ;
  }

  xroot = ( A + B ) * 0.5 ; 
  xnew  = POISON_VAL ;

  // iterate 
  for(it=0;;++it) {
    ff = f(xroot) ; df = f_prime(xroot) ;
    if ( fabs(df) < TINY ) {
      err = 2 ;
      return A-1. ;
    }
    
    if(  ( fabs(ff) < tol )  || ( fabs(xnew-xroot) < 2.* ( tol + 2.*macheps*fabs(ff) ) ) ) {
      return xroot ;
    }

    xnew = xroot - ff/df ;

    if ( xnew > B || xnew < A ) {
      xnew = ( A + B ) * 0.5 ;
      ff   = f(xnew) ;
      if ( ff * fa > 0. ) {
	xa = xnew ;
	fa = ff;
      } else {
	xb = xnew ;
	fb = ff ;
      }
    }
    xroot = xnew ;
  }

  err = 1 ;
  return A-1 ;
}

#endif
