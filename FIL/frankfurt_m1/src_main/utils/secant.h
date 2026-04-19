#ifndef SECANT_HH
#define SECANT_HH

#include <cmath>
#include <limits>

#ifdef DBG_SECANT_
#include <iostream>
#endif

static constexpr const double TINY_ = 1e-90;
static constexpr const int maxiter = 100;

template< typename F>
  static int zero_secant(double const& x0, double const& x1, double const& t, double& xroot, F& f) {

  double xn{x1}, xnm1{x0};
  double tol, tmp ;
  
  double fn{ f(x1) }, fnm1{ f(x0) } ;

  int niter{0} ;

  static constexpr double macheps = std::numeric_limits<double>::epsilon() ;
  
  auto const fp = [] ( double const& x1, double const& x0, double const& f1, double const& f0 ) {
    return ( f1 - f0 ) / ( x1 - x0 + TINY_ ) ; } ;

  do {

    #ifdef DBG_SECANT_
    std::cout << "iter xn xnm1 fn fp\t" << xn << ' ' << xnm1 << ' ' << fn << ' ' << fp(xn,xnm1,fn,fnm1) << std::endl ;
    #endif
    
    tmp = xn ;
    xn = xnm1 - fn/fp(tmp,xnm1,fn,fnm1) ; 
    
    xnm1 = tmp ;
    fn = f(xn) ;
    fnm1 = f(xnm1) ; 
    
    if ( niter > maxiter )
      return -1 ;
    tol = 6. * macheps * fn + t ;
    niter ++ ;
  } while ( ( std::fabs(xn - xnm1) > tol ) ) ;

  xroot = xn ;
  return 0 ;

}


#endif 
