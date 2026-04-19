#ifndef HH_PEGASUS
#define HH_PEGASUS

#include <cmath>
#include <limits>
#include <tuple>

#ifdef DBG_PEGASUS_
#include<iostream>
#endif 

namespace m1_rootfinders {

template <typename F>
inline std::tuple<double,int> zero_pegasus(double a, double b, const double t, F& f ) {

  double xroot,fa, fb, ff, tol{t} ;
  int s{0} ; 
  
  int j = 0;
  fb = f(b); fa = f(a) ;
  constexpr double macheps = std::numeric_limits<double>::epsilon();

  if ( fabs(fb) < 2. * macheps  ) {
    return std::make_tuple(b,0) ; }
  else if ( fabs(fa) < 2. * macheps ) {
    return std::make_tuple(a,0) ; }
  
  for(;;) {

    xroot = (fa * b - fb*a ) / ( fa - fb ) ;
    ff = f(xroot) ;

#ifdef DBG_PEGASUS_
    std::cout << "iter,  a, b, c, fc:\t" << j << "  " << a << " "  << b << " " << " " << xroot << " " << ff << std::endl;
    #endif 

    tol = 6. * macheps * fabs(ff) + t ;
    
    if ( fabs(b-a) < tol )
      return std::make_tuple(xroot,j) ;
    
    if ( (ff*fb) > 0 )
      {
	b  = xroot ;
	fb = ff ;
	if( s==1 )
	  fa *= fb/(fb+ff) ;
	s = 1;

      } else if ( (ff * fa) > 0. ) {
      a = xroot ;
      fa = ff ;
      if(s == 2)
	fb *= fa/(fa+ff) ;
      s = 2;
    } else {
      return std::make_tuple(xroot,j) ;
    }
    j++; 
  } 
  
}

}
#endif
