#ifndef HH_ILLINOIS
#define HH_ILLINOIS

#include <cmath>
#include <limits>
#include <tuple>

#ifdef DBG_ILLINOIS_
#include<iostream>
#endif 

namespace m1_rootfinders {
  
template <typename F>
inline std::tuple<double,int> zero_illinois(double a, double b, const double t, F& f) {

  double xroot,fa, fb, ff, tol{t} ;
  int s{0} ; 
  
  int j = 0;
  fb = f(b); fa = f(a) ;
  constexpr double macheps = std::numeric_limits<double>::epsilon();

  if ( fabs(fb) < t ) {
    return std::make_tuple(b,0) ; }
  else if ( fabs(fa) < t ) {
    return std::make_tuple(a,0) ; }
  else if ( fa * fb > 0. ) {
    if ( fabs(fb) > fabs(fa) ) {
      xroot = a ;
      ff = fa ; } else {
      xroot = b ;
      ff = fb ; }

    tol = t ;
    for(int i=0; i<2; i++) {
      tol *= 10 ;
      if ( fabs(ff) < t ) {
	return std::make_tuple(xroot,0) ;
      }
    }

    return std::make_tuple(-1,-1) ;

  }
  
  
  for(j=0;;j++) {

    xroot = (fa*b - fb*a ) / ( fa - fb ) ;
    ff = f(xroot) ;

#ifdef DBG_ILLINOIS_
    std::cout << "iter,  a, b, fa, fb, c, fc:\t" << j << "  " << a << " "  << b << " " << fa << ' ' << fb << ' ' << " " << xroot << " " << ff << std::endl;
    #endif 

    tol = 2. * macheps * fabs(ff) + t ;
    
    if ( fabs(b-a) < 2. * tol )
      return std::make_tuple(xroot,j) ;
    
    if ( (ff*fb) > 0 )
      {
	b  = xroot ;
	fb = ff ;
	if( s==1 )
	  fa *= 0.5 ;
	s = 1;

      } else if ( (ff*fa) > 0. ) {
      a = xroot ;
      fa = ff ;
      if(s == 2)
	fb *= 0.5 ;
      s = 2;
    } else {
      return std::make_tuple(xroot,j) ;
    }
    

  } 
  
}


}
#endif
