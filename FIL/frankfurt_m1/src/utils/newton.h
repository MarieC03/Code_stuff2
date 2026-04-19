#ifndef HH_M1_NEWTON_RAPHSON
#define HH_M1_NEWTON_RAPHSON

#include <cmath>
#include <limits>
#include <tuple>

#ifdef DEBUG_NR_
#include<iostream>
#endif 

namespace m1_rootfinders {
  template <typename T, class functor, class criterion>
  static double rootfind_constrained_newton_raphson(T const& a, T const& b, functor const& func, criterion const& stopif, unsigned long& iter) noexcept  
  {
    unsigned long const maxiter = iter ;
    
    auto fl = func(a) ; auto fh = func(b) ;
    auto xl = a; auto xh = b ;
    auto x0 = 0.5*(xh+xl);
    auto f0 = func(x0);
#ifdef DEBUG_NR_
    std::cout << f0.first << std::endl ;
#endif 
    iter = 0 ; 
    auto x1 = x0 ;
    do {
      x0 = x1 ;
      f0 = func(x1) ;
      x1 = x0 - f0.first / f0.second ;
      // if we exit the allowed range use bisection
      if ( x1 < xl || x1 > xh ) { 
	if( f0.first*fh.first > 0 ) {
	  xh = x1 ;
	  fh = f0 ;
	} else {
	  xl = x1 ;
	  fl = f0 ;
	}
      }
      
#ifdef DEBUG_NR_
      std::cout << "Iter: " << iter << '\n' ;
      std::cout << "x " << x0 << ", " << f0.first/f0.second << std::endl ;
      std::cout << "dx " << std::fabs(x0-x1) << std::endl ; 
#endif 	
      iter ++ ;
    } while( !stopif(x0,x1) && iter < maxiter  ) ;
    
    return x1 ; 
  }

}
#endif
