#ifndef ROOTFIND_HH_MARGHERITA_
#define ROOTFIND_HH_MARGHERITA_

#include "utils_rootfind.hh"

#include <stdlib.h>
#include <utility>
#include <tuple>
#include <limits>

#ifdef DEBUG_
#include <iostream>
#endif 

namespace rootfind {    

  /** Newton Raphson root polishing algorithm. 
   *  @param x0 : Initial guess, a good one goes a long way!
   *  @param a  : Low end of the bracket, feel free to set to a very high negative value for unconstrained newton raphson
   *  @param b  : High end of the bracket, feel free to set to a very high value for unconstrained newton raphson
   *  @param func  : callale, should return a std::pair<T,T> containing the function evaluation as first and its derivative as second elements respectively 
   *  @param tol : Tolerance
   *  @param iter: initially set to the max iteration count, upon return contains the number of iterations the code went through 
   *  @return : a T which results from a Newton Rapshon step s.t. xk-1, xk satisfy the stopif criterion or whatever the last computed value is if iter > maxiter 
   * Note: The function never throws. The user is responsible to check for failure by verifying that iter < maxiter.
   */ 
  template <typename T, class functor>
  static inline T rootfind_newton_raphson(T const& a, T const& b, functor const& func, T const& tol, unsigned long& iter) noexcept  
    {
      unsigned long const maxiter = iter ;

      constexpr const static T macheps = std::numeric_limits<T>::epsilon() ;
      
      auto const fa = func(a) ; auto const fb = func(b) ;

      auto x0 = ( fa.first * a - fb.first * b ) / ( fa.first - fb.first ) ;
      auto f0 = func(x0) ; 

      iter = 0 ; 
      auto x1 = x0 ;
      
      do {
	x0 = x1 ;
	f0 = func(x1) ;
	x1 = x0 - f0.first / f0.second ;

	if ( std::fabs(x0-x1) <= tol + 2. * macheps * std::fabs(x1) )
	  return x1 ;

	if ( x1 < a || x1 > b ) { 
	  if ( f0.first * fa.first < 0 )
	    x1 = ( fa.first * a - f0.first*x0 ) / ( fa.first - f0.first ) ;
	  else
	    x1 = ( fb.first * b - f0.first*x0 ) / ( fb.first -  f0.first ) ;
	}
	
	iter ++ ;
      } while(  iter < maxiter  ) ;
      
      return x1 ; 
    }


  /** Brent rootfinding. Adapted from code by John Burkardt taken from FORTRAN77 version by Richard Brent. 
   *  @param x0 : low end of the bracket
   *  @param x1 : high end of the bracket
   *  @param func: callable, can be a lambda closure, must take exactly one argumen
   *  @param t: absolute tolerance on root.
   *  @param iter: initially set to the max iteration count, upon return contains the number of iterations the code went through 
   *  @return : a std::pair<T,T> which encloses the root and satisfies the stopif criterion or signaling_NaN() if the precondition f0*f1 < 0 isn't met
   * Note: The function never throws. The user is responsible to check for failure by verifying that iter < maxiter.
   * Provided the root is initially bracketed this algorithm is guaranteed to converge. 
   */ 
  template <typename T, class function>
  static inline T rootfind_brent(T const a, T const b, function & func, const T& t, unsigned long& iter) noexcept 
  {
    T c;
    T d;
    T e;
    T fa;
    T fb;
    T fc;
    T m;
    T p;
    T q;
    T r;
    T s;
    T sa;
    T sb;
    T tol;

    T constexpr static const macheps = std::numeric_limits<T>::epsilon() ;
    
    sa = a; sb = b ;
    
    fa = func(sa); fb = func(sb) ;

    // now we begin
    unsigned long const maxiter = iter ;
    iter = 0;

    c = sa;
    fc = fa;
    e = sb-sa ;
    d = e ;
    
    do {
      if (std::fabs(fc) < std::fabs(fb)) {
	sa = sb;
	sb = c;
	c = sa;
	fa = fb;
	fb = fc;
	fc = fa;
      }

      tol = 2. * macheps * std::fabs(sb) + t ;
      m   = 0.5 * ( c - sb ) ;
      
      if ( std::fabs(m) < tol || fb == 0.0) {
	break;
      }
      
      if (std::fabs(e) < tol || std::fabs(fa) <= std::fabs(fb)) {
	e = m;
	d = e;
      } else {
	s = fb / fa;
	
	if (sa == c) {
	  p = 2.0 * m * s;
	  q = 1.0 - s;
	} else {
	  q = fa / fc;
	  r = fb / fc;
	  p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
	  q = (q - 1.0) * (r - 1.0) * (s - 1.0);
	}
	
	if (0.0 < p) {
	  q = -q;
	} else {
	  p = -p;
	}
	
	s = e;
	e = d;
	
	if (2.0 * p < 3.0 * m * q - std::fabs(tol * q) &&
	    p < std::fabs(0.5 * s * q)) {
	  d = p / q;
      } else {
	  e = m;
	  d = e;
	}
      }
      sa = sb;
      fa = fb;
      
      if (tol < std::fabs(d)) {
	sb = sb + d;
      } else if (0.0 < m) {
	sb = sb + tol;
      } else {
	sb = sb - tol;
      }
      
      fb = func(sb);
      
      if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)) {
	c = sa;
	fc = fa;
	e = sb - sa;
	d = e;
      }
      iter ++ ; 
    } while (   iter<maxiter  ) ;
    
    return sb ;
  }
  
  template <typename T, class function, class criterion>
  static std::pair<T,T> rootfind_illinois(T x0, T x1, function & func, criterion const& stopif, unsigned long& iter) noexcept
  {
    // first check whether we already have a solution
    if ( stopif(x0,x1, 1.) )
      {
        iter = 0 ;
        return std::make_pair(x0,x1) ;
      }

    auto f0 { func(x0) }, f1{ func(x1) } ;

    if ( utils::sign(f1) * utils::sign(f0) > 0 ) {
      return std::make_pair( std::numeric_limits<T>::quiet_NaN() , std::numeric_limits<T>::quiet_NaN() ) ;
    }
    int side { 0 } ;
    unsigned long const maxiter = iter ;
    iter = 0 ;
    T xroot, froot ;
    do {

      xroot = ( f0*x1 - f1*x0 ) / ( f0 - f1 ) ;
      froot = func(xroot) ;

      if ( froot * f1 > 0 ) {
        x1 = xroot ;
        f1 = froot ;
        if ( side == 1 )
          f0 *= 0.5 ;
        side = 1;
      } else if ( froot * f0 > 0 ) {
        x0 = xroot ;
        f0 = froot ;
        if ( side == 2 )
          f1 *= 0.5 ;
        side = 2 ;
      } else {
        return std::make_pair(xroot,xroot) ;
      }

      iter ++ ;
    } while( !stopif(x0,x1,froot) && iter < maxiter ) ;

    return std::make_pair(x0,x1) ;
  }
  
}
#endif 
