#ifndef HH_ITP
#define HH_ITP

#include <cmath>
#include <limits>

#ifdef DBG_ITP_
#include<iostream>
#endif 

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename F>
inline double zero_itp(double a, double b, const double k1, const double k2, const int n0, const double t, F& f) {

  double xh, xf, xitp, ya,yb,yitp, r, delta, sigma, tol{t}, xt; 
  
  int const nh = ceil(log2((b-a)/2./t)) ;
  int const nmax = nh + n0  ; 
  int j = 0;
  yb = f(b); ya = f(a) ;
  constexpr double macheps = std::numeric_limits<double>::epsilon();
  
  do {

    #ifdef DBG_ITP_
    std::cout << "iter,  a, b, c, fc:\t" << j << "  " << a << " "  << b << " " << " " << (b+a)*0.5 << " " << f(0.5*(b+a)) << std::endl;
    #endif 
    
    xh = 0.5*(a+b)  ; r = t * ( 1<<(nmax-j) ) - 0.5*(b-a) ; delta = k1*pow(b-a,k2) ;

    xf = ( yb*a - ya*b ) / ( yb -ya ) ;

    sigma = sgn( xh - xf ) ;

    if( delta <= std::fabs(xh - xf) )
      xt = xf + sigma*delta ;
    else
      xt = xh ;

    if ( std::fabs( xt - xh ) <= r )
      xitp = xt ;
    else
      xitp = xh - sigma*r ;

    yitp = f(xitp) ;

    if ( yitp * ya < 0. ){
      b = xitp ; yb = yitp ;
    }  else if ( yitp * ya > 0. ) {
      a = xitp; ya = yitp ;
    } else {
      a = xitp; b = xitp ;
    }
    tol = 6.*macheps * std::fabs(yitp) + t ;
    j++ ;
  } while ( (b-a) > tol ) ;

  return 0.5 * (b+a) ;
  
}

#endif
