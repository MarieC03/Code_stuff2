#include "Margherita_M1.h"
#include "constexpr_utils.hh"

#ifndef BLACK_BODY_M1
#define BLACK_BODY_M1

template<int n> 
inline __attribute__((always_inline)) 
double black_body_mev(double const &T, double const& eta, int const g=1) {
    using namespace Margherita_constants ; 
    using namespace M1_Constants;
    static_assert( (n==0) || (n==1) ) ;
    double const bb = g * 4.*pi*clite/::int_pow<3>(hc_mevcm) * ::int_pow<3+n>(T) * Fermi_Dirac<2+n>::get(eta) ; 
    return ::max(bb, 1.e-30)  ;
} 

template< int n> 
inline __attribute__((always_inline)) 
double black_body_cactus(double const &T, double const& eta, int const g=1) {
    using namespace Margherita_constants ; 
    using namespace M1_Constants;
    static_assert( (n==0) || (n==1) ) ;
    double const bb = black_body_mev<n>(T,eta,g) * mev_to_erg * EPSGF * RHOGF / clite;
    return ::max(bb, 1.e-30)  ;
} 

#endif 
