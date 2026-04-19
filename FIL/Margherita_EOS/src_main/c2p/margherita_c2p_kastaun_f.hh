////
//// This file is part of Margherita, the light-weight EOS framework
////
////  Copyright (C) 2017, Elias Roland Most
////                      <emost@th.physik.uni-frankfurt.de>
////
////  This program is free software: you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation, either version 3 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
////  You should have received a copy of the GNU General Public License
////  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "../Cold/cold_pwpoly.hh"
#include "../Cold/cold_pwpoly_implementation.hh"
//#include "../utils/ND-newton-raphson.cc"
#include "constexpr_utils.hh"
#include "../utils/rootfind.hh"

#include <type_traits>

template <typename eos>
class C2P_MHD_kastaun_f {
  friend eos;

 public:
  static constexpr bool fix_conservatives = true ;
  static constexpr int numcons = 9;
  static constexpr int numprims = 10;
  static constexpr int num_aux =  4; // v r_sq W_lorentz dmu 
  static constexpr unsigned long maxiters = 2000 ;
  typedef const metric_c metric;
  typedef std::array<double, numcons> cons;
  typedef std::array<double, numprims> prims;

  
  typedef double type;
  static constexpr double tol = 1.e-15;  // tolerance on hW
  
  static inline void update_main_primitives(const type &mu, cons *__restrict CONS, prims *__restrict PRIMS,
					    const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
					    double *__restrict aux_vars, double &W) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;

    auto const chi = 1. / ( 1. + mu * B2tilde ) ;
    auto const rsq = stilde_sq * ::int_pow<2>(chi) + mu * chi * ( 1. + chi ) * ::int_pow<2>(SdotBtilde) ;
    
    auto const q = (*CONS)[TAUENERGY] - 0.5 * B2tilde - 0.5 * ::int_pow<2>(mu)*::int_pow<2>(chi)*(B2tilde * stilde_sq-::int_pow<2>(SdotBtilde)) ;
    
    auto const v0_sq = stilde_sq / ( ::int_pow<2>(eos::c2p_h_min) + stilde_sq ) ;

    auto const vsq = std::min(::int_pow<2>(mu)*rsq,v0_sq) ;
    
    W = 1. / std::sqrt(1. - vsq) ;
    
    (*PRIMS)[RHOB] = (*CONS)[RHOSTAR] / W ;
    
    (*PRIMS)[EPS]  = W * ( q - mu*rsq ) + vsq * ::int_pow<2>(W) / ( 1. + W ) ;
    
    // Enforce table bounds
    (*PRIMS)[RHOB] = std::min( eos::eos_rhomax, std::max( eos::eos_rhomin, (*PRIMS)[RHOB] ) ) ;
    typename eos::error_type eos_error ;
    
    const auto epsrange =
      eos::eps_range__rho_ye((*PRIMS)[RHOB], (*PRIMS)[YE], eos_error);
    (*PRIMS)[EPS] = min(epsrange[1], max(epsrange[0], (*PRIMS)[EPS]));
    
    (*PRIMS)[PRESSURE] = eos::press__eps_rho_ye((*PRIMS)[EPS], (*PRIMS)[RHOB], (*PRIMS)[YE], eos_error) ;
    // handle error ! (TODO)
    
    auto const a = (*PRIMS)[PRESSURE] / ( (*PRIMS)[RHOB] * ( 1. + (*PRIMS)[EPS] ) ) ;
    
    auto const vb = (1. + a) * ( 1. + q - mu * rsq ) ;
    auto const va = (1. + a) * ( 1. + (*PRIMS)[EPS] ) / W ;

    auto const v = std::max(va,vb) ;
    aux_vars[2] = W ;
    aux_vars[0] = v; aux_vars[1] = rsq ;
    
  }

  static inline type evaluate(const type &mu, cons *__restrict CONS,
                              prims *__restrict PRIMS, const double &stilde_sq,
                              const double &B2tilde, const double &SdotBtilde,
                              double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    double W ;
    update_main_primitives(mu, CONS, PRIMS, stilde_sq, B2tilde, SdotBtilde, aux_vars,W) ;
    #ifdef DEBUG
    std::cout << "\n";
    std::cout << "PRIMS[PRESSURE] :" << (*PRIMS)[PRESSURE] << std::endl;
    std::cout << "PRIMS[EPS] :" << (*PRIMS)[EPS] << std::endl;
    std::cout << "residual :" << mu - 1. / ( aux_vars[0] + mu*aux_vars[1] ) << std::endl;
    ;
    std::cout << "mu :" << mu << std::endl;
    std::cout << "W :" << W << std::endl;
#endif
    return mu - 1. / ( aux_vars[0] + mu*aux_vars[1] ) ;
  }


  static inline void compute_braketing_interval(
      double &A, double &B, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    
    A = 0 ;
    auto const fa = [&] (double const& mu) {
		      auto const chi = 1. / ( 1. + mu * B2tilde ) ;
		      auto const rsq = stilde_sq * ::int_pow<2>(chi) + mu * chi * ( 1. + chi ) * ( ::int_pow<2>(SdotBtilde) );
		      return mu*std::sqrt(::int_pow<2>(eos::c2p_h_min) + rsq) -1. ;
		    };

    auto const dr2 = [&] (double const& mu) {
		       auto const chi = 1. / ( 1. + mu * B2tilde ) ;
		       auto const rsq = stilde_sq * ::int_pow<2>(chi) + mu * chi * ( 1. + chi ) * ( ::int_pow<2>(SdotBtilde) );
		       return -2.*rsq*B2tilde / ::int_pow<3>(B2tilde*mu +1.) + ( chi + ::int_pow<2>(chi) - mu*B2tilde/::int_pow<2>(B2tilde*mu + 1.) - 2.*mu*B2tilde/::int_pow<3>(B2tilde*mu+1.) )
			 * ::int_pow<2>(SdotBtilde)  ;
		     };

    auto const dfa = [&] (double const& mu) {
		       auto const chi = 1. / ( 1. + mu * B2tilde ) ;
		       auto const rsq = stilde_sq * ::int_pow<2>(chi) + mu * chi * ( 1. + chi ) * (  ::int_pow<2>(SdotBtilde) );
		       return std::sqrt( ::int_pow<2>(eos::c2p_h_min) + rsq)  + mu / ( 2. * std::sqrt(::int_pow<2>(eos::c2p_h_min) + rsq) ) * dr2(mu) ; 
		     };

    auto const fun = [&] (double const& mu) {
			return std::make_pair(fa(mu), dfa(mu) ) ;
		     };

    auto const stopif = [&] (double const& x0, double const& x1, double const& f)
			{
			  return std::fabs(x0-x1) < 1e-15 + std::numeric_limits<double>::epsilon() * std::fabs(x0) ;
			} ;
    
    unsigned long iters{ 15 } ;
    
    B = rootfind::rootfind_newton_raphson(0., 1./eos::c2p_h_min, fun, 1.e-15, iters ) ;
       
    if( iters == maxiters ) {
      iters *= 1e05 ;
      //B = rootfind::rootfind_brent(0.,1./eos::c2p_h_min, fa, 1.e-15, iters);
      auto b = rootfind::rootfind_illinois(0.,1./eos::c2p_h_min, fa, stopif, iters);
      B= 0.5 * ( b.first + b.second ) ;
    }
    
  }
  
  template <typename F_t>
  static inline double find_root(const double &A, const double &B, F_t &F) {
    double xm, xp; 
    unsigned long iters {maxiters} ;
    auto const stopif = [&] (double const& x0, double const& x1, double const& f)
			{
			  return std::fabs(x0-x1) < 1e-15 + std::numeric_limits<double>::epsilon() * std::fabs(x0) ;
			} ;
    std::tie(xm, xp) = rootfind::rootfind_illinois<>(A,B,F, stopif, iters);
    //return rootfind::rootfind_brent<>(A,B,F,1.e-15,iters) ;
    F.aux_vars[3] = std::fabs(xm - xp) ;
    return 0.5 * ( xm + xp ) ;
  }
  
  
  // Eq 69 in Kastaun et al 
  template <typename F_t>
  static inline double error(double &mu, F_t &F) {
    using namespace Margherita_C2P_conventions;
    /*
    const auto h = (1. + (*F.PRIMS)[EPS] + (*F.PRIMS)[PRESSURE] / (*F.PRIMS)[RHOB]) ;
    F.gamma_lorentz = 1. / mu / h ;
    return F.error_estimate(1./mu, F.gamma_lorentz);
    */
    return F.aux_vars[3] * ::int_pow<2>(F.aux_vars[2]) / mu ;
  }
};
