/*
 * =====================================================================================
 *
 *       Filename:  M1_find_betaeq.hh
 *
 *        Version:  1.0
 *       Revision:  none
 *       Compiler:  gcc
 *
 *       Author  :  Carlo Musolino
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <tuple>
#include "M1.hh"
#include "../margherita.hh"

#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_multiroots.h"

#include <limits>

#ifndef __H_M1_ROOTFIND_
#define __H_M1_ROOTFIND_



struct fluid_state {
  CCTK_REAL rho;
  CCTK_REAL temp;
  CCTK_REAL yle; // electronic lepton fraction
  CCTK_REAL ylmu; // muonic lepton fraction
  CCTK_REAL u ;
//Harry changed from 3 to 5
  std::array<double,5> tau_n;
  std::array<bool,M1_interactions::NINTERACTIONS> which_interactions;
};


template<typename eos>
static int compute_deviation_from_betaeq_rootfinder(const gsl_vector * x, void * p, gsl_vector * f);

template<typename eos>
static std::tuple<CCTK_REAL,CCTK_REAL,CCTK_REAL,bool,size_t> compute_T_ye_ymu_betaeq(
								std::array<CCTK_REAL,2*NUMSPECIES> const& Urad, 
								std::array<CCTK_REAL,REQ_FLUID_VARS> const& Ufluid, 
								std::array<CCTK_REAL, NUMSPECIES> const & __tau_n, 
								std::array<bool,M1_interactions::NINTERACTIONS> const & which_interactions) {
  using namespace Margherita_constants ;
  using namespace Margherita_M1_EAS ;
  struct fluid_state current_fluid_state;
//Harry: missing NNUX_BAR = 9, no need:?
  constexpr const int ENUE=0, ENUE_BAR=1, ENUMU=2, ENUMU_BAR=3, ENUX=4, NNUE=5, NNUE_BAR=6, NNUMU=7, NNUMU_BAR=8;
  current_fluid_state.rho  = Ufluid[RHO] ;
  current_fluid_state.temp = Ufluid[TEMP];
  current_fluid_state.yle   = Ufluid[YE] + mnuc_Msun / Ufluid[RHO] * ( Urad[NNUE] - Urad[NNUE_BAR] )  ;
  if (eos::use_muonic_eos) {
  	current_fluid_state.ylmu   = Ufluid[YMU] + mnuc_Msun / Ufluid[RHO] * ( Urad[NNUMU] - Urad[NNUMU_BAR] )  ;
  }
  //current_fluid_state.yl   = Ufluid[YE] ;

  constexpr const double macheps = std::numeric_limits<double>::epsilon() ;
  
  typename eos::error_type error ;
  double dummy;
  if (eos::use_muonic_eos) {
  	const double eps = eos::eps__temp_rho_yle_ymu(current_fluid_state.temp, current_fluid_state.rho, current_fluid_state.yle, 
						current_fluid_state.ylmu, error);
  	current_fluid_state.u = Urad[ENUE] + Urad[ENUE_BAR] + Urad[ENUMU] + Urad[ENUMU_BAR] + Urad[ENUX] + Ufluid[RHO]*(1. + eps) ;
  } else {
  	const double eps = eos::eps__temp_rho_yle_ymu(current_fluid_state.temp, current_fluid_state.rho, current_fluid_state.yle, 
						dummy, error);
  	current_fluid_state.u = Urad[ENUE] + Urad[ENUE_BAR] + Urad[ENUX] + Ufluid[RHO]*(1. + eps) ;
  }
 
  current_fluid_state.which_interactions=which_interactions ;
  
  current_fluid_state.tau_n = __tau_n ; 
  
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  gsl_multiroot_function F ;
  F.f = &(compute_deviation_from_betaeq_rootfinder<eos>) ;
  F.n = 3 ; // if I used 2, it seg fault. the size of x and f??
  //F.n = 2 ; // what is 2 meaning??  FIXME
  F.params = static_cast<void*>(&current_fluid_state) ;

  int status;
  size_t i, iter=0 ;

  gsl_vector *x = gsl_vector_alloc(3); // ye, ymu, temp

  gsl_vector_set (x,0,Ufluid[YE]) ;
  gsl_vector_set (x,1,Ufluid[YMU]) ;
  gsl_vector_set (x,2,Ufluid[TEMP]) ;

  T = gsl_multiroot_fsolver_hybrid;
  s = gsl_multiroot_fsolver_alloc(T,3) ; // Znu
  //s = gsl_multiroot_fsolver_alloc(T,2) ; 

  gsl_error_handler_t * old_err_hand = gsl_set_error_handler_off() ; 
  
  gsl_multiroot_fsolver_set(s,&F,x) ;

  do {

    iter++ ;
    status = gsl_multiroot_fsolver_iterate(s) ;

    #ifdef DEBUG_BETAEQ
     printf ("iter = %3u x = % .3f % .3f "
          "f(x) = % .5e % .5e\n",
          iter,

// Harry:  s->x  : Teq, Yeeq, Ymueq
// Harry: s->f : function (1), (2) , (3) in perego 2019
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2));
    #endif

     if(status) break ;
    
     status = gsl_multiroot_test_delta(s->dx , s->x, M1_Constants::MY_ROOTFINDING_TOL, 6.*macheps ) ;
  } while (status == GSL_CONTINUE && iter < M1_Constants::MAX_ITER_ROOTFINDING ) ;

  
  if( status!=GSL_SUCCESS ){
    #ifdef  DEBUG_BETAEQ
    printf("Failure!\n");
#endif
    gsl_vector_free(x); gsl_multiroot_fsolver_free(s);
    return std::make_tuple(0.,0.,0,false,iter) ;
  }
  
  
  double const ye_eq = gsl_vector_get (s->x, 0);
  double const ymu_eq = gsl_vector_get (s->x, 1);
  double const temp_eq = gsl_vector_get (s->x, 2);
  
  gsl_vector_free(x); gsl_multiroot_fsolver_free(s);
  
  return std::make_tuple(ye_eq,ymu_eq,temp_eq,true,iter) ;
}


template<typename eos>
static int compute_deviation_from_betaeq_rootfinder(const gsl_vector * x, void * p, gsl_vector * f){
  using namespace Margherita_M1_EAS ;
  using FG = Fugacities<CCTK_REAL>;
  // read params
  double rho = ((struct fluid_state *) p)->rho;
  double temp = ((struct fluid_state *) p)->temp;
  double yle = ((struct fluid_state *) p)->yle;
  double ylmu = ((struct fluid_state *) p)->ylmu;
  double u = ((struct fluid_state *) p)->u;
  std::array<double,5> tau_n = ((struct fluid_state *) p)->tau_n;
  //std::array<double,3> tau_n = ((struct fluid_state *) p)->tau_n;
  std::array<bool,M1_interactions::NINTERACTIONS> which_interactions = ((struct fluid_state *) p) -> which_interactions; 
  double x0 = gsl_vector_get (x,0); // ye
  double x1 = gsl_vector_get (x,1); // ymu
  double x2 = gsl_vector_get (x,2); // temp
  // eos call to get eps
  typename eos::error_type error ;
  double dummy; 
  const double eps = eos::eps__temp_rho_yle_ymu(x2, rho, x0, x1, error);
  const double e = rho * ( 1. + eps ) ;
  // get weak rates
  auto F = Fugacities<CCTK_REAL>(rho,x2,x0,x1,tau_n) ;
  EAS<CCTK_REAL> eas{F, which_interactions} ; 
  // error handling in case eos throws 
  if(error.any()) {
    gsl_vector_set(f, 2, std::numeric_limits<double>::quiet_NaN()) ;
    gsl_vector_set(f, 1, std::numeric_limits<double>::quiet_NaN()) ;
    gsl_vector_set(f, 0, std::numeric_limits<double>::quiet_NaN()) ;
    return GSL_EINVAL;
  }


  // compute f 
  // Harry:  5spec,  nux * 2 + numu + numu_bar;  3 spec, nux * 4
  if (eos::use_muonic_eos) {
        // ! Note: Z_nu (in the notation of Perego+2019) is already undensitized ( missing a factor of m_nuc / rho )
        const double y0 = x0 + eas.get_neutrino_density<NUE,0>() - eas.get_neutrino_density<NUE_BAR,0>() - yle ;
        const double y1 = x1 + eas.get_neutrino_density<NUMU,0>() - eas.get_neutrino_density<NUMU_BAR,0>() - ylmu ;
        const double y2 = e +  ( eas.get_neutrino_density<NUE,1>() + eas.get_neutrino_density<NUE_BAR,1>() +
                                 eas.get_neutrino_density<NUMU,1>() + eas.get_neutrino_density<NUMU_BAR,1>() + 
				 2. * eas.get_neutrino_density<NUX,1>() ) - u ;
  	// set f 
  	gsl_vector_set (f,0,y0) ;
  	gsl_vector_set (f,1,y1) ;
  	gsl_vector_set (f,2,y2) ;
  } else {
  	// ! Note: Z_nu (in the notation of Perego+2019) is already undensitized ( missing a factor of m_nuc / rho )
  	const double y0 = x0 + eas.get_neutrino_density<NUE,0>() - eas.get_neutrino_density<NUE_BAR,0>() - yle ;
  	const double y1 = 0.0 ;
  	const double y2 = e +  ( eas.get_neutrino_density<NUE,1>() + 
        	            eas.get_neutrino_density<NUE_BAR,1>() + 4. * eas.get_neutrino_density<NUX,1>() ) - u ;
  	// set f 
  	gsl_vector_set (f,0,y0) ;
  	gsl_vector_set (f,1,y1) ;
  	gsl_vector_set (f,2,y2) ;
  } // use_muonic_eos


  return GSL_SUCCESS ;
}



#endif

