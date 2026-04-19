/**
 * @file driver_M1_conserv_to_prims.cc
 *
 * This file is part of frankfurt_m1
 *
 * Fill the primitive variable GFs.
 * For radiation this is trivial, since it
 * only requires a division by sqrtgamma.
 * At this stage we also check atmosphere conditions
 * and enforce causality. We also compute the neutrino
 * number normalisation Gamma and find the primitive n.
 * Moreover symmetry BCs are applied here -- only 
 * equatorial supported for now.
 * Finally this is where we fill the eddington factor GFs.
 * 
 * @author Carlo Musolino
 */

#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <tuple>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "frankfurt_m1.h" /* Generic #define's and function prototypes */
#include "m1_closure.hh"
#include "Margherita_M1.h"

// Currently this routine computes the closure everywhere.
// This is done so that we  can store a consistent value of z
// as well as to get the factor for the neutrino number c2p
// which depends on fluid frame radiation moments.
// This however might be inefficient especially since we recompute
// the closure in add_E_extrinsic_curvature_... at driver_evaluate_M1_rhs.C

extern "C" void frankfurt_m1_conservs_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  
  using namespace frankfurt_m1 ;
  
  int levelnumber = GetRefinementLevel(cctkGH);

  //CCTK_VInfo(CCTK_THORNSTRING,"==========================================================");
  if ( verbosity>1 )
    CCTK_VInfo(CCTK_THORNSTRING,"***** Iteration %d Lev: %d, Radiation c2p *****",cctk_iteration,levelnumber);
  // CCTK_VInfo(CCTK_THORNSTRING,"==========================================================");
  CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2] ;
  
#pragma omp parallel for collapse(3) 
  for(int i=0;i<cctk_lsh[0];i++) for(int j=0; j<cctk_lsh[1];j++) for(int k=0;k<cctk_lsh[2];k++){
	
	CCTK_INT index3d = CCTK_GFINDEX3D(cctkGH,i,j,k);
	// Initialise metric object 
	metric_c gamma( {gxx[index3d], gxy[index3d], gxz[index3d], gyy[index3d], gyz[index3d], gzz[index3d]},
			{betax[index3d],betay[index3d],betaz[index3d]},alp[index3d] ) ;

	auto const sqrtgamma = gamma.sqrtdet ;
	// Initialise closure object 
	closure_t cl( {zvecx[index3d],zvecy[index3d],zvecz[index3d]}, &gamma ) ;
	
	CCTK_REAL F_sq, E_sq;

	const double rL = sqrt( x[index3d] * x[index3d] + y[index3d] * y[index3d] + z[index3d] * z[index3d] ) ;
	const double E_atm = ( rL > r_atmo_damping ) ? std::max( 1e-20 , E_atmo * (r_atmo_damping / rL)*(r_atmo_damping/rL) ) : E_atmo ;
	const double N_atm = ( rL > r_atmo_damping ) ? std::max( 1e-20 , N_atmo * (r_atmo_damping / rL)*(r_atmo_damping/rL) ) : N_atmo ;
	
	double E_thresh = E_atm * sqrtgamma * ( 1. + M1_atmo_tol ) ;
	double N_thresh = N_atm * sqrtgamma * ( 1. + M1_atmo_tol ) ;

	bool nue_c2p_failed{false}, nue_bar_c2p_failed{false}, numu_c2p_failed{false}, numu_bar_c2p_failed{false}, nux_c2p_failed{false}; 

        if (use_5_spec_m1) {
  		if(
  		   (Enue_star[index3d] > E_thresh) and
  		   (Nnue_star[index3d] > N_thresh)
  		   ) {
  		  Enue[index3d]   = Enue_star[index3d] / sqrtgamma; 
  		  Fnue_x[index3d] = Ft_nue_x[index3d] / sqrtgamma; Fnue_y[index3d] = Ft_nue_y[index3d] / sqrtgamma; Fnue_z[index3d] = Ft_nue_z[index3d] / sqrtgamma;
  		  // limit flux to physical regime
  		  std::array<CCTK_REAL,3> Fl {  Fnue_x[index3d], Fnue_y[index3d], Fnue_z[index3d] } ;
  		  F_sq = gamma.square_3_lower<0,3>(Fl);
  		  E_sq = Enue[index3d]*Enue[index3d] ;
  		  if( F_sq > E_sq ) {
  		    auto const factor = std::sqrt(E_sq/F_sq)*0.999999 ;
  		    Fnue_x[index3d] *= factor ;
  		    Fnue_y[index3d] *= factor ;
  		    Fnue_z[index3d] *= factor ;
  		  }
  		    
  		  // update the closure to store it and the N normalisation
  		  // Gamma
  		  zeta_nue[index3d] = cl.update_closure<m1_closure_f>( Enue[index3d],
  								     { Fnue_x[index3d],
  								       Fnue_y[index3d],
  								       Fnue_z[index3d] }, zeta_nue[index3d], true) ;
  		  
  		  
  		  Nnue[index3d]   = Nnue_star[index3d] / sqrtgamma / cl.Gamma ;
  		  eps_nue[index3d] = cl.get_average_energy(Nnue_star[index3d]) * Margherita_constants::mnuc_MeV ;
  	
  		  auto const P = cl.get_pressure() ;
  		  PUPnue_xx[index3d] = P[0] ; PUPnue_xy[index3d] = P[1] ; PUPnue_xz[index3d] = P[2] ;
  		  PUPnue_yy[index3d] = P[3] ; PUPnue_yz[index3d] = P[4] ; PUPnue_zz[index3d] = P[5] ;
  		  
  		} else {
  		  nue_c2p_failed = true ;
  		}
  	
  		if(
  		   (Enue_bar_star[index3d] > E_thresh ) and
  		   (Nnue_bar_star[index3d] > N_thresh ) 
  		   ) {
  	
  		  Enue_bar[index3d] = Enue_bar_star[index3d] / sqrtgamma ;
  		  Fnue_bar_x[index3d] = Ft_nue_bar_x[index3d] / sqrtgamma; Fnue_bar_y[index3d] = Ft_nue_bar_y[index3d] / sqrtgamma; Fnue_bar_z[index3d] = Ft_nue_bar_z[index3d] / sqrtgamma;
  		  // limit flux to physical regime
  		  std::array<CCTK_REAL,3> Fl {  Fnue_bar_x[index3d], Fnue_bar_y[index3d], Fnue_bar_z[index3d] } ;
  		  F_sq = gamma.square_3_lower<0,3>(Fl);
  		  E_sq = Enue_bar[index3d]*Enue_bar[index3d]  ;
  		  if( F_sq > E_sq ) {
  		    auto const factor = std::sqrt(E_sq/F_sq)*0.999999 ;
  		    Fnue_bar_x[index3d] *= factor ;
  		    Fnue_bar_y[index3d] *= factor ;
  		    Fnue_bar_z[index3d] *= factor ;
  		  }
  		    
  		  // update the closure to store it and the N normalisation
  		  // Gamma
  		  zeta_nue_bar[index3d] = cl.update_closure<m1_closure_f>( Enue_bar[index3d],
  								     { Fnue_bar_x[index3d],
  								       Fnue_bar_y[index3d],
  								       Fnue_bar_z[index3d] }, zeta_nue_bar[index3d], true) ;
  		  Nnue_bar[index3d]   = Nnue_bar_star[index3d] / sqrtgamma / cl.Gamma ;
  		  eps_nue_bar[index3d] = cl.get_average_energy(Nnue_bar_star[index3d]) * Margherita_constants::mnuc_MeV ;
  		  auto const P = cl.get_pressure() ;
  		  PUPnue_bar_xx[index3d] = P[0] ; PUPnue_bar_xy[index3d] = P[1] ; PUPnue_bar_xz[index3d] = P[2] ;
  		  PUPnue_bar_yy[index3d] = P[3] ; PUPnue_bar_yz[index3d] = P[4] ; PUPnue_bar_zz[index3d] = P[5] ;
  		} else {
  		  nue_bar_c2p_failed=true;
  		}

                if(
                   (Enumu_star[index3d] > E_thresh ) and
                   (Nnumu_star[index3d] > N_thresh )
                   ) {
                  
                  Enumu[index3d] = Enumu_star[index3d] / sqrtgamma ; 
                  Fnumu_x[index3d] = Ft_numu_x[index3d] / sqrtgamma; Fnumu_y[index3d] = Ft_numu_y[index3d] / sqrtgamma; Fnumu_z[index3d] = Ft_numu_z[index3d] / sqrtgamma;
                  // limit flux to physical regime
                  std::array<CCTK_REAL,3> Fl {  Fnumu_x[index3d], Fnumu_y[index3d], Fnumu_z[index3d] } ;
                  F_sq = gamma.square_3_lower<0,3>(Fl);
                  E_sq = Enumu[index3d]*Enumu[index3d]  ;
                  if( F_sq > E_sq ) { 
                    auto const factor = std::sqrt(E_sq/F_sq)*0.999999;
                    Fnumu_x[index3d] *= factor ;
                    Fnumu_y[index3d] *= factor ;
                    Fnumu_z[index3d] *= factor ;
                  }
                    
                  // update the closure to store it and the N normalisation
                  // Gamma
                  zeta_numu[index3d] = cl.update_closure<m1_closure_f>( Enumu[index3d],
                                                                     { Fnumu_x[index3d],
                                                                       Fnumu_y[index3d],
                                                                       Fnumu_z[index3d] }, zeta_numu[index3d], true) ;
                  Nnumu[index3d]   = Nnumu_star[index3d] / sqrtgamma / cl.Gamma ;
                  eps_numu[index3d] = cl.get_average_energy(Nnumu_star[index3d]) * Margherita_constants::mnuc_MeV ;
                  auto const P = cl.get_pressure() ;
                  PUPnumu_xx[index3d] = P[0] ; PUPnumu_xy[index3d] = P[1] ; PUPnumu_xz[index3d] = P[2] ;
                  PUPnumu_yy[index3d] = P[3] ; PUPnumu_yz[index3d] = P[4] ; PUPnumu_zz[index3d] = P[5] ;
                } else {
                  numu_c2p_failed=true;
                }

                if(
                   (Enumu_bar_star[index3d] > E_thresh ) and
                   (Nnumu_bar_star[index3d] > N_thresh )
                   ) {

                  Enumu_bar[index3d] = Enumu_bar_star[index3d] / sqrtgamma ; 
                  Fnumu_bar_x[index3d] = Ft_numu_bar_x[index3d] / sqrtgamma; Fnumu_bar_y[index3d] = Ft_numu_bar_y[index3d] / sqrtgamma; Fnumu_bar_z[index3d] = Ft_numu_bar_z[index3d] / sqrtgamma;
                  // limit flux to physical regime
                  std::array<CCTK_REAL,3> Fl {  Fnumu_bar_x[index3d], Fnumu_bar_y[index3d], Fnumu_bar_z[index3d] } ;
                  F_sq = gamma.square_3_lower<0,3>(Fl) + 1.0e-45;
                  E_sq = Enumu_bar[index3d]*Enumu_bar[index3d]  ;
                  if( F_sq > E_sq ) {
                    auto const factor = std::sqrt(E_sq/F_sq)*0.999 ;
                    Fnumu_bar_x[index3d] *= factor ;
                    Fnumu_bar_y[index3d] *= factor ;
                    Fnumu_bar_z[index3d] *= factor ;
                  }

                  // update the closure to store it and the N normalisation
                  // Gamma
                  zeta_numu_bar[index3d] = cl.update_closure<m1_closure_f>( Enumu_bar[index3d],
                                                                     { Fnumu_bar_x[index3d],
                                                                       Fnumu_bar_y[index3d],
                                                                       Fnumu_bar_z[index3d] }, zeta_numu_bar[index3d], true) ;
                  Nnumu_bar[index3d]   = Nnumu_bar_star[index3d] / sqrtgamma / cl.Gamma ;
                  eps_numu_bar[index3d] = cl.get_average_energy(Nnumu_bar_star[index3d]) * Margherita_constants::mnuc_MeV ;
                  auto const P = cl.get_pressure() ;
                  PUPnumu_bar_xx[index3d] = P[0] ; PUPnumu_bar_xy[index3d] = P[1] ; PUPnumu_bar_xz[index3d] = P[2] ;
                  PUPnumu_bar_yy[index3d] = P[3] ; PUPnumu_bar_yz[index3d] = P[4] ; PUPnumu_bar_zz[index3d] = P[5] ;
                } else {
                  numu_bar_c2p_failed=true;
                }

  		if(
  		   (Enux_star[index3d] > E_thresh  ) and
  		   (Nnux_star[index3d] > N_thresh  )
  		   ) {
  	
  		   Enux[index3d] = Enux_star[index3d] / sqrtgamma;
  		   Fnux_x[index3d] = Ft_nux_x[index3d] / sqrtgamma; Fnux_y[index3d] = Ft_nux_y[index3d] / sqrtgamma; Fnux_z[index3d] = Ft_nux_z[index3d] / sqrtgamma;
  		   // limit flux to physical regime
  		  std::array<CCTK_REAL,3> Fl {  Fnux_x[index3d], Fnux_y[index3d], Fnux_z[index3d] } ;
  		  F_sq = gamma.square_3_lower<0,3>(Fl) + 1.0e-45;
  		  E_sq = Enux[index3d]*Enux[index3d] ; 
  		  if( F_sq > E_sq ) {
  		    auto const factor = std::sqrt(E_sq/F_sq)*0.999;
  		    Fnux_x[index3d] *= factor ;
  		    Fnux_y[index3d] *= factor ;
  		    Fnux_z[index3d] *= factor ;
  		  }
  		    
  		  // update the closure to store it and the N normalisation
  		  // Gamma
  		   zeta_nux[index3d] = cl.update_closure<m1_closure_f>( Enux[index3d],
  								      { Fnux_x[index3d],
  									Fnux_y[index3d],
  									Fnux_z[index3d] }, zeta_nux[index3d], true) ;
  		   Nnux[index3d]   = Nnux_star[index3d] / sqrtgamma / cl.Gamma ;
  		   eps_nux[index3d] = cl.get_average_energy(Nnux_star[index3d]) * Margherita_constants::mnuc_MeV ;
  		   auto const P = cl.get_pressure() ;
  		   PUPnux_xx[index3d] = P[0] ; PUPnux_xy[index3d] = P[1] ; PUPnux_xz[index3d] = P[2] ;
  		   PUPnux_yy[index3d] = P[3] ; PUPnux_yz[index3d] = P[4] ; PUPnux_zz[index3d] = P[5] ;
  		} else {
  		  nux_c2p_failed = true ;
  		}
  	
  	
  		// excision 
  		if(alp[index3d] < lapse_excision){
  		  nue_c2p_failed = nue_bar_c2p_failed = numu_c2p_failed = numu_bar_c2p_failed = nux_c2p_failed = true ;
  		}
  	
  		
  		// Now we set to atmo 
  		if(nue_c2p_failed) {
  		  zeta_nue[index3d] = 0. ;
  		  Nnue[index3d]   = N_atm ;
  		  Nnue_star[index3d] = N_atm * sqrtgamma ;
  		  Enue[index3d]   = E_atm ;
  		  Enue_star[index3d] = E_atm * sqrtgamma ;
  		  Fnue_x[index3d] = 0. ; Fnue_y[index3d] = 0. ; Fnue_z[index3d] = 0. ;
  		  Ft_nue_x[index3d] = 0. ; Ft_nue_y[index3d] = 0. ; Ft_nue_z[index3d] = 0. ;
  		  eps_nue[index3d] = 0.;
  		  PUPnue_xx[index3d] = 0. ; PUPnue_xy[index3d] = 0. ; PUPnue_xz[index3d] = 0. ;
  		  PUPnue_yy[index3d] = 0. ; PUPnue_yz[index3d] = 0. ; PUPnue_zz[index3d] = 0. ;
  		}
  	
  		if(nue_bar_c2p_failed) {
  		  zeta_nue_bar[index3d] = 0. ;
  		  Nnue_bar[index3d]   = N_atm ;
  		  Nnue_bar_star[index3d] = N_atm * sqrtgamma ;
  		  Enue_bar[index3d]   = E_atm ;
  		  Enue_bar_star[index3d] = E_atm * sqrtgamma ;
  		  Fnue_bar_x[index3d] = 0. ; Fnue_bar_y[index3d] = 0. ; Fnue_bar_z[index3d] = 0. ;
  		  Ft_nue_bar_x[index3d] = 0. ; Ft_nue_bar_y[index3d] = 0. ; Ft_nue_bar_z[index3d] = 0. ;
  		  eps_nue_bar[index3d] = 0.;
  		  PUPnue_bar_xx[index3d] = 0. ; PUPnue_bar_xy[index3d] = 0. ; PUPnue_bar_xz[index3d] = 0. ;
  		  PUPnue_bar_yy[index3d] = 0. ; PUPnue_bar_yz[index3d] = 0. ; PUPnue_bar_zz[index3d] = 0. ;
  		}

                if(numu_c2p_failed) {
                  zeta_numu[index3d] = 0. ;
                  Nnumu[index3d]   = N_atm ;
                  Nnumu_star[index3d] = N_atm * sqrtgamma ;
                  Enumu[index3d]   = E_atm ;
                  Enumu_star[index3d] = E_atm * sqrtgamma ;
                  Fnumu_x[index3d] = 0. ; Fnumu_y[index3d] = 0. ; Fnumu_z[index3d] = 0. ;
                  Ft_numu_x[index3d] = 0. ; Ft_numu_y[index3d] = 0. ; Ft_numu_z[index3d] = 0. ;
                  eps_numu[index3d] = 0.;
                  PUPnumu_xx[index3d] = 0. ; PUPnumu_xy[index3d] = 0. ; PUPnumu_xz[index3d] = 0. ;
                  PUPnumu_yy[index3d] = 0. ; PUPnumu_yz[index3d] = 0. ; PUPnumu_zz[index3d] = 0. ;
                }

                if(numu_bar_c2p_failed) {
                  zeta_numu_bar[index3d] = 0. ;
                  Nnumu_bar[index3d]   = N_atm ;
                  Nnumu_bar_star[index3d] = N_atm * sqrtgamma ;
                  Enumu_bar[index3d]   = E_atm ;
                  Enumu_bar_star[index3d] = E_atm * sqrtgamma ;
                  Fnumu_bar_x[index3d] = 0. ; Fnumu_bar_y[index3d] = 0. ; Fnumu_bar_z[index3d] = 0. ;
                  Ft_numu_bar_x[index3d] = 0. ; Ft_numu_bar_y[index3d] = 0. ; Ft_numu_bar_z[index3d] = 0. ;
                  eps_numu_bar[index3d] = 0.;
                  PUPnumu_bar_xx[index3d] = 0. ; PUPnumu_bar_xy[index3d] = 0. ; PUPnumu_bar_xz[index3d] = 0. ;
                  PUPnumu_bar_yy[index3d] = 0. ; PUPnumu_bar_yz[index3d] = 0. ; PUPnumu_bar_zz[index3d] = 0. ;
                }
  	
  		if(nux_c2p_failed) {
  		  zeta_nux[index3d] = 0. ;
  		  Nnux[index3d]   = N_atm ;
  		  Nnux_star[index3d]   = N_atm * sqrtgamma ;
  		  Enux[index3d]   = E_atm ;
  		  Enux_star[index3d] = E_atm * sqrtgamma ;
  		  Fnux_x[index3d] = 0. ; Fnux_y[index3d] = 0. ; Fnux_z[index3d] = 0. ;
  		  Ft_nux_x[index3d] = 0. ; Ft_nux_y[index3d] = 0. ; Ft_nux_z[index3d] = 0. ;
  		  eps_nux[index3d] = 0.;
  		  PUPnux_xx[index3d] = 0. ; PUPnux_xy[index3d] = 0. ; PUPnux_xz[index3d] = 0. ;
  		  PUPnux_yy[index3d] = 0. ; PUPnux_yz[index3d] = 0. ; PUPnux_zz[index3d] = 0. ;
  		}	

	} else { // 3 spec

  		if(
  		   (Enue_star[index3d] > E_thresh) and
  		   (Nnue_star[index3d] > N_thresh)
  		   ) {
  		  Enue[index3d]   = Enue_star[index3d] / sqrtgamma; 
  		  Fnue_x[index3d] = Ft_nue_x[index3d] / sqrtgamma; Fnue_y[index3d] = Ft_nue_y[index3d] / sqrtgamma; Fnue_z[index3d] = Ft_nue_z[index3d] / sqrtgamma;
  		  // limit flux to physical regime
  		  std::array<CCTK_REAL,3> Fl {  Fnue_x[index3d], Fnue_y[index3d], Fnue_z[index3d] } ;
  		  F_sq = gamma.square_3_lower<0,3>(Fl) ;
  		  E_sq = Enue[index3d]*Enue[index3d] ;
  		  if( F_sq > E_sq ) {
  		    auto const factor = std::sqrt(E_sq/F_sq)*0.999999 ;
  		    Fnue_x[index3d] *= factor ;
  		    Fnue_y[index3d] *= factor ;
  		    Fnue_z[index3d] *= factor ;
  		  }
  		    
  		  // update the closure to store it and the N normalisation
  		  // Gamma
  		  zeta_nue[index3d] = cl.update_closure<m1_closure_f>( Enue[index3d],
  								     { Fnue_x[index3d],
  								       Fnue_y[index3d],
  								       Fnue_z[index3d] }, zeta_nue[index3d], true) ;
  		  
  		  
  		  Nnue[index3d]   = Nnue_star[index3d] / sqrtgamma / cl.Gamma ;
  		  eps_nue[index3d] = cl.get_average_energy(Nnue_star[index3d]) * Margherita_constants::mnuc_MeV ;
  	
  		  auto const P = cl.get_pressure() ;
  		  PUPnue_xx[index3d] = P[0] ; PUPnue_xy[index3d] = P[1] ; PUPnue_xz[index3d] = P[2] ;
  		  PUPnue_yy[index3d] = P[3] ; PUPnue_yz[index3d] = P[4] ; PUPnue_zz[index3d] = P[5] ;
  		  
  		} else {
  		  nue_c2p_failed = true ;
  		}
  	
  		if(
  		   (Enue_bar_star[index3d] > E_thresh ) and
  		   (Nnue_bar_star[index3d] > N_thresh ) 
  		   ) {
  	
  		  Enue_bar[index3d] = Enue_bar_star[index3d] / sqrtgamma ;
  		  Fnue_bar_x[index3d] = Ft_nue_bar_x[index3d] / sqrtgamma; Fnue_bar_y[index3d] = Ft_nue_bar_y[index3d] / sqrtgamma; Fnue_bar_z[index3d] = Ft_nue_bar_z[index3d] / sqrtgamma;
  		  // limit flux to physical regime
  		  std::array<CCTK_REAL,3> Fl {  Fnue_bar_x[index3d], Fnue_bar_y[index3d], Fnue_bar_z[index3d] } ;
  		  F_sq = gamma.square_3_lower<0,3>(Fl) ;
  		  E_sq = Enue_bar[index3d]*Enue_bar[index3d]  ;
  		  if( F_sq > E_sq ) {
  		    auto const factor = std::sqrt(E_sq/F_sq)*0.999999 ;
  		    Fnue_bar_x[index3d] *= factor ;
  		    Fnue_bar_y[index3d] *= factor ;
  		    Fnue_bar_z[index3d] *= factor ;
  		  }
  		    
  		  // update the closure to store it and the N normalisation
  		  // Gamma
  		  zeta_nue_bar[index3d] = cl.update_closure<m1_closure_f>( Enue_bar[index3d],
  								     { Fnue_bar_x[index3d],
  								       Fnue_bar_y[index3d],
  								       Fnue_bar_z[index3d] }, zeta_nue_bar[index3d], true) ;
  		  Nnue_bar[index3d]   = Nnue_bar_star[index3d] / sqrtgamma / cl.Gamma ;
  		  eps_nue_bar[index3d] = cl.get_average_energy(Nnue_bar_star[index3d]) * Margherita_constants::mnuc_MeV ;
  		  auto const P = cl.get_pressure() ;
  		  PUPnue_bar_xx[index3d] = P[0] ; PUPnue_bar_xy[index3d] = P[1] ; PUPnue_bar_xz[index3d] = P[2] ;
  		  PUPnue_bar_yy[index3d] = P[3] ; PUPnue_bar_yz[index3d] = P[4] ; PUPnue_bar_zz[index3d] = P[5] ;
  		} else {
  		  nue_bar_c2p_failed=true;
  		}
  	
  		if(
  		   (Enux_star[index3d] > E_thresh  ) and
  		   (Nnux_star[index3d] > N_thresh  )
  		   ) {
  	
  		   Enux[index3d] = Enux_star[index3d] / sqrtgamma;
  		   Fnux_x[index3d] = Ft_nux_x[index3d] / sqrtgamma; Fnux_y[index3d] = Ft_nux_y[index3d] / sqrtgamma; Fnux_z[index3d] = Ft_nux_z[index3d] / sqrtgamma;
  		   // limit flux to physical regime
  		  std::array<CCTK_REAL,3> Fl {  Fnux_x[index3d], Fnux_y[index3d], Fnux_z[index3d] } ;
  		  F_sq = gamma.square_3_lower<0,3>(Fl) ;
  		  E_sq = Enux[index3d]*Enux[index3d] ; 
  		  if( F_sq > E_sq ) {
  		    auto const factor = std::sqrt(E_sq/F_sq)*0.999999;
  		    Fnux_x[index3d] *= factor ;
  		    Fnux_y[index3d] *= factor ;
  		    Fnux_z[index3d] *= factor ;
  		  }
  		    
  		  // update the closure to store it and the N normalisation
  		  // Gamma
  		   zeta_nux[index3d] = cl.update_closure<m1_closure_f>( Enux[index3d],
  								      { Fnux_x[index3d],
  									Fnux_y[index3d],
  									Fnux_z[index3d] }, zeta_nux[index3d], true) ;
  		   Nnux[index3d]   = Nnux_star[index3d] / sqrtgamma / cl.Gamma ;
  		   eps_nux[index3d] = cl.get_average_energy(Nnux_star[index3d]) * Margherita_constants::mnuc_MeV ;
  		   auto const P = cl.get_pressure() ;
  		   PUPnux_xx[index3d] = P[0] ; PUPnux_xy[index3d] = P[1] ; PUPnux_xz[index3d] = P[2] ;
  		   PUPnux_yy[index3d] = P[3] ; PUPnux_yz[index3d] = P[4] ; PUPnux_zz[index3d] = P[5] ;
  		} else {
  		  nux_c2p_failed = true ;
  		}
  	
  	
  		// excision 
  		if(alp[index3d] < lapse_excision){
  		  nue_c2p_failed = nue_bar_c2p_failed = nux_c2p_failed = true ;
  		}
  	
  		
  		// Now we set to atmo 
  		if(nue_c2p_failed) {
  		  zeta_nue[index3d] = 0. ;
  		  Nnue[index3d]   = N_atm ;
  		  Nnue_star[index3d] = N_atm * sqrtgamma ;
  		  Enue[index3d]   = E_atm ;
  		  Enue_star[index3d] = E_atm * sqrtgamma ;
  		  Fnue_x[index3d] = 0. ; Fnue_y[index3d] = 0. ; Fnue_z[index3d] = 0. ;
  		  Ft_nue_x[index3d] = 0. ; Ft_nue_y[index3d] = 0. ; Ft_nue_z[index3d] = 0. ;
  		  eps_nue[index3d] = 0.;
  		  PUPnue_xx[index3d] = 0. ; PUPnue_xy[index3d] = 0. ; PUPnue_xz[index3d] = 0. ;
  		  PUPnue_yy[index3d] = 0. ; PUPnue_yz[index3d] = 0. ; PUPnue_zz[index3d] = 0. ;
  		}
  	
  		if(nue_bar_c2p_failed) {
  		  zeta_nue_bar[index3d] = 0. ;
  		  Nnue_bar[index3d]   = N_atm ;
  		  Nnue_bar_star[index3d] = N_atm * sqrtgamma ;
  		  Enue_bar[index3d]   = E_atm ;
  		  Enue_bar_star[index3d] = E_atm * sqrtgamma ;
  		  Fnue_bar_x[index3d] = 0. ; Fnue_bar_y[index3d] = 0. ; Fnue_bar_z[index3d] = 0. ;
  		  Ft_nue_bar_x[index3d] = 0. ; Ft_nue_bar_y[index3d] = 0. ; Ft_nue_bar_z[index3d] = 0. ;
  		  eps_nue_bar[index3d] = 0.;
  		  PUPnue_bar_xx[index3d] = 0. ; PUPnue_bar_xy[index3d] = 0. ; PUPnue_bar_xz[index3d] = 0. ;
  		  PUPnue_bar_yy[index3d] = 0. ; PUPnue_bar_yz[index3d] = 0. ; PUPnue_bar_zz[index3d] = 0. ;
  		}
  	
  		if(nux_c2p_failed) {
  		  zeta_nux[index3d] = 0. ;
  		  Nnux[index3d]   = N_atm ;
  		  Nnux_star[index3d]   = N_atm * sqrtgamma ;
  		  Enux[index3d]   = E_atm ;
  		  Enux_star[index3d] = E_atm * sqrtgamma ;
  		  Fnux_x[index3d] = 0. ; Fnux_y[index3d] = 0. ; Fnux_z[index3d] = 0. ;
  		  Ft_nux_x[index3d] = 0. ; Ft_nux_y[index3d] = 0. ; Ft_nux_z[index3d] = 0. ;
  		  eps_nux[index3d] = 0.;
  		  PUPnux_xx[index3d] = 0. ; PUPnux_xy[index3d] = 0. ; PUPnux_xz[index3d] = 0. ;
  		  PUPnux_yy[index3d] = 0. ; PUPnux_yz[index3d] = 0. ; PUPnux_zz[index3d] = 0. ;
  		}	

                // using 3 spec --> numu numubar zeros! 
                  zeta_numu[index3d] = 0. ;
                  Nnumu[index3d]   = N_atm ;
                  Nnumu_star[index3d]   = N_atm * sqrtgamma ;
                  Enumu[index3d]   = E_atm ;
                  Enumu_star[index3d] = E_atm * sqrtgamma ;
                  Fnumu_x[index3d] = 0. ; Fnumu_y[index3d] = 0. ; Fnumu_z[index3d] = 0. ;
                  Ft_numu_x[index3d] = 0. ; Ft_numu_y[index3d] = 0. ; Ft_numu_z[index3d] = 0. ;
                  eps_numu[index3d] = 0.;
                  PUPnumu_xx[index3d] = 0. ; PUPnumu_xy[index3d] = 0. ; PUPnumu_xz[index3d] = 0. ;
                  PUPnumu_yy[index3d] = 0. ; PUPnumu_yz[index3d] = 0. ; PUPnumu_zz[index3d] = 0. ;

                  zeta_numu_bar[index3d] = 0. ;
                  Nnumu_bar[index3d]   = N_atm ;
                  Nnumu_bar_star[index3d]   = N_atm * sqrtgamma ;
                  Enumu_bar[index3d]   = E_atm ;
                  Enumu_bar_star[index3d] = E_atm * sqrtgamma ;
                  Fnumu_bar_x[index3d] = 0. ; Fnumu_bar_y[index3d] = 0. ; Fnumu_bar_z[index3d] = 0. ;
                  Ft_numu_bar_x[index3d] = 0. ; Ft_numu_bar_y[index3d] = 0. ; Ft_numu_bar_z[index3d] = 0. ;
                  eps_numu_bar[index3d] = 0.;
                  PUPnumu_bar_xx[index3d] = 0. ; PUPnumu_bar_xy[index3d] = 0. ; PUPnumu_bar_xz[index3d] = 0. ;
                  PUPnumu_bar_yy[index3d] = 0. ; PUPnumu_bar_yz[index3d] = 0. ; PUPnumu_bar_zz[index3d] = 0. ;

        } // use_5_spec_m1
      } // for loop 
  return;
}
