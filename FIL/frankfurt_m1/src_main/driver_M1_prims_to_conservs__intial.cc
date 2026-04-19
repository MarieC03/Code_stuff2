/**
 * @file driver_M1_prims_to_conservs__initial.cc
 *
 * This file is part of frankfurt_m1
 *
 *  Set conservatives based on initial data.
 *
 * @author  Carlo Musolino
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
#include "Margherita_M1.h"

extern "C" void frankfurt_m1_prims_to_conservs__initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;



  int levelnumber = GetRefinementLevel(cctkGH);
  CCTK_REAL * restrict velx = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,0)]);
  CCTK_REAL * restrict vely = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,1)]);
  CCTK_REAL * restrict velz = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,2)]);
  //CCTK_VInfo(CCTK_THORNSTRING,"==========================================================");
  CCTK_VInfo(CCTK_THORNSTRING,"***** Initial computation of the conservative variables *****",cctk_iteration,levelnumber);
  // CCTK_VInfo(CCTK_THORNSTRING,"==========================================================");
  CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2] ; 


  // first we calculate the atmo quantities 
#pragma omp parallel for collapse(3) 
  for(int i=0;i<cctk_lsh[0];i++) for(int j=0; j<cctk_lsh[1];j++) for(int k=0;k<cctk_lsh[2];k++){
	CCTK_INT index3d = CCTK_GFINDEX3D(cctkGH,i,j,k);
	

	metric_c gamma( {gxx[index3d], gxy[index3d], gxz[index3d], gyy[index3d], gyz[index3d], gzz[index3d]},
			{betax[index3d],betay[index3d],betaz[index3d]}, alp[index3d] ) ;

	auto const sqrtgamma = gamma.sqrtdet ;
	
	CCTK_REAL const F32 = 3.15143479496164;
	
	Nnue_star[index3d] = Nnue[index3d] * sqrtgamma ;
	Enue_star[index3d] = Enue[index3d] * sqrtgamma ;
	Ft_nue_x[index3d]  = Fnue_x[index3d] * sqrtgamma ;
	Ft_nue_y[index3d]  = Fnue_y[index3d] * sqrtgamma ;
	Ft_nue_z[index3d]  = Fnue_z[index3d] * sqrtgamma ;
	temp_nue[index3d]  = temp[index3d] ;
	eps_nue[index3d]  = F32 * temp[index3d];
	
	Nnua_star[index3d] = Nnua[index3d] * sqrtgamma ;
	Enua_star[index3d] = Enua[index3d] * sqrtgamma ;
	Ft_nua_x[index3d]  = Fnua_x[index3d] * sqrtgamma ;
	Ft_nua_y[index3d]  = Fnua_y[index3d] * sqrtgamma ;
	Ft_nua_z[index3d]  = Fnua_z[index3d] * sqrtgamma ;
	temp_nua[index3d]  = temp[index3d];
	eps_nua[index3d]   = F32*temp[index3d];
	
	Nnux_star[index3d] = Nnux[index3d] * sqrtgamma ;
	Enux_star[index3d] = Enux[index3d] * sqrtgamma ;
	Ft_nux_x[index3d]  = Fnux_x[index3d] * sqrtgamma ;
	Ft_nux_y[index3d]  = Fnux_y[index3d] * sqrtgamma ;
	Ft_nux_z[index3d]  = Fnux_z[index3d] * sqrtgamma ;
	temp_nux[index3d]  = temp[index3d];
	eps_nux[index3d]   = F32*temp[index3d];
	
	// Now this is an ugly solution to make sure that at t=0 we still get something for tau
	// In particular, tau will be set to the analytic estimate from Deaton+ 2013 (which sucks)
	// Anyway since we start with everything = 0 zeta would be 1 either way... I just hope
	// that we'll rapidly converge to a more sensible value for zeta --> tau which allows the code
	// to run well in the cold regions ! 
	zeta_e[index3d] = 0.; zeta_a[index3d] = 0.; zeta_x[index3d] = 0.;
	
	
      } // for loop 
  return;
}
