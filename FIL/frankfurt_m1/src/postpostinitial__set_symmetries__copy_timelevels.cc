//-------------------------------------------------
// Stuff to run right after initial data is set up
//-------------------------------------------------
/**
 * @file postpostinitial__set_symmetries__copy_timelevels.Ccc
 *
 * This file is part of frankfurt_m1
 *
 *  Register evolved variables with MoL
 *
 * @author  Carlo Musolino
 */ 
#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "frankfurt_m1.h"

extern "C" void frankfurt_m1_PostPostInitial_Set_Symmetries__Copy_Timelevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  #pragma omp simd 
  for( int index=0; index < cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]; index++) {
	Nnue_star_p[index] = Nnue_star[index] ;
	Nnue_star_p_p[index] = Nnue_star[index] ;
	Nnue_bar_star_p[index] = Nnue_bar_star[index] ;
	Nnue_bar_star_p_p[index] = Nnue_bar_star[index] ;
	Nnumu_star_p[index] = Nnumu_star[index] ;
	Nnumu_star_p_p[index] = Nnumu_star[index] ;
	Nnumu_bar_star_p[index] = Nnumu_bar_star[index] ;
	Nnumu_bar_star_p_p[index] = Nnumu_bar_star[index] ;
	Nnux_star_p[index] = Nnux_star[index] ;
	Nnux_star_p_p[index] = Nnux_star[index] ;

	Enue_star_p[index] = Enue_star[index] ;
	Enue_star_p_p[index] = Enue_star[index] ;
	Enue_bar_star_p[index] = Enue_bar_star[index] ;
	Enue_bar_star_p_p[index] = Enue_bar_star[index] ;
	Enumu_star_p[index] = Enumu_star[index] ;
	Enumu_star_p_p[index] = Enumu_star[index] ;
	Enumu_bar_star_p[index] = Enumu_bar_star[index] ;
	Enumu_bar_star_p_p[index] = Enumu_bar_star[index] ;
	Enux_star_p[index] = Enux_star[index] ;
	Enux_star_p_p[index] = Enux_star[index] ;

	Ft_nue_x_p[index] = Ft_nue_x[index] ;
	Ft_nue_y_p[index] = Ft_nue_y[index] ;
	Ft_nue_z_p[index] = Ft_nue_z[index] ;
	Ft_nue_x_p_p[index] = Ft_nue_x[index] ;
	Ft_nue_y_p_p[index] = Ft_nue_y[index] ;
	Ft_nue_z_p_p[index] = Ft_nue_z[index] ;

	Ft_nue_bar_x_p[index] = Ft_nue_bar_x[index] ;
	Ft_nue_bar_y_p[index] = Ft_nue_bar_y[index] ;
	Ft_nue_bar_z_p[index] = Ft_nue_bar_z[index] ;
	Ft_nue_bar_x_p_p[index] = Ft_nue_bar_x[index] ;
	Ft_nue_bar_y_p_p[index] = Ft_nue_bar_y[index] ;
	Ft_nue_bar_z_p_p[index] = Ft_nue_bar_z[index] ;

        Ft_numu_x_p[index] = Ft_numu_x[index] ;
        Ft_numu_y_p[index] = Ft_numu_y[index] ;
        Ft_numu_z_p[index] = Ft_numu_z[index] ;
        Ft_numu_x_p_p[index] = Ft_numu_x[index] ;
        Ft_numu_y_p_p[index] = Ft_numu_y[index] ;
        Ft_numu_z_p_p[index] = Ft_numu_z[index] ;

        Ft_numu_bar_x_p[index] = Ft_numu_bar_x[index] ;
        Ft_numu_bar_y_p[index] = Ft_numu_bar_y[index] ;
        Ft_numu_bar_z_p[index] = Ft_numu_bar_z[index] ;
        Ft_numu_bar_x_p_p[index] = Ft_numu_bar_x[index] ;
        Ft_numu_bar_y_p_p[index] = Ft_numu_bar_y[index] ;
        Ft_numu_bar_z_p_p[index] = Ft_numu_bar_z[index] ;

	Ft_nux_x_p[index] = Ft_nux_x[index] ;
	Ft_nux_y_p[index] = Ft_nux_y[index] ;
	Ft_nux_z_p[index] = Ft_nux_z[index] ;
	Ft_nux_x_p_p[index] = Ft_nux_x[index] ;
	Ft_nux_y_p_p[index] = Ft_nux_y[index] ;
	Ft_nux_z_p_p[index] = Ft_nux_z[index] ;
	
      }
};

  

