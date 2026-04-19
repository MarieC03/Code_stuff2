/**
 * @file driver_select_BCs.cc
 *
 * This file is part of frankfurt_m1
 *
 *  Select Boundary Conditions and SYNC.
 *
 * @author  Carlo Musolino
 */ 
#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include "frankfurt_m1.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include <string>

extern "C" void frankfurt_m1_select_BCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  std::string my_bc ;

  if( CCTK_Equals(Radiation_BC, "outflow") ) {
    my_bc = "flat" ;
  } else if ( CCTK_Equals(Radiation_BC, "frozen") ) {
    my_bc = "static" ;
  } else {
    my_bc = "none" ;
  }

  int ierr = 0;
  auto const frankfurt_m1_register_BCs_var = [&] (std::string const& vname, std::string const& bc_type)
					     {
					       ierr += Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES,
									       cctk_nghostzones[0], -1, vname.c_str(), bc_type.c_str() ) ;
					       
					     } ;

  //if (use_5_spec_m1) {
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enue","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enue_bar","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enumu","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enumu_bar","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enux","none") ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnue","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnue_bar","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnumu","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnumu_bar","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnux","none") ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_x","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_bar_x","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnumu_x","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnumu_bar_x","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnux_x","none") ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_y","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_bar_y","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnumu_y","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnumu_bar_y","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnux_y","none") ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_z","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_bar_z","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnumu_z","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnumu_bar_z","none") ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Fnux_z","none") ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enue_star", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enue_bar_star", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enumu_star", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enumu_bar_star", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Enux_star", my_bc) ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnue_star", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnue_bar_star", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnumu_star", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnumu_bar_star", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Nnux_star", my_bc) ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_x", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_bar_x", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_numu_x", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_numu_bar_x", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nux_x", my_bc) ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_y", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_bar_y", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_numu_y", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_numu_bar_y", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nux_y", my_bc) ;
         
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_z", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_bar_z", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_numu_z", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_numu_bar_z", my_bc) ;
         frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nux_z", my_bc) ;

  //} else {
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Enue","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Enue_bar","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Enux","none") ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Nnue","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Nnue_bar","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Nnux","none") ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_x","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_bar_x","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnux_x","none") ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_y","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_bar_y","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnux_y","none") ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_z","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnue_bar_z","none") ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Fnux_z","none") ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Enue_star", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Enue_bar_star", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Enux_star", my_bc) ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Nnue_star", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Nnue_bar_star", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Nnux_star", my_bc) ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_x", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_bar_x", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nux_x", my_bc) ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_y", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_bar_y", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nux_y", my_bc) ;

  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_z", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nue_bar_z", my_bc) ;
  //       frankfurt_m1_register_BCs_var("frankfurt_m1::Ft_nux_z", my_bc) ;
  //} // use_5_spec_m1

  if(ierr) {
    CCTK_WARN(0, "Problems registering BCs") ;
  }
  
}

