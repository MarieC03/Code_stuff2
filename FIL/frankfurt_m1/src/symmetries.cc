/**
 * @file driver_select_BCs.cc
 *
 * This file is part of frankfurt_m1
 *
 * Register symmetries.
 *
 * @author  Carlo Musolino
 */ 
#include "cctk.h"
#include "frankfurt_m1.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"


extern "C" void frankfurt_m1_set_symmetries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VError(VERR_DEF_PARAMS,"ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);
  
  if( verbosity ) 
    CCTK_VInfo(CCTK_THORNSTRING,"Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);
  
  int sym[3];
  int ierr=0;
  
  sym[0] = 1; sym[1] = 1; sym[2] = 1;  
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Enue_star");
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Enue_bar_star");
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Enumu_star");
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Enumu_bar_star");
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Enux_star");
  
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Nnue_star");
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Nnue_bar_star");
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Nnumu_star");
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Nnumu_bar_star");
  ierr += SetCartSymVN(cctkGH,sym,"frankfurt_m1::Nnux_star");  

  sym[0] = -1; sym[1] = 1; sym[2] = 1;
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nue_x");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nue_bar_x");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_numu_x");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_numu_bar_x");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nux_x");

  sym[0] = 1; sym[1] = -1; sym[2] = 1;
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nue_y");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nue_bar_y");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_numu_y");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_numu_bar_y");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nux_y");

  sym[0] = 1; sym[1] = 1; sym[2] = -1;
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nue_z");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nue_bar_z");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_numu_z");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_numu_bar_z");
  ierr += SetCartSymVN(cctkGH, sym,"frankfurt_m1::Ft_nux_z");

  if ( ierr ) {
    CCTK_WARN(0, "Problems registering Symmetries") ;
  }
  

}

