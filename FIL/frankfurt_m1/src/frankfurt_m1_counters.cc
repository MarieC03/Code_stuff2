/**
 * @file frankfurt_m1_counters.cc
 *
 * This file is part of frankfurt_m1
 *
 *  Handle counters for RK substeps and RHS evaluation.
 *
 * @author  Carlo Musolino
 */ 
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <cstring>

extern "C" void frankfurt_m1_init_rhs_counter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (use_5_spec_m1) {
     *frankfurt_m1_rhs_species = 5 ;
  } else {
     *frankfurt_m1_rhs_species = 3 ;
  }

}

extern "C" void frankfurt_m1_decrease_rhs_counter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  #ifdef DBG_COUNTERS_
  CCTK_VInfo(CCTK_THORNSTRING,"Species counter: %d", *frankfurt_m1_rhs_species) ;
  #endif 
  (*frankfurt_m1_rhs_species) --;

}

extern "C" void frankfurt_m1_init_RK(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *RK_counter = 2 ;

  if (CCTK_IsFunctionAliased("EnableProlongating"))
  {
    EnableProlongating(0);
  }
  else
  {
    CCTK_WARN(CCTK_WARN_ALERT, "Cannot disable prolongation as function"
              " \"EnableProlongating\" is not provided by any thorn!");
  }

}

extern "C" void frankfurt_m1_decrease_RK(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  #ifdef DBG_COUNTERS_
  CCTK_VInfo(CCTK_THORNSTRING,"Species counter: %d", *frankfurt_m1_rhs_species) ;
  #endif 
  (*RK_counter) --;

  if ( !(*RK_counter) ) 
  {
    if (CCTK_IsFunctionAliased("EnableProlongating"))
    {
      EnableProlongating(1);
    }
    else
    {
      CCTK_WARN(CCTK_WARN_ALERT, "Cannot enable prolongation as function"
                " \"EnableProlongating\" is not provided by any thorn!");
    }
  }

}
