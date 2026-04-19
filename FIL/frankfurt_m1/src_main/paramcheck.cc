/**
 * @file paramcheck.cc
 *
 * This file is part of frankfurt_m1
 *
 * @author  Carlo Musolino
 */ 
#include "cctk.h"
#include <cstdio>
#include <string>
#include <cmath>
#include <cstddef>
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"


extern "C" void frankfurt_m1_paramcheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::vel[0]")    < 2  or 
      CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::vel[1]")    < 2  or 
      CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::vel[2]")    < 2  or 
      CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::w_lorentz") < 2  or 
      CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric")        < 2  or
      CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift")         < 2  or
      CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse")         < 2  ) {
        CCTK_ERROR("The current M1 implementation requires at least one past timelevel for Hydro and Metric.") ;
  }
}
