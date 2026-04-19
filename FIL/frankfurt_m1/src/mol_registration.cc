/**
 * @file mol_registration.cc
 *
 * This file is part of frankfurt_m1
 *
 *  Register evolved variables with MoL
 *
 *  NB: all variables are registered as constrained 
 *      because the evolution is handled internally 
 *      within frankfurt_m1. We just use MoL for 
 *      rotating the timelevels.
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


extern "C" void frankfurt_m1_register_vars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0, group, rhs, rhs2 ;

  // Register evolution & RHS gridfunction groups

  CCTK_VInfo(CCTK_THORNSTRING, "Hello from frankfurt_m1, your radiation thorn!");
  CCTK_VInfo(CCTK_THORNSTRING, "Atmosphere density is set to: %1.5g",rho_atm);
  CCTK_VInfo(CCTK_THORNSTRING, "Neutrino atmosphere energy is set to: %1.5g",E_atmo);
  
  /* ALL OTHER EVOLVED VARIABLES (nu*_star, Enu*, F*_x, F*_y, F*_z) */
  group = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_nue");
  rhs   = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_nue_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs) ; 

  group = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_nue_bar");
  rhs   = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_nue_bar_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs) ; 

  group = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_numu");
  rhs   = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_numu_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs) ; 

  group = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_numu_bar");
  rhs   = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_numu_bar_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs) ; 

  group = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_nux");
  rhs = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_vectors_nux_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs) ; 

  group = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_scalars");
  rhs = CCTK_GroupIndex("frankfurt_m1::m1_conservatives_scalars_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs) ; 

  ierr += MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[0]")) ;
  ierr += MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[1]")) ;
  ierr += MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[2]")) ;

  ierr += MoLRegisterConstrained(CCTK_VarIndex("HydroBase::w_lorentz")) ;
  

    // lapse, metric, curv (Note that MoL does not register SaveAndRestore variables if they are registered otherwise by other thorns, so the below registration actually happens only when using fixed spacetime.)
  if(register_ADM_SandR){
  	group = CCTK_GroupIndex("ADMBase::metric");
	ierr+=MoLRegisterSaveAndRestoreGroup(group);
  	group = CCTK_GroupIndex("ADMBase::lapse");
	ierr+=MoLRegisterSaveAndRestoreGroup(group);
  	group = CCTK_GroupIndex("ADMBase::shift");
	ierr+=MoLRegisterSaveAndRestoreGroup(group);
  	group = CCTK_GroupIndex("ADMBase::curv");
	ierr+=MoLRegisterSaveAndRestoreGroup(group);
  }

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
