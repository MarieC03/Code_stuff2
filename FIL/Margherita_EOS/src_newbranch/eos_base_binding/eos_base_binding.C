
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "EOS_Base.h"
#include "../Margherita_EOS.h"

#define VERR_DEF_PARAMS __LINE__, __FILE__, CCTK_THORNSTRING

template<typename eos>
double Margherita_EOSB_Pressure(double rho, double eps)
{
  typename eos::error_type error;
  double epsL,ye,ymu;
  double temp=0.;
  return eos::press_eps_yle_ymu__beta_eq__rho_temp
    (epsL, ye, ymu, rho, temp, error );
}

template<typename eos>
double Margherita_EOSB_Eps(double rho, double press)
{
  typename eos::error_type error;
  double epsL,ye,ymu;
  double temp=0.;
  double pressL = eos::press_eps_yle_ymu__beta_eq__rho_temp
    (epsL, ye, ymu, rho, temp, error );
  return epsL;
}

template<typename eos>
double Margherita_EOSB_Rho(double press, double eps)
{
  typename eos::error_type error;
  double temp=0.;
  return eos::rho__beta_eq__press_temp(press,temp, error );
}


//These are the only two functions we need.

extern "C"
CCTK_REAL Margherita_EOSB_Pressure_Tab(CCTK_REAL rho, CCTK_REAL eps){
// return Margherita_EOSB_Pressure<EOS_Tabulated>(rho,eps);
   double epsL;
   typename Hot_Slice::error_t error;
   return Hot_Slice::press_cold_eps_cold__rho(epsL,rho,error);
}

extern "C"
CCTK_REAL Margherita_EOSB_Eps_Tab(CCTK_REAL rho, CCTK_REAL press){
 //return Margherita_EOSB_Eps<EOS_Tabulated>(rho,press);
   double epsL;
   typename Hot_Slice::error_t error;
   Hot_Slice::press_cold_eps_cold__rho(epsL,rho,error);
   return epsL;
}

extern "C"
CCTK_REAL Margherita_EOSB_Rho_Tab(CCTK_REAL eps, CCTK_REAL press){
// return Margherita_EOSB_Rho<EOS_Tabulated>(press,eps);
   typename Hot_Slice::error_t error;
   return Hot_Slice::rho__press_cold(press,error);
}

extern "C"
CCTK_REAL Margherita_EOSB_Pressure_Poly(CCTK_REAL rho, CCTK_REAL eps){
 return Margherita_EOSB_Pressure<EOS_Polytropic>(rho,eps);
}

extern "C"
CCTK_REAL Margherita_EOSB_Eps_Poly(CCTK_REAL rho, CCTK_REAL press){
 return Margherita_EOSB_Eps<EOS_Polytropic>(rho,press);
}

extern "C"
CCTK_REAL Margherita_EOSB_Rho_Poly(CCTK_REAL eps, CCTK_REAL press){
 return Margherita_EOSB_Rho<EOS_Polytropic>(press,eps);
}

extern "C"
CCTK_REAL Margherita_EOSB_Pressure_ColdTable(CCTK_REAL rho, CCTK_REAL eps){
 return Margherita_EOSB_Pressure<EOS_Hybrid_ColdTable>(rho,eps);
}

extern "C"
CCTK_REAL Margherita_EOSB_Eps_ColdTable(CCTK_REAL rho, CCTK_REAL press){
 return Margherita_EOSB_Eps<EOS_Hybrid_ColdTable>(rho,press);
}

extern "C"
CCTK_REAL Margherita_EOSB_Rho_ColdTable(CCTK_REAL eps, CCTK_REAL press){
 return Margherita_EOSB_Rho<EOS_Hybrid_ColdTable>(press,eps);
}

extern "C"
CCTK_REAL Margherita_EOSB_Pressure_LepTab(CCTK_REAL rho, CCTK_REAL eps){
// return Margherita_EOSB_Pressure<EOS_Tabulated>(rho,eps);
   double epsL;
   typename Hot_Slice::error_t error;
   return Hot_Slice::press_cold_eps_cold__rho(epsL,rho,error);
}

extern "C"
CCTK_REAL Margherita_EOSB_Eps_LepTab(CCTK_REAL rho, CCTK_REAL press){
 //return Margherita_EOSB_Eps<EOS_Tabulated>(rho,press);
   double epsL;
   typename Hot_Slice::error_t error;
   Hot_Slice::press_cold_eps_cold__rho(epsL,rho,error);
   return epsL;
}

extern "C"
CCTK_REAL Margherita_EOSB_Rho_LepTab(CCTK_REAL eps, CCTK_REAL press){
// return Margherita_EOSB_Rho<EOS_Tabulated>(press,eps);
   typename Hot_Slice::error_t error;
   return Hot_Slice::rho__press_cold(press,error);
}



extern "C"
void Margherita_EOS_Base_Binding(CCTK_ARGUMENTS){
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle = EOS_RegisterMethod("2D_Polytrope");
  if(CCTK_Equals(eos_type, "polytropic")){
    EOS_RegisterPressure( handle, Margherita_EOSB_Pressure_Poly);
    EOS_RegisterSpecificIntEnergy(handle, Margherita_EOSB_Eps_Poly);
    EOS_RegisterRestMassDens(handle, Margherita_EOSB_Rho_Poly);
  }
  else 
   if(CCTK_Equals(eos_type, "tabulated")){
    EOS_RegisterPressure( handle, Margherita_EOSB_Pressure_Tab);
    EOS_RegisterSpecificIntEnergy(handle, Margherita_EOSB_Eps_Tab);
    EOS_RegisterRestMassDens(handle, Margherita_EOSB_Rho_Tab);
   }
   else
   if(CCTK_Equals(eos_type, "tabulated_cold")){
    EOS_RegisterPressure( handle, Margherita_EOSB_Pressure_ColdTable);
    EOS_RegisterSpecificIntEnergy(handle, Margherita_EOSB_Eps_ColdTable);
    EOS_RegisterRestMassDens(handle, Margherita_EOSB_Rho_ColdTable);
   }
   else
   if(CCTK_Equals(eos_type, "tabulated_leptonic")){
    EOS_RegisterPressure( handle, Margherita_EOSB_Pressure_LepTab);
    EOS_RegisterSpecificIntEnergy(handle, Margherita_EOSB_Eps_LepTab);
    EOS_RegisterRestMassDens(handle, Margherita_EOSB_Rho_LepTab);
   }
   else
    CCTK_VError(VERR_DEF_PARAMS,"ERROR.  Problem in Margherita EOS_Base binding!\n");
}
