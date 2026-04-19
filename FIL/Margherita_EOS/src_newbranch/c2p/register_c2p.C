//
// This file is part of Margherita, the light-weight EOS framework
//
//  Copyright (C) 2017, Elias Roland Most
//                      <emost@th.physik.uni-frankfurt.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "../Margherita_EOS.h"
#include "../margherita_c2p.hh"

#include <algorithm>
#include <iostream>
#include <array>
#include <cassert>

using namespace Margherita_C2P_Registration;

//C2P Limits - declaration
  double Margherita_C2P_Limits::tau_min;
  double Margherita_C2P_Limits::lorentz_max ;
  double Margherita_C2P_Limits::z_max ;
  double Margherita_C2P_Limits::psi6_bh ;
  double Margherita_C2P_Limits::eps_max ;

  template<typename eos>
  double C2P_Hydro<eos>::k_max;


template<typename cold_eos>
inline void Margherita_register_eos_limits(std::array<double,MAX_NUM_LIMITS> limits){


   cold_eos::c2p_rho_atm =  std::max(cold_eos::eos_rhomin, limits[RHOB_ATM]);
   cold_eos::c2p_ye_atm  =  0.;
   cold_eos::c2p_ymu_atm =  0.;
   cold_eos::c2p_temp_atm = std::min(std::max(cold_eos::eos_tempmin, limits[TEMP_ATM]), cold_eos::eos_tempmax);

   typename cold_eos::error_type error;
   cold_eos::c2p_eps_atm = cold_eos::eps__temp_rho_yle_ymu(cold_eos::c2p_temp_atm, cold_eos::c2p_rho_atm,cold_eos::c2p_ye_atm,cold_eos::c2p_ymu_atm,error);

   cold_eos::c2p_eps_max = limits[EPS_MAX];
};

template<>
inline void Margherita_register_eos_limits<EOS_Tabulated>(std::array<double,MAX_NUM_LIMITS> limits){

      EOS_Tabulated::c2p_rho_atm =  std::min(std::max(EOS_Tabulated::eos_rhomin, limits[RHOB_ATM]), EOS_Tabulated::eos_rhomax);
      EOS_Tabulated::c2p_temp_atm =  std::min(std::max(EOS_Tabulated::eos_tempmin, limits[TEMP_ATM]), EOS_Tabulated::eos_tempmax);
      EOS_Tabulated::c2p_ye_atm =  std::min(std::max(EOS_Tabulated::eos_yemin, limits[YE_ATM]), EOS_Tabulated::eos_yemax);
      EOS_Tabulated::c2p_ymu_atm = 0.;

      if(EOS_Tabulated::atm_beta_eq){
        double yeL,epsL,ymuL;
        EOS_Tabulated::error_type error;
        EOS_Tabulated::press_eps_yle_ymu__beta_eq__rho_temp(epsL,yeL,ymuL,EOS_Tabulated::c2p_rho_atm,EOS_Tabulated::c2p_temp_atm, error);
        EOS_Tabulated::c2p_ye_atm = yeL;
      }

      typename EOS_Tabulated::error_type error;
      EOS_Tabulated::c2p_eps_atm = EOS_Tabulated::eps__temp_rho_yle_ymu(EOS_Tabulated::c2p_temp_atm, EOS_Tabulated::c2p_rho_atm,
									EOS_Tabulated::c2p_ye_atm,EOS_Tabulated::c2p_ymu_atm,error);

      EOS_Tabulated::c2p_eps_max = limits[EPS_MAX];
};


template<>
inline void Margherita_register_eos_limits<EOS_Leptonic>(std::array<double,MAX_NUM_LIMITS> limits){

      EOS_Leptonic::c2p_rho_atm  =  std::min(std::max(EOS_Leptonic::eos_rhomin,  limits[RHOB_ATM]), EOS_Leptonic::eos_rhomax);
      EOS_Leptonic::c2p_temp_atm =  std::min(std::max(EOS_Leptonic::eos_tempmin, limits[TEMP_ATM]), EOS_Leptonic::eos_tempmax);
      EOS_Leptonic::c2p_ye_atm  =  std::min(std::max(EOS_Leptonic::eos_ylemin, limits[YE_ATM]), EOS_Leptonic::eos_ylemax);
      EOS_Leptonic::c2p_ymu_atm  =  std::min(std::max(EOS_Leptonic::eos_ymumin, limits[YMU_ATM]), EOS_Leptonic::eos_ymumax);

      if(EOS_Leptonic::atm_muonic_beta_eq){
        double yleL,ymuL,epsL;
        EOS_Leptonic::error_type error;
        EOS_Leptonic::press_eps_yle_ymu__beta_eq__rho_temp(epsL, yleL, ymuL,
                                                           EOS_Leptonic::c2p_rho_atm, EOS_Leptonic::c2p_temp_atm,
                                                           error);
        EOS_Leptonic::c2p_ye_atm = yleL;
        EOS_Leptonic::c2p_ymu_atm = ymuL;

      }
      typename EOS_Leptonic::error_type error;
      EOS_Leptonic::c2p_eps_atm = EOS_Leptonic::eps__temp_rho_yle_ymu(EOS_Leptonic::c2p_temp_atm, EOS_Leptonic::c2p_rho_atm, 
                                              EOS_Leptonic::c2p_ye_atm, EOS_Leptonic::c2p_ymu_atm,error);

      EOS_Leptonic::c2p_eps_max = limits[EPS_MAX];
};

//General routine, but polytropic and tabulated need special treatments
template<typename eos>
void Margherita_register_limits(std::array<double,MAX_NUM_LIMITS> limits){


//C2P Limits
  Margherita_C2P_Limits::tau_min = limits[TAUENERGY_MIN];
  Margherita_C2P_Limits::lorentz_max = limits[LORENTZ_MAX];
  Margherita_C2P_Limits::psi6_bh = limits[PSI6_BH];
  Margherita_C2P_Limits::eps_max = limits[EPS_MAX];


  auto lorentz_sq = limits[LORENTZ_MAX]*limits[LORENTZ_MAX];
  Margherita_C2P_Limits::z_max = sqrt(lorentz_sq -1.);
  auto v_max_sq = (lorentz_sq - 1.) / lorentz_sq;
  C2P_Hydro<eos>::k_max = 2.*sqrt(v_max_sq)/(1.+ v_max_sq);

//C2P Table limits

  Margherita_register_eos_limits<eos>(limits);

  std::cout<< "This is Margherita, your EOS framework!" <<std::endl;
  std::cout<< "Please check the following values I have set for the atmosphere treatment:" <<std::endl;
  std::cout<< "RHOB_ATM: " << eos::c2p_rho_atm <<std::endl;
  std::cout<< "TEMP_ATM: " << eos::c2p_temp_atm <<std::endl;
  std::cout<< "YE_ATM: " << eos::c2p_ye_atm <<std::endl;
  std::cout<< "YMU_ATM: " << eos::c2p_ymu_atm <<std::endl;
  std::cout<< "EPS_ATM: " << eos::c2p_eps_atm<<std::endl;
  std::cout<< "EPS_MIN: " << eos::c2p_eps_min<<std::endl;
  std::cout<< "P_MAX: " << eos::c2p_press_max<<std::endl;
  std::cout<< "H_MIN: " << eos::c2p_h_min<<std::endl;
  std::cout<< "H_MAX: " << eos::c2p_h_max<<std::endl;
  std::cout<< "K_MAX: " << C2P_Hydro<eos>::k_max << std::endl;
  std::cout<< "Z_MAX: " << Margherita_C2P_Limits::z_max << std::endl;

   //abort(); 
}

#ifndef STANDALONE

extern "C"
void Margherita_Register_Atmosphere_C2P(CCTK_ARGUMENTS){
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
   
  // THE ORDER HERE MATTERS!!
  std::array<double,MAX_NUM_LIMITS> limitsL{{
								     rho_atm,
								     temp_atm,
								     ye_atm,
								     ymu_atm,
								     eps_max,
								     tauenergy_min,
								     gamma_max,
								     psi6_threshold		
								  }}; 

  //Register this here
  if(CCTK_Equals(eos_type, "tabulated")){
     EOS_Tabulated::atm_beta_eq = atmosphere_is_beta_eq;
     Margherita_register_limits<EOS_Tabulated>(limitsL);
  } else if (CCTK_Equals(eos_type, "tabulated_leptonic")){
     EOS_Leptonic::atm_muonic_beta_eq = atmosphere_is_beta_eq;
     EOS_Leptonic::atm_beta_eq = atmosphere_is_beta_eq;
     Margherita_register_limits<EOS_Leptonic>(limitsL);
  }else{
    if(CCTK_Equals(eos_type, "polytropic")){
     Margherita_register_limits<EOS_Polytropic>(limitsL);
    } else
    if(CCTK_Equals(eos_type, "tabulated_cold")){
     Margherita_register_limits<EOS_Hybrid_ColdTable>(limitsL);
    } else assert(!"EOS unknown");
  }
 
  // passed
  //std::cout<< "finishing Margherita_Register_Atmosphere_C2P " << std::endl;
}

#endif



