//
//
// Quick and dirty test to benchmark the inversion process
// and check it for consistency. This file is anything but clean!
//
//  Copyright (C) 2021, Carlo Musolino
//  Based on routines by Elias Roland Most
//

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

// Let's have some debug output
#define STANDALONE
#define COLDTABLE_SETUP

#define DEBUG

#include "../../3D_Table/readtable.cc"
#include "../M1.hh"


int main() {

  //

  std::cout << std::setprecision(15);

  // EOS
  std::string table_name =
  std::string("/zhome/academic/HLRS/xfp/xfpmusol/EOS_Tables/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");

  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;
  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, recompute_mu_nu);

  std::cout << "Read Table" << std::endl;


  double rho = 7.582629253507815e14*Margherita_constants::RHOGF;
  double ye = 0.0234693877551020;
  double temp = 0.1;//16.852590447507527;

  const auto tau_init= compute_analytic_opacity(rho);

  std::cout << "TAU: " << tau_init <<std::endl;

  std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
  auto tau_e = tau_n;

  auto F= Fugacities(rho,temp,ye, std::move(tau_n), std::move(tau_e));

  auto eas = EAS(std::move(F));
  eas.calc_eas()

 std::cout <<std::endl;
 std::cout << "Q_cgs" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  eas.Q_cgs(i) << " , " ;
 std::cout <<std::endl;
  std::cout << "Q_cactus" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  eas.Q_cactus(i) << " , " ;
 std::cout <<std::endl;
 
 std::cout << "kappa_s" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  eas.kappa_s_cactus(i) << " , " ;
 std::cout <<std::endl;

  std::cout << "kappa_a_cactus" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  eas.kappa_a_cgs(i) << " , " ;
 std::cout <<std::endl;
  std::cout << "eps_fluid_cgs" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  eas.eps_fluid_cgs(i) << " , " ;
 std::cout <<std::endl;
  std::cout << "eps_leak_cgs" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  eas.eps_leak_cgs(i)  << " , " ;
 std::cout <<std::endl;
  std::cout << "eps_cgs" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  eas.eps_cgs(i) << " , " ;
 std::cout <<std::endl;

EOS_Tabulated::error_type error;
auto range = EOS_Tabulated::eps_range__rho_ye(rho,ye,error);

std::cout <<std::endl;
std::cout<<"rho*eps_min: "<< rho*range[0] <<std::endl;


 
  return 0;
}
