//
//
// Quick and dirty test to benchmark the inversion process
// and check it for consistency. This file is anything but clean!
//
//  Copyright (C) 2017, Elias Roland Most
//                      <emost@th.physik.uni-frankfurt.de>
//

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <random>

// Let's have some debug output
#define STANDALONE
#define COLDTABLE_SETUP

//#define DEBUG

#include "../FIL_M1_headers.h"
#include "../../3D_Table/readtable.cc"
#include "../M1.hh"

constexpr int num_rho = 400;
constexpr int num_ye = 200;
constexpr int num_temp = 250;

int main() {

  //

  std::cout << std::setprecision(15);

  // EOS
  std::string table_name =
    std::string("/Users/carlomusolino/numrel/EOS_Tables/scollapse/LS220_3335_rho391_temp163_ye66.h5");
		//"/Users/carlomusolino/numrel/EOS_Tables/EOSs/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  // Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5 LS220_3335_rho391_temp163_ye66.h5 BHB_lpEOS_rho234_temp180_ye60_version_1.02_20140422.h5 BHB_lEOS_rho234_temp180_ye60_version_1.02_20140422.h5

  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;
  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, recompute_mu_nu);

  std::cout << "Read Table" << std::endl;




  // pseudo random number generator
  std::random_device rd;  // seed
  std::mt19937 gen(rd()); // generator

  // rho uniform distrib
  std::uniform_real_distribution<> dist_rho(
      std::log(EOS_Tabulated::eos_rhomin), std::log(EOS_Tabulated::eos_rhomax));
  // ye uniform distrib
  std::uniform_real_distribution<> dist_ye(EOS_Tabulated::eos_yemin,
                                           EOS_Tabulated::eos_yemax);
  // temp uniform distrib
  std::uniform_real_distribution<> dist_temp(
      std::log(EOS_Tabulated::eos_tempmin),
      std::log(EOS_Tabulated::eos_tempmax));

  //Number of NaNs
  int nan_count = 0;



  for (int nn_rho = 0; nn_rho < num_rho; nn_rho++)
    for (int nn_ye = 0; nn_ye < num_ye; nn_ye++) 
      for (int nn_temp = 0; nn_temp < num_temp; nn_temp++) {


  double rho = exp(dist_rho(gen));
  double temp = exp(dist_temp(gen));
  double ye = dist_ye(gen);

  const auto tau_init= compute_analytic_opacity(rho);
  std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
  //std::array<double,3> tau_n {{ 1.e5, 1.e5 ,1.e5}};
  auto F= Fugacities<double>(rho,temp,ye, std::move(tau_n));
    
  auto eas = EAS(F);
  eas.calc_eas();
  

  bool nan_found=false;

  for(int i=0;i<3;++i){
    nan_found = nan_found || !std::isfinite(eas.R_cactus(i));
    nan_found = nan_found || !std::isfinite(eas.Q_cactus(i));
    nan_found = nan_found || !std::isfinite(eas.kappa_a_cactus(i));
    nan_found = nan_found || !std::isfinite(eas.kappa_s_cactus(i));
    nan_found = nan_found || !std::isfinite(eas.kappa_n_cactus(i));
  }

  nan_count += nan_found;
  if(nan_found){
//	std::cout << log(rho) << "\t" << log(temp) << "\t" << ye <<std::endl;
	}

  }

  std::cout << "NaN count: MargheritaLeak has produced " << nan_count << " NaNs." <<std::endl; 
 
  return 0;
}
