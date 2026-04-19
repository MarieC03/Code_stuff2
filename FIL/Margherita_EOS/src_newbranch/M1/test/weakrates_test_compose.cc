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

// Let's have some debug output
#define STANDALONE
#define COLDTABLE_SETUP

#define DEBUG

#include "../../3D_Table/readtable.cc"
#include "../leakage.hh"


int main() {

  //

  std::cout << std::setprecision(15);

  // EOS
  std::string table_name =
    std::string("/Users/emost/Downloads/eoscompose.h5");
 //   std::string("/Users/emost/Downloads/"
//		    "LS220_234r_136t_50y_analmu_20091212_SVNr26.h5");
//                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  //std::string("/Users/emost/tmp/EOS/new/eos_tables/DD2_margherita_09-May-2017.h5");
  constexpr bool use_energy_shift = true;

  EOS_Tabulated::readtable_compose(table_name.c_str());

  std::cout << "Read Table" << std::endl;


  double rho = 7.582629253507815e14*Margherita_constants::RHOGF;
  double ye = 0.234693877551020;
  double temp = 16.852590447507527;

  const auto tau_init= compute_analytic_opacity(rho);

  std::cout << "TAU: " << tau_init <<std::endl;

  std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
  auto tau_e = tau_n;

  auto F= Fugacities(rho,temp,ye, std::move(tau_n), std::move(tau_e));

  auto lk = Leakage(std::move(F));
  lk.calc_leakage();

 std::cout << "R_cgs" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  lk.R_emission_cgs(i) << " , ";
 std::cout <<std::endl;
 std::cout << "Q_cgs" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  lk.Q_emission_cgs(i) << " , " ;
 std::cout <<std::endl;

 std::cout << "R_diff" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  lk.R_diff_cgs(i) << " , ";
 std::cout <<std::endl;
 std::cout << "Q_diff" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  lk.Q_diff_cgs(i) << " , " ;
 std::cout <<std::endl;


 std::cout << "R_eff_cgs" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  lk.R_cgs(i) << " , ";
 std::cout <<std::endl;
 std::cout << "Q_eff_cgs" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  lk.Q_cgs(i) << " , " ;
 std::cout <<std::endl;


 std::cout << "Cactus" <<std::endl;
 std::cout << "R_cactus" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  lk.R_cactus(i) << " , " ; 
 std::cout <<std::endl;
 std::cout << "Q_cactus" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout << lk.Q_cactus(i) << " , " ;
 std::cout <<std::endl;

 std::cout << "Cactus normalised" <<std::endl;
 std::cout << "R_cactus" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout <<  lk.R_cactus(i)/rho << " , " ; 
 std::cout <<std::endl;
 std::cout << "Q_cactus" <<std::endl;
  for(int i=0;i<3;++i)
    std::cout << lk.Q_cactus(i)/rho << " , " ;
 std::cout <<std::endl;



EOS_Tabulated::error_type error;
auto range = EOS_Tabulated::eps_range__rho_ye(rho,ye,error);

std::cout <<std::endl;
std::cout<<"rho*eps_min: "<< rho*range[0] <<std::endl;


 
  return 0;
}
