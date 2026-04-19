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
#include <fstream>
#include <string>
#include <vector>

// Let's have some debug output
#define STANDALONE
#define COLDTABLE_SETUP

#define DEBUG

enum M1_interactions { BETA=0,
		       PLASMON_DECAY,
		       BREMS,
		       PAIR,
		       NINTERACTIONS };


#include "../../3D_Table/readtable.cc"
#include "../M1.hh"

template<size_t N>
using F = Fermi_Dirac<N> ;

template<size_t N, size_t M>
using FR = Fermi_Dirac_Ratio<N,M> ;

void inline generate_logspace(const double start, const double finish, const int numpoints, std::vector<double>& out)
{
  out.clear();
  double li = std::log10(start); double le = std::log10(finish); double delta = (le-li)/numpoints;
  double x = li;
  for(int i=0;i<numpoints;i++) {
    out.push_back(std::pow(10,x));
    x += delta;
  }			   
}
int main() {

  //
  using FG=Fugacities<double>;
  
  std::cout << std::setprecision(15);

  std::cout << "F<3>(0)  " << F<3>::get(0.0) << std::endl ;
  std::cout << "FR<3,2>(0)  " << FR<3,2>::get(0.0) << std::endl ; 
  // EOS
  std::string table_name = std::string("/home/relastro-shared/EOS/scollapse/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  //std::string("/zhome/academic/HLRS/xfp/xfpmusol/EOS_Tables/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");

  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;
  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, recompute_mu_nu);

  std::cout << "Read Table" << std::endl;


  double rho = 7.582629253507815e14*Margherita_constants::RHOGF;//7.582629253507815e14*Margherita_constants::RHOGF;
  double ye = 0.234693877551020;
  double temp = 16.852590447507527;

  rho = 0.00119667;
  temp = 29.63349846;
  ye = 0.10156642;
    
  
  std::cout << "rho: " << rho << std::endl;
  std::cout << "ye: " << ye << std::endl;

  std::cout << "temp: " << temp << std::endl;
  
  const auto tau_init= Margherita_M1_EAS::EAS<double>::compute_analytic_opacity(rho); 
  std::cout << "Tau_init: " << tau_init << std::endl ;
  
  std::array<double,3> tau_n {{0.2, 2.1, 0.2}}; //{{tau_init,tau_init,tau_init}};
  auto tau_e = tau_n;
  
  auto F= Fugacities<double>(rho,temp,ye, std::move(tau_n) );
  std::array<bool, NINTERACTIONS> wi = { true, true, true, true } ;  
  auto eas = Margherita_M1_EAS::EAS<double>(F, wi );
  eas.calc_eas();

 std::cout << "Eta_e:  " << F.eta[0] << std::endl ;
 std::cout << "Ftau_e: " << F.ftau[0] << std::endl ;



 
  return 0;
}
