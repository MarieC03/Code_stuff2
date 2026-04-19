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
#include "../FIL_M1_headers.h"
#include "../../3D_Table/readtable.cc"
#include "../M1.hh"

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

  std::ofstream nue_file;
  std::ofstream nua_file;
  std::ofstream nux_file;

  nue_file.open("nue_weakrates.txt");
  nua_file.open("nua_weakrates.txt");
  nux_file.open("nux_weakrates.txt");

  std::array<std::ofstream*,3> outfiles {{&nue_file, &nua_file, &nux_file}};
  // EOS
  std::string table_name = std::string("/Users/carlomusolino/numrel/EOS_Tables/EOSs/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  //std::string("/zhome/academic/HLRS/xfp/xfpmusol/EOS_Tables/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");

  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;
  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, recompute_mu_nu);

  std::cout << "Read Table" << std::endl;


  double rho = 7.582629253507815e14*Margherita_constants::RHOGF;
  double ye = 0.234693877551020;
  //double temp = 0;//16.852590447507527;
  std::vector<double> temps ;
  const int  npoints = 10000;
  generate_logspace(0.1,200,npoints,temps);
  
  std::cout << "rho: " << rho << std::endl;
  std::cout << "ye: " << ye << std::endl;

  for(auto& temp : temps){
    std::cout << "temp: " << temp << std::endl;
    
    const auto tau_init= compute_analytic_opacity(rho);

    std::cout << "Tau_init: " << tau_init << std::endl ;

    std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
    auto tau_e = tau_n;
    
    auto F= Fugacities<double>(rho,temp,ye, std::move(tau_n));
    
    auto eas = EAS(F);
    eas.calc_eas();

    // Output to files
    for(int i=0;i<3;++i){
      (*outfiles[i]) << temp << '\t' ;
      (*outfiles[i]) << eas.Q_cgs(i) << '\t' ;
      (*outfiles[i]) << eas.kappa_s_cgs(i) << '\t' ;
      (*outfiles[i]) << eas.kappa_a_cgs(i) << '\t' ;
      (*outfiles[i]) << eas.eps2_fluid_cgs(i) << '\t' ;
      (*outfiles[i]) << eas.eps2_leak_cgs(i) << '\t' ;
      (*outfiles[i]) << eas.eps_cgs(i) << '\t' ;
      (*outfiles[i]) << F.eta[i] << '\t' ;
      (*outfiles[i]) << F.eta[FG::PROTON] << '\t' ;
      (*outfiles[i]) << F.eta[FG::NEUTRON] << '\t' ;
      (*outfiles[i]) << F.ftau[i] << '\t' ;
      (*outfiles[i]) << F.tau_n[i] << '\t' ;
      (*outfiles[i]) << std::endl;
    }
  
  }

  nue_file.close();
  nua_file.close();
  nux_file.close();
  
  EOS_Tabulated::error_type error;
  auto range = EOS_Tabulated::eps_range__rho_ye(rho,ye,error);
  
  std::cout <<std::endl;
  std::cout<<"rho*eps_min: "<< rho*range[0] <<std::endl;


 
  return 0;
}
