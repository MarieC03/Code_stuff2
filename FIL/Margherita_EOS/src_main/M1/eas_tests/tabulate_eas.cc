#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <memory>
#include <utility>

// Let's have some debug output
#define STANDALONE
#define COLDTABLE_SETUP
#define __STANDALONE_M1_EAS_TEST_
//#define _TESTING_ILEAS_FORMULA__
//#define DEBUG
#include "../FIL_M1_headers.h"
#include "../../3D_Table/readtable.cc"
#include "../M1.hh"

#include "eas_write_table.cc"


std::string const table_name = std::string("/Users/carlomusolino/numrel/EOS_Tables/EOSs/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");


std::string const Q_name = "Q_DDw_mode.dat";
std::string const ka_name = "ka_DD2_mode.dat";
std::string const ks_name = "ks_DD2_mode.dat";


double constexpr const RHO_GF = 1.61887093132742e-18;

decltype(auto) read_data_from_file(std::string const& fname) {

  std::ifstream file(fname.c_str()) ;
  
  std::string line;
  int nlines=0;
  std::vector<double> data;
  
  while( std::getline(file,line) ){
    nlines++;
    data.push_back(std::stod(line.c_str()));
  }
  return std::make_tuple(data,nlines);
}


int main(){

  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;
  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, recompute_mu_nu);

  size_t const npoints = 8000 ;

  double rhomin{1e08*RHO_GF}, rhomax{1e15*RHO_GF};
  double tempmin{1e-01}, tempmax{1e02};
  double yemin{0.01}, yemax{0.45};

  double rho{2e11*RHO_GF}, temp{4}, ye{0.33};
  //double rho{0.00128},temp{0.1},ye{0.02};

  //double rho=1e-05;
  //double temp=40;
  //double ye=5e-02;
  
  {
    std::ofstream kappa_s{"kappa_s_rho.dat"};
    std::ofstream kappa_a{"kappa_a_rho.dat"};
    std::ofstream kPP{"kappa_pairprocs_rho.dat"};
    std::ofstream kNNB{"kappa_brems_rho.dat"};
    std::ofstream kPD("kappa_plasmon_decay_rho.dat"); 
    std::ofstream kappa_abs{"kappa_charged_currents_rho.dat"};
    std::ofstream BB{"BB.dat"};
    

    
    double h{(rhomax-rhomin)/npoints};
    for(int i=0; i<npoints; i++) {
      
      double const rhoL  = rhomin + i*h;
      double const tempL = temp;
      double const yeL   = ye;

      if(!(i%100) )
	std::cout << rhoL/RHO_GF << std::endl ;
      // calc eas with Margherita at yeL tempL rhoL 
      const auto tau_init= compute_analytic_opacity(rhoL);
      std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
      auto F= Fugacities<double>(rhoL,tempL,yeL, std::move(tau_n));
      auto eas = EAS<double>(F);
      eas.calc_eas();
      
      //save data

      kappa_a << rhoL/RHO_GF << ' ' << eas.kappa_a_cgs(NUE) << " " << eas.kappa_a_cgs(NUA) << " " << eas.kappa_a_cgs(NUX) << "\n";
      kappa_s << rhoL/RHO_GF << ' ' << eas.kappa_s_cgs(NUE) << " " << eas.kappa_s_cgs(NUA) << " " << eas.kappa_s_cgs(NUX) << "\n";
      BB << eas.BB_cactus(NUE) << " " << eas.BB_cactus(NUA) << " " << eas.BB_cactus(NUX) << '\n';
      kNNB << rhoL/RHO_GF << ' ' << eas.kappa_a_brems_cgs(NUE) << ' ' << eas.kappa_a_brems_cgs(NUA) << ' ' << eas.kappa_a_brems_cgs(NUX) << '\n';
      kPP << rhoL/RHO_GF << ' ' << eas.kappa_a_epannihil_cgs(NUE) << ' ' << eas.kappa_a_epannihil_cgs(NUA) << ' ' << eas.kappa_a_epannihil_cgs(NUX) << '\n';
      kPD << rhoL/RHO_GF << ' ' << eas.kappa_a_plasmon_cgs(NUE) << ' ' << eas.kappa_a_plasmon_cgs(NUA) << ' ' << eas.kappa_a_plasmon_cgs(NUX) << '\n' ;
      kappa_abs << rhoL/RHO_GF << ' ' << eas.kappa_abs_cgs(NUE) << ' ' << eas.kappa_abs_cgs(NUA) << ' ' << eas.kappa_abs_cgs(NUX) << '\n';
    }
  }
  {
    std::ofstream kappa_s{"kappa_s_temp.dat"};
    std::ofstream kappa_a{"kappa_a_temp.dat"};
    std::ofstream kPP{"kappa_pairprocs_temp.dat"};
    std::ofstream kNNB{"kappa_brems_temp.dat"};
    std::ofstream kPD("kappa_plasmon_decay_temp.dat"); 
    std::ofstream kappa_abs{"kappa_charged_currents_temp.dat"};
    std::ofstream BB{"BB.dat"};
    
    
    
    double h{(tempmax-tempmin)/npoints};

    for(int i=0; i<npoints; i++) {
      
      double const rhoL  = rho;
      double const tempL = tempmin + i*h;
      double const yeL   = ye;
      
      // calc eas with Margherita at yeL tempL rhoL 
      const auto tau_init= compute_analytic_opacity(rhoL);
      std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
      auto F= Fugacities<double>(rhoL,tempL,yeL, std::move(tau_n));
      auto eas = EAS<double>(F);
      eas.calc_eas();
      
      //save data
      kappa_a << tempL << ' ' << eas.kappa_a_cgs(NUE) << " " << eas.kappa_a_cgs(NUA) << " " << eas.kappa_a_cgs(NUX) << "\n";
      kappa_s << tempL << ' ' << eas.kappa_s_cgs(NUE) << " " << eas.kappa_s_cgs(NUA) << " " << eas.kappa_s_cgs(NUX) << "\n";
      BB << eas.BB_cactus(NUE) << " " << eas.BB_cactus(NUA) << " " << eas.BB_cactus(NUX) << '\n';
      kNNB << tempL << ' ' << eas.kappa_a_brems_cgs(NUE) << ' ' << eas.kappa_a_brems_cgs(NUA) << ' ' << eas.kappa_a_brems_cgs(NUX) << '\n';
      kPP << tempL << ' ' << eas.kappa_a_epannihil_cgs(NUE) << ' ' << eas.kappa_a_epannihil_cgs(NUA) << ' ' << eas.kappa_a_epannihil_cgs(NUX) << '\n';
      kPD << tempL << ' ' << eas.kappa_a_plasmon_cgs(NUE) << ' ' << eas.kappa_a_plasmon_cgs(NUA) << ' ' << eas.kappa_a_plasmon_cgs(NUX) << '\n' ;
      kappa_abs << tempL << ' ' << eas.kappa_abs_cgs(NUE) << ' ' << eas.kappa_abs_cgs(NUA) << ' ' << eas.kappa_abs_cgs(NUX) << '\n';
    }
  }
    {
    std::ofstream kappa_a{"kappa_a_ye.dat"};
    std::ofstream kappa_s{"kappa_s_ye.dat"};
    std::ofstream kPP{"kappa_pairprocs_ye.dat"};
    std::ofstream kNNB{"kappa_brems_ye.dat"};
    std::ofstream kPD("kappa_plasmon_decay_ye.dat"); 
    std::ofstream kappa_abs{"kappa_charged_currents_ye.dat"};
    std::ofstream BB{"BB.dat"};
    
    
    
    double h{(yemax-yemin)/npoints};

    for(int i=0; i<npoints; i++) {
      
      double const rhoL  = rho;
      double const tempL = temp;//tempmin + i*h;
      double const yeL   = yemin + i*h;
      
      // calc eas with Margherita at yeL tempL rhoL 
      const auto tau_init= compute_analytic_opacity(rhoL);
      std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
      auto F= Fugacities<double>(rhoL,tempL,yeL, std::move(tau_n));
      auto eas = EAS<double>(F);
      eas.calc_eas();
      
      //save data
      kappa_a << yeL << ' ' << eas.kappa_a_cgs(NUE) << " " << eas.kappa_a_cgs(NUA) << " " << eas.kappa_a_cgs(NUX) << "\n";
      kappa_s << yeL << ' ' << eas.kappa_s_cgs(NUE) << " " << eas.kappa_s_cgs(NUA) << " " << eas.kappa_s_cgs(NUX) << "\n";
      BB << eas.BB_cactus(NUE) << " " << eas.BB_cactus(NUA) << " " << eas.BB_cactus(NUX) << '\n';
      kNNB << yeL << ' ' << eas.kappa_a_brems_cgs(NUE) << ' ' << eas.kappa_a_brems_cgs(NUA) << ' ' << eas.kappa_a_brems_cgs(NUX) << '\n';
      kPP << yeL << ' ' << eas.kappa_a_epannihil_cgs(NUE) << ' ' << eas.kappa_a_epannihil_cgs(NUA) << ' ' << eas.kappa_a_epannihil_cgs(NUX) << '\n';
      kPD << yeL << ' ' << eas.kappa_a_plasmon_cgs(NUE) << ' ' << eas.kappa_a_plasmon_cgs(NUA) << ' ' << eas.kappa_a_plasmon_cgs(NUX) << '\n' ;
      kappa_abs << yeL << ' ' << eas.kappa_abs_cgs(NUE) << ' ' << eas.kappa_abs_cgs(NUA) << ' ' << eas.kappa_abs_cgs(NUX) << '\n';
    }
  }

}


