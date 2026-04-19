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


std::string const Q_name = "Q_GR1D_mode.dat";
std::string const ka_name = "ka_GR1D_mode.dat";
std::string const ks_name = "ks_GR1D_mode.dat";





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

  std::vector<double> rho,temp,ye;
  size_t num_ye,num_temp,num_rho;

  std::tie(ye,num_ye) = read_data_from_file("/Users/carlomusolino/M1Code/FIL_M1/tests/ye.dat");
  std::tie(temp,num_temp) = read_data_from_file("/Users/carlomusolino/M1Code/FIL_M1/tests/temp.dat");
  std::tie(rho,num_rho) = read_data_from_file("/Users/carlomusolino/M1Code/FIL_M1/tests/rho.dat");

  size_t const npoints = 800 ;
  
  std::ofstream kappa_a{ka_name.c_str()};
  std::ofstream kappa_s{ks_name.c_str()};
  std::ofstream Q{Q_name.c_str()};

  std::ofstream QPP{"Q_pairprocs.dat"};
  std::ofstream QL{"Q_K.dat"};
  std::ofstream QNNB{"Q_brems.dat"};
  std::ofstream QB{"Q_beta.dat"};

  std::ofstream kappa_abs{"kappa_charged_currents.dat"};
  std::ofstream BB{"BB.dat"};

  std::ofstream eps2leak{"eps2_leak.dat"};
  std::ofstream eps2fluid{"eps2_fluid.dat"};
   
  for(int i=0; i<npoints; i++) {
	
    double const rhoL = rho[i];
    double const tempL = temp[i];
    double const yeL = ye[i];
    
    // calc eas with Margherita at yeL tempL rhoL 
    const auto tau_init= compute_analytic_opacity(rhoL);
    std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
    auto F= Fugacities<double>(rhoL,tempL,yeL, std::move(tau_n));
    auto eas = EAS<double>(F);
    eas.calc_eas();
    
    //save data
    Q << eas.Q_cactus(NUE) << " " << eas.Q_cactus(NUA) << " " << eas.Q_cactus(NUX) << "\n" ;
    kappa_a << eas.kappa_a_cactus(NUE) << " " << eas.kappa_a_cactus(NUA) << " " << eas.kappa_a_cactus(NUX) << "\n";
    kappa_s << eas.kappa_s_cactus(NUE) << " " << eas.kappa_s_cactus(NUA) << " " << eas.kappa_s_cactus(NUX) << "\n";
    BB << eas.BB_cactus(NUE) << " " << eas.BB_cactus(NUA) << " " << eas.BB_cactus(NUX) << '\n';
    QB << eas.Q_beta_cactus(NUE) << ' ' << eas.Q_beta_cactus(NUA) << ' ' << eas.Q_beta_cactus(NUX) << '\n';
    QNNB << eas.Q_brems_cactus(NUE) << ' ' << eas.Q_brems_cactus(NUA) << ' ' << eas.Q_brems_cactus(NUX) << '\n';
    QPP << eas.Q_pp_cactus(NUE) << ' ' << eas.Q_pp_cactus(NUA) << ' ' << eas.Q_pp_cactus(NUX) << '\n';
    QL << eas.Q_K_cactus(NUE) << ' ' << eas.Q_K_cactus(NUA) << ' ' << eas.Q_K_cactus(NUX) << '\n' ;
    kappa_abs << eas.kappa_abs_cactus(NUE) << ' ' << eas.kappa_abs_cactus(NUA) << ' ' << eas.kappa_abs_cactus(NUX) << '\n';

    eps2leak << eas.eps2_leak_cgs(NUE) << ' ' << eas.eps2_leak_cgs(NUA) << ' ' << eas.eps2_leak_cgs(NUX) << '\n' ;
    eps2fluid << eas.eps2_fluid_cgs(NUE) << ' ' << eas.eps2_fluid_cgs(NUA) << ' ' << eas.eps2_fluid_cgs(NUX) << '\n';
  }

}


