//
//
// Carlo Musolino 
//
//

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono> 

// Let's have some debug output
#define STANDALONE
//#define DEBUG_BETAEQ
#define COLDTABLE_SETUP
using CCTK_REAL = double ;

#include "../FIL_M1_headers.h"
#include "../../3D_Table/readtable.cc"
/*
#include "../M1.hh"
#include "../M1_find_betaeq.hh"
#include "../../Margherita_EOS.h"
*/
#include "../M1.hh"
#include "../fugacities.hh"


int main() {
  using namespace std ;
  using FG = Fugacities<double> ;


  std::string table_name =
    std::string("/Users/carlomusolino/numrel/EOS_Tables/EOSs/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  //                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  //std::string("/Users/emost/tmp/EOS/new/eos_tables/DD2_margherita_09-May-2017.h5");
  constexpr bool use_energy_shift = true;

  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, false);

  std::cout << "Read Table" << std::endl;
  
  double rho = 1e-04  ;
  double yl = 0.09848139;//0.234693877551020;
  double temp = 28.66098856;
  
  double tau_n = compute_analytic_opacity(rho) ;



  std::array<double,3> tau = {tau_n,tau_n,tau_n};

  auto F = Fugacities<double>(rho,temp,yl,tau) ;
  
  
  cout << tau_n << endl ;
  cout << F.ftau[NUE] + F.ftau[NUA] + F.ftau[NUX] << endl ;
}
