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

enum M1_interactions { BETA=0,
		       PLASMON_DECAY,
		       BREMS,
		       PAIR,
		       NINTERACTIONS };

// Let's have some debug output
#define STANDALONE
#define DEBUG_BETAEQ
#define DEBUG_FOR_REAL
#define COLDTABLE_SETUP

using CCTK_REAL = double ;

#include "../FIL_M1_headers.h"
#include "../../3D_Table/readtable.cc"
#include "../M1.hh"
#include "../M1_find_betaeq.hh"
#include "../../Margherita_EOS.h"


int main() {

  //

  std::cout << std::setprecision(15);
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  
  // EOS
  std::string table_name =
    std::string("/zhome/academic/HLRS/xfp/xfpmusol/EOS_Tables/scollapse/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
//                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  //std::string("/Users/emost/tmp/EOS/new/eos_tables/DD2_margherita_09-May-2017.h5");
  constexpr bool use_energy_shift = true;

  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, false);

  std::cout << "Read Table" << std::endl;

  double rho = 0.00111869;//7.582629253507815e14*Margherita_constants::RHOGF;
  double yl = 0.1;//0.09848139;//0.234693877551020;
  double temp = 0.01;//28.66098856;

  EOS_Tabulated::error_type error ;
  double eps = EOS_Tabulated::eps__temp_rho_ye(temp, rho, yl, error);
/*  
  double rho = exp(-8.06717343536804);
  double temp = exp(-3.18252892776548);
  double ye =      0.140581339415548;
  */
  
  static constexpr const double NOT_USED = 0. ;
  std::array<double,REQ_FLUID_VARS> Ufluid = { rho, yl, temp} ;
  std::array<bool,M1_interactions::NINTERACTIONS> wi{true} ;
  std::array<double,MAXNUMVARS> Urad = 
      {NOT_USED, 8.0984561e-15, 0.0, 0.0, 0.0,
       NOT_USED, 9.13910255e-15, 0.0, 0.0, 0.0,
       NOT_USED, 8.60386684e-15, 0.0, 0.0, 0.0,
       //Jthick, ZVECX, ZVECY, ZVECY
       NOT_USED, 0.0, 0.0, 0.0};

  double ye_eq,T_eq;
  bool success;
  size_t niter ;
  
  auto t1 = high_resolution_clock::now();
  std::tie(ye_eq,T_eq,success,niter) = compute_T_ye_betaeq<EOS_Tabulated>(Urad,Ufluid,wi) ;
  auto t2 = high_resolution_clock::now();
  std::cout << "Here is what we started with (minim): " << Ufluid[YE] << "\t" << Ufluid[TEMP] << std::endl ;
  std::cout << "Here is where we ended up (minim): " << ye_eq << "\t" << T_eq << std::endl;

  duration<double, std::milli> ms_double_1 = t2 - t1;
  std::cout << "Function minimization: " << ms_double_1.count() << "ms" << std::endl ;
  

  return 0;
}
