#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

// Let's have some debug output
#define STANDALONE
#define COLDTABLE_SETUP
#define PWPOLY_SETUP
//#define DEBUG

typedef double CCTK_REAL;

#define C2P_TYPE_1 C2P_MHD
#define C2P_TYPE_2 C2P_MHD_KASTAUN

#include "../../src/margherita.hh"
#include "../../src/3D_Table/readtable.cc"

#include "../../src/margherita_c2p.hh"
#include "../../src/c2p/register_c2p.C"

#include <benchmark/benchmark.h>

using namespace Margherita_C2P_conventions ; 

static void BM_kc2p(benchmark::State& state) {
    std::string table_name =
//      std::string("/Users/emost/Downloads/"
//                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
    //std::string("/Users/emost/Downloads/FTYSS_EOS_rho220_temp180_ye65_version2.0_20171101.h5");
    std::string("/home/relastro-shared/EOS/scollapse/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5") ;
  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;

  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, false);

  std::cout << "Read Table" << std::endl;

  // Try flat space

  std::array<double, 6> metric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};
  std::array<double, 6> invmetric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};

  metric_c METRIC(metric, invmetric, 1.0);

  // Set primitives
  std::array<double, C2P_MHD<EOS_Tabulated>::numprims> PRIMS{{
      1.234e-4, // RHOB
      0.0,       // PRESSURE
      0.078,     // ZVECX
      -0.18,     // ZVECY
      0.003,     // ZVECZ
      1.e-0,      // EPS
      0.0,       // CS2
      0.0,       // ENTROPY
      100.0,       // TEMP
      0.1        // Ye
  }};

  std::array<double, 3> Bvec = {{
      1.e-4, 1.e-5,
      0.5e-7 // 1.e-7
  }};

  /*  std::array<double, 3> Bvec = {{
                                  0.,
                                  0.,
                                  0.//1.e-7
                                  }} ;

  */
  std::array<double, Margherita_C2P_Registration::MAX_NUM_LIMITS> limitsL{{
      1.0e-10, // rho_atm,
      0.0,     // temp_atm,
      -1.,     // ye_atm,
      1.e2,    // eps_max,
      1.e-14,  // tauenergy_min,
      10.,     // gamma_max,
      100      /// psi6_threshold
  }};
  EOS_Tabulated::atm_beta_eq = true;

  Margherita_register_limits<EOS_Tabulated>(limitsL);

  std::array<double, C2P_MHD<EOS_Tabulated>::numprims> PRIMSN;

  EOS_Tabulated::error_type error; // TODO How do we handle this?

  PRIMS[EPS] = EOS_Tabulated::eps__temp_rho_ye(PRIMS[TEMP],PRIMS[RHOB],PRIMS[YE],error);

  double h;
  auto epsrange =
      EOS_Tabulated::eps_range__rho_ye(PRIMS[RHOB], PRIMS[YE], error);
  PRIMS[EPS] = std::max(epsrange[0], std::min(epsrange[1], PRIMS[EPS]));
  // 4. Compute a by making a pressure call
  // Update all vars here, cs2, temp, etc..
  PRIMS[PRESSURE] = EOS_Tabulated::press_h_csnd2_temp_entropy__eps_rho_ye(
      h, PRIMS[CS2], PRIMS[TEMP], PRIMS[ENTROPY], // out all but temp (inout)
      PRIMS[EPS], PRIMS[RHOB], PRIMS[YE], error); // in

  // Declare conservatives
  std::array<double, C2P_MHD<EOS_Tabulated>::numcons> CONS;
  // Need to compute Lorentz factor and zvec_low
  std::array<double, 3> zvec_lowL = METRIC.lower_index<ZVECX, NUM_PRIMS>(PRIMS);
  const double z2L = PRIMS[ZVECX] * zvec_lowL[0] + PRIMS[ZVECY] * zvec_lowL[1] +
                     PRIMS[ZVECZ] * zvec_lowL[2];
  const double WL = sqrt(1. + z2L);


  // Compute conservatives
  
  for(auto _ : state) {
    
    CONS = C2P_MHD<EOS_Tabulated>::compute_conservatives(PRIMS, METRIC, Bvec);
  // Fun test with hydro conservatives
  /*  auto CONS2 = C2P_Hydro<EOS_Polytropic>::compute_conservatives(PRIMS,
      METRIC);
      for(int j=0; j<C2P_Hydro<EOS_Polytropic>::numcons; ++j) {
      std::cout << CONS[j] << "  ,  " << CONS2[j] <<std::endl;
      CONS[j]=CONS2[j];
      }
  */
    PRIMSN[YE] = PRIMS[YE];
    PRIMSN[TEMP] = PRIMS[TEMP];
    ////////    Invert conservatives
  error_bits_t error_bits;
    Margherita_Con2Prim<C2P_TYPE_2, EOS_Tabulated>(CONS, PRIMSN, METRIC,
						 error_bits);
  }

  benchmark::ClobberMemory() ;
}

static void BM_pc2p(benchmark::State& state) {
  std::string table_name =
//      std::string("/Users/emost/Downloads/"
//                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
    //std::string("/Users/emost/Downloads/FTYSS_EOS_rho220_temp180_ye65_version2.0_20171101.h5");
    std::string("/home/relastro-shared/EOS/scollapse/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5") ;
  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;

  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, false);

  std::cout << "Read Table" << std::endl;

  // Try flat space

  std::array<double, 6> metric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};
  std::array<double, 6> invmetric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};

  metric_c METRIC(metric, invmetric, 1.0);

  // Set primitives
  std::array<double, C2P_MHD<EOS_Tabulated>::numprims> PRIMS{{
      1.234e-4, // RHOB
      0.0,       // PRESSURE
      0.078,     // ZVECX
      -0.18,     // ZVECY
      0.003,     // ZVECZ
      1.e-0,      // EPS
      0.0,       // CS2
      0.0,       // ENTROPY
      100.0,       // TEMP
      0.1        // Ye
  }};

  std::array<double, 3> Bvec = {{
      1.e-4, 1.e-5,
      0.5e-7 // 1.e-7
  }};

  /*  std::array<double, 3> Bvec = {{
                                  0.,
                                  0.,
                                  0.//1.e-7
                                  }} ;

  */
  std::array<double, Margherita_C2P_Registration::MAX_NUM_LIMITS> limitsL{{
      1.0e-10, // rho_atm,
      0.0,     // temp_atm,
      -1.,     // ye_atm,
      1.e2,    // eps_max,
      1.e-14,  // tauenergy_min,
      10.,     // gamma_max,
      100      /// psi6_threshold
  }};
  EOS_Tabulated::atm_beta_eq = true;

  Margherita_register_limits<EOS_Tabulated>(limitsL);

  std::array<double, C2P_MHD<EOS_Tabulated>::numprims> PRIMSN;

  EOS_Tabulated::error_type error; // TODO How do we handle this?

  PRIMS[EPS] = EOS_Tabulated::eps__temp_rho_ye(PRIMS[TEMP],PRIMS[RHOB],PRIMS[YE],error);

  double h;
  auto epsrange =
      EOS_Tabulated::eps_range__rho_ye(PRIMS[RHOB], PRIMS[YE], error);
  PRIMS[EPS] = std::max(epsrange[0], std::min(epsrange[1], PRIMS[EPS]));
  // 4. Compute a by making a pressure call
  // Update all vars here, cs2, temp, etc..
  PRIMS[PRESSURE] = EOS_Tabulated::press_h_csnd2_temp_entropy__eps_rho_ye(
      h, PRIMS[CS2], PRIMS[TEMP], PRIMS[ENTROPY], // out all but temp (inout)
      PRIMS[EPS], PRIMS[RHOB], PRIMS[YE], error); // in

  // Declare conservatives
  std::array<double, C2P_MHD<EOS_Tabulated>::numcons> CONS;
  // Need to compute Lorentz factor and zvec_low
  std::array<double, 3> zvec_lowL = METRIC.lower_index<ZVECX, NUM_PRIMS>(PRIMS);
  const double z2L = PRIMS[ZVECX] * zvec_lowL[0] + PRIMS[ZVECY] * zvec_lowL[1] +
                     PRIMS[ZVECZ] * zvec_lowL[2];
  const double WL = sqrt(1. + z2L);



  for( auto _: state) {
      // Compute conservatives
    CONS = C2P_MHD<EOS_Tabulated>::compute_conservatives(PRIMS, METRIC, Bvec);
  // Fun test with hydro conservatives
  /*  auto CONS2 = C2P_Hydro<EOS_Polytropic>::compute_conservatives(PRIMS,
      METRIC);
      for(int j=0; j<C2P_Hydro<EOS_Polytropic>::numcons; ++j) {
      std::cout << CONS[j] << "  ,  " << CONS2[j] <<std::endl;
      CONS[j]=CONS2[j];
      }
  */
    PRIMSN[YE] = PRIMS[YE];
    PRIMSN[TEMP] = PRIMS[TEMP];
  ////////    Invert conservatives
    error_bits_t error_bits;
    Margherita_Con2Prim<C2P_TYPE_1, EOS_Tabulated>(CONS, PRIMSN, METRIC,
						 error_bits);
  }
}


BENCHMARK(BM_kc2p) ;
BENCHMARK(BM_pc2p) ;

BENCHMARK_MAIN() ; 
