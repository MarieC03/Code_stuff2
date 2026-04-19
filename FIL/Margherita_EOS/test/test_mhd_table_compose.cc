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
#define PWPOLY_SETUP
#define DEBUG

#define C2P_TYPE C2P_MHD_NEWMAN

#include "../src/margherita.hh"
#include "../src/3D_Table/readtable.cc"

#include "../src/margherita_c2p.hh"
#include "../src/c2p/register_c2p.C"

constexpr int MAXIT =
    1; // Need to sample many inversions for performance measurement
// Don't forget to disable STANDALONE when benchmarking over many iterations..

// This is from
// https://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c
#include <sys/time.h>
typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
  struct timeval now;
  gettimeofday(&now, NULL);
  return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

using namespace Margherita_C2P_conventions;
int main() {

  std::cout << std::setprecision(15);

  // EOS
  std::string table_name =
      std::string("/Users/emost/Downloads/eoscompose.h5");
//                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
//  std::string("/Users/emost/tmp/EOS/new/eos_tables/DD2_margherita_09-May-2017.h5");
  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;

  EOS_Tabulated::readtable_compose(table_name.c_str());

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

  double secs = 0.;

  // Start timing
  timestamp_t t0 = get_timestamp();
  // Invert conservatives for MAXIT times to benchmark the inversion
  for (int i = 0; i < MAXIT; ++i) {

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
    Margherita_Con2Prim<C2P_TYPE, EOS_Tabulated>(CONS, PRIMSN, METRIC,
                                                 error_bits);
  }
  ////////
  // Stop timing
  timestamp_t t1 = get_timestamp();
  // Add to global timer
  secs += (t1 - t0) / 1000000.0L;
  secs /= MAXIT;

  const std::array<std::string, NUM_PRIMS> names{
      {"RHOB", "PRESSURE", "ZVECX", "ZVECY", "ZVECZ", "EPS", "CS2", "ENTROPY",
       "TEMP", "YE"}};

  std::cout << std::endl;
  std::cout << "Specific enthalpy * lorentz: " << h * WL << std::endl;
  std::cout << std::endl;
  // Output old and new primitives
  for (int i = 0; i < NUM_PRIMS; ++i) {
    std::cout << names[i] << ": " << PRIMS[i] << " , " << PRIMSN[i]
              << std::endl;
  }

  // Output average timing
  std::cout << "\n Average time[s]: " << secs << std::endl;

  return 0;
}
