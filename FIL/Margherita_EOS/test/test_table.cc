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
// Table debugging output
//#define DEBUG

#include "../src/margherita.hh"
#include "../src/3D_Table/readtable.cc"

#include "../src/margherita_c2p.hh"
#include "../src/c2p/register_c2p.C"

constexpr int MAXIT =
    100000; // Need to sample many inversions for performance measurement
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

  //

  std::cout << std::setprecision(15);

  // EOS
  std::string table_name =
//    std::string("/Users/emost/Downloads/"
//                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  std::string("/Users/emost/tmp/EOS/new/eos_tables/DD2_margherita_09-May-2017.h5");
  constexpr bool use_energy_shift = true;

  EOS_Tabulated::readtable(table_name.c_str(), use_energy_shift, false);

  std::cout << "Read Table" << std::endl;

  // Try flat space

  std::array<double, 6> metric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};
  std::array<double, 6> invmetric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};

  metric_c METRIC(metric, invmetric, 1.0);

  // Set primitives
  std::array<double, C2P_Hydro<EOS_Tabulated>::numprims> PRIMS{{
      1.234e-13, // RHOB
      0.0,      // PRESSURE
      0.178,    // ZVECX
      -0.18,    // ZVECY
      0.103,    // ZVECZ
      5.e-1,    // EPS
      0.0,      // CS2
      0.0,      // ENTROPY
      0.0,      // TEMP
      0.3       // Ye
  }};

  std::array<double, C2P_Hydro<EOS_Tabulated>::numprims> PRIMSN;

  EOS_Tabulated::error_type error; // TODO How do we handle this?
  double h;

  PRIMS[TEMP] = 0.0;

  // 4. Compute a by making a pressure call
  // Update all vars here, cs2, temp, etc..
  PRIMS[PRESSURE] = EOS_Tabulated::press_h_csnd2_temp_entropy__eps_rho_ye(
      h, PRIMS[CS2], PRIMS[TEMP], PRIMS[ENTROPY], // out all but temp (inout)
      PRIMS[EPS], PRIMS[RHOB], PRIMS[YE], error); // in

  // Declare conservatives
  std::array<double, C2P_Hydro<EOS_Tabulated>::numcons> CONS;

  double secs = 0.;

  // We need to register limits
  //
  std::array<double, Margherita_C2P_Registration::MAX_NUM_LIMITS> limitsL{{
      1.0e-14, // rho_atm,
      0.0,     // temp_atm,
      1.0,     // ye_atm,
      2.,      // eps_max,
      1.e-14,  // tauenergy_min,
      10.,     // gamma_max,
      100      /// psi6_threshold
  }};
  Margherita_register_limits<EOS_Tabulated>(limitsL);

  // Invert conservatives for MAXIT times to benchmark the inversion
  for (int i = 0; i < MAXIT; ++i) {

    // Compute conservatives
    CONS = C2P_Hydro<EOS_Tabulated>::compute_conservatives(PRIMS, METRIC);
    PRIMSN[YE] = PRIMS[YE];
    PRIMSN[TEMP] = PRIMS[TEMP];

    // Start timing
    timestamp_t t0 = get_timestamp();
    ////////    Invert conservatives

    error_bits_t error_bits;
    Margherita_Con2Prim<C2P_Hydro, EOS_Tabulated>(CONS, PRIMSN, METRIC,
                                                  error_bits);

    ////////
    // Stop timing
    timestamp_t t1 = get_timestamp();

    // Add to global timer
    secs += (t1 - t0) / 1000000.0L;
  }
  secs /= MAXIT;

  const std::array<std::string, NUM_PRIMS> names{
      {"RHOB", "PRESSURE", "ZVECX", "ZVECY", "ZVECZ", "EPS", "CS2", "ENTROPY",
       "TEMP", "YE"}};

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
