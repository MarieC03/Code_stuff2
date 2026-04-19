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
// Table debugging output
#define DEBUG

#include "../src/margherita_c2p.hh"
#include "../src/Margherita_EOS.h"
#include "../src/margherita.hh"
#include "../src/3D_Table/readtable.cc"

using namespace Margherita_C2P_conventions;
int main() {

  //

  std::cout << std::setprecision(15);

  // EOS
  std::string table_name =
      std::string("/Users/emost/Downloads/"
                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  constexpr bool use_energy_shift = true;

  EOS_Tabulated::readtable(table_name.c_str(), use_energy_shift, false);

  std::cout << "Read Table" << std::endl;

  // Try flat space

  std::array<double, 6> metric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};
  std::array<double, 6> invmetric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};

  metric_c METRIC(metric, invmetric, 1.0);

  // Set primitives
  std::array<double, NUM_PRIMS> PRIMS{{
      1.234e-14, // RHOB
      0.0,       // PRESSURE
      0.178,     // ZVECX
      -0.18,     // ZVECY
      0.103,     // ZVECZ
      5.e-1,     // EPS
      0.0,       // CS2
      0.0,       // ENTROPY
      0.01,      // TEMP
      0.3        // Ye
  }};

  std::array<double, NUM_PRIMS> PRIMSN;

  EOS_Tabulated::error_type error; // TODO How do we handle this?
  double h;

  // 4. Compute a by making a pressure call
  // Update all vars here, cs2, temp, etc..
  PRIMS[PRESSURE] = EOS_Tabulated::press_eps_ye__beta_eq__rho_temp(
      PRIMS[EPS], PRIMS[YE], PRIMS[RHOB], PRIMS[TEMP], error);

  const std::array<std::string, NUM_PRIMS> names{
      {"RHOB", "PRESSURE", "ZVECX", "ZVECY", "ZVECZ", "EPS", "CS2", "ENTROPY",
       "TEMP", "YE"}};

  std::cout << std::endl;
  // Output old and new primitives
  for (int i = 0; i < NUM_PRIMS; ++i) {
    std::cout << names[i] << ": " << PRIMS[i] << std::endl;
  }

  return 0;
}
