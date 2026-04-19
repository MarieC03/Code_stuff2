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
#define DEBUG

#define PWPOLY_SETUP

#include "../src/Cold/cold_pwpoly_implementation.hh"
#include "../src/Hybrid/hybrid_implementation.hh"
#include "../src/Margherita_EOS.h"

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
  constexpr int num_pieces = 1;
  constexpr double entropy_min = 1.0e-9;
  constexpr double gamma_th = 2.0;

  // First create Arrays..
  Cold_PWPoly::num_pieces = num_pieces;
  Cold_PWPoly::rhomax = 1.;
  Cold_PWPoly::rhomin = 1.e-20;

  EOS_Polytropic::gamma_th_m1 = (gamma_th - 1.);
  EOS_Polytropic::entropy_min = (entropy_min);

  // Set EOS
  // Ideal fluid without cold part
  Cold_PWPoly::k_tab[0] = 0.0;
  Cold_PWPoly::gamma_tab[0] = 2.0;
  Cold_PWPoly::rho_tab[0] = 0.0;
  Cold_PWPoly::eps_tab[0] = 0.0;
  Cold_PWPoly::P_tab[0] = 0.0;

  // Try flat space

  std::array<double, 6> metric{{1.e4, 0.0, -10.0, 50.0, 1., 50.0}};

  metric_c METRIC(metric);
  std::cout << "sqrtdet: " << METRIC.sqrtdet << std::endl;

  if (std::isnan(METRIC.sqrtdet)) {
    std::cout << "Metric is not positive definite" << std::endl;
    abort();
  }

  // Set primitives
  std::array<double, C2P_MHD<EOS_Polytropic>::numprims> PRIMS{{
      1.234e-12, // RHOB
      0.0,      // PRESSURE
      0.078,    // ZVECX
      -0.18,    // ZVECY
      0.003,    // ZVECZ
      1.e-1,    // EPS
      0.0,      // CS2
      0.0,      // ENTROPY
      0.0,      // TEMP
      0.0       // Ye
  }};

  std::array<double, 3> Bvec = {{
      1.e-6, 4.e-5, 1.e-7
  }};

  /*  std::array<double, 3> Bvec = {{
                                  0.,
                                  0.,
                                  0.//1.e-7
                                  }} ;

  */

  // We need to register limits
  //
  std::array<double, Margherita_C2P_Registration::MAX_NUM_LIMITS> limitsL{{
      1.0e-15, // rho_atm,
      0.0,     // temp_atm,
      0.0,     // ye_atm,
      1.e1,   // eps_max,
      0.,  // tauenergy_min,
      20.,     // gamma_max,
      1000000  /// psi6_threshold
  }};
  Margherita_register_limits<EOS_Polytropic>(limitsL);

  std::array<double, C2P_MHD<EOS_Polytropic>::numprims> PRIMSN;

  EOS_Polytropic::error_type error; // TODO How do we handle this?
  double h;
  // 4. Compute a by making a pressure call
  // Update all vars here, cs2, temp, etc..
  PRIMS[PRESSURE] = EOS_Polytropic::press_h_csnd2_temp_entropy__eps_rho_ye(
      h, PRIMS[CS2], PRIMS[TEMP], PRIMS[ENTROPY], // out all but temp (inout)
      PRIMS[EPS], PRIMS[RHOB], PRIMS[YE], error); // in

  // Declare conservatives
  std::array<double, C2P_MHD<EOS_Polytropic>::numcons> CONS;
  // Need to compute Lorentz factor and zvec_low
  std::array<double, 3> zvec_lowL = METRIC.lower_index<ZVECX, NUM_PRIMS>(PRIMS);
  const double z2L = PRIMS[ZVECX] * zvec_lowL[0] + PRIMS[ZVECY] * zvec_lowL[1] +
                     PRIMS[ZVECZ] * zvec_lowL[2];
  const double WL = sqrt(1. + z2L);

  double secs = 0.;

  double residual =-1;
  // Invert conservatives for MAXIT times to benchmark the inversion
  for (int i = 0; i < MAXIT; ++i) {

    // Compute conservatives
    CONS = C2P_MHD<EOS_Polytropic>::compute_conservatives(PRIMS, METRIC, Bvec);
    // Fun test with hydro conservatives
    /*  auto CONS2 = C2P_Hydro<EOS_Polytropic>::compute_conservatives(PRIMS,
      METRIC);
      for(int j=0; j<C2P_Hydro<EOS_Polytropic>::numcons; ++j) {
        std::cout << CONS[j] << "  ,  " << CONS2[j] <<std::endl;
        CONS[j]=CONS2[j];
      }
  */
    // Start timing
    timestamp_t t0 = get_timestamp();
    ////////    Invert conservatives

    error_bits_t error_bits;
    residual = Margherita_Con2Prim<C2P_MHD, EOS_Polytropic>(CONS, PRIMSN, METRIC,
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
  std::cout << "Specific enthalpy * lorentz: " << h * WL << std::endl;
  std::cout << "residual: " << residual << std::endl;
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
