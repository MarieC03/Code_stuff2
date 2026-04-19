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
#include <random>

// Let's have some debug output
#define STANDALONE
//#define DEBUG
//#define OUTPUT
//#define TAUNEGATIVE

#include "../src/3D_Table/readtable.cc"

#include "../src/margherita_c2p.hh"
#include "../src/c2p/register_c2p.C"

constexpr int num_rho = 300;
constexpr int num_temp = 100;

constexpr double Ye = 0.2;

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

template <typename T> static inline double L1(T &A, T &B) {
  double sum = 0;
  for (int i = 0; i < A.size(); ++i) {
    sum += std::fabs(A[i] - B[i]) / std::max(std::fabs(A[i]), std::fabs(B[i]));
  }
  return sum / A.size();
}

using namespace Margherita_C2P_conventions;
int main() {

  std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(16);

  // EOS
  std::string table_name =
      std::string("/home/astro/most/"
                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;

  EOS_Tabulated::readtable(table_name.c_str(), use_energy_shift, false);


  // pseudo random number generator
  std::random_device rd;  // seed
  std::mt19937 gen(rd()); // generator

  // rho uniform distrib
  std::uniform_real_distribution<> dist_rho(
      std::log(EOS_Tabulated::eos_rhomin), std::log(EOS_Tabulated::eos_rhomax));
  // ye uniform distrib
  std::uniform_real_distribution<> dist_ye(EOS_Tabulated::eos_yemin,
                                           EOS_Tabulated::eos_yemax);
  // temp uniform distrib
  std::uniform_real_distribution<> dist_temp(
      std::log(EOS_Tabulated::eos_tempmin),
      std::log(EOS_Tabulated::eos_tempmax));

  std::uniform_real_distribution<> dist_v(-1., 1.);

  // Try flat space

  std::array<double, 6> metric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};
  std::array<double, 6> invmetric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};

  metric_c METRIC(metric, invmetric, 1.0);


  /*  std::array<double, 3> Bvec = {{
                                  0.,
                                  0.,
                                  0.//1.e-7
                                  }} ;

  */
  std::array<double, Margherita_C2P_Registration::MAX_NUM_LIMITS> limitsL{{
      1.0e-15, // rho_atm,
      0.0,     // temp_atm,
      -1.,     // ye_atm,
      1.e99,   // eps_max,
      1.e-14,  // tauenergy_min,
      10.,     // gamma_max,
      100      /// psi6_threshold
  }};
  EOS_Tabulated::atm_beta_eq = true;

  Margherita_register_limits<EOS_Tabulated>(limitsL);

  for (int nn_rho = 0; nn_rho < num_rho; nn_rho++)
    for (int nn_temp = 0; nn_temp < num_temp; nn_temp++) {

      std::array<double, C2P_Hydro<EOS_Tabulated>::numprims> PRIMSN;
      // Set primitives
      

       std::array<double, 3> Bvec = {};

      std::array<double, C2P_Hydro<EOS_Tabulated>::numprims> PRIMS{{
          1.234e-14, // RHOB
          0.0,       // PRESSURE
          0.078,     // ZVECX
          -0.08,     // ZVECY
          0.03,      // ZVECZ
          1.e-5,     // EPS
          0.0,       // CS2
          0.0,       // ENTROPY
          0.0,       // TEMP
          0.2        // Ye
      }};

      // We would like to iterate over rho, temp
      PRIMS[RHOB] = exp(dist_rho(gen));
      PRIMS[TEMP] = exp(dist_temp(gen));
      PRIMS[ZVECX] = dist_v(gen);
      PRIMS[ZVECY] = dist_v(gen);
      PRIMS[ZVECZ] = dist_v(gen);
      

      EOS_Tabulated::error_type error; // TODO How do we handle this?

      PRIMS[PRESSURE] = EOS_Tabulated::press__temp_rho_ye(
          PRIMS[TEMP], PRIMS[RHOB], PRIMS[YE], error);
      PRIMS[EPS] = EOS_Tabulated::eps_csnd2_entropy__temp_rho_ye(
          PRIMS[CS2], PRIMS[ENTROPY], PRIMS[TEMP], PRIMS[RHOB], PRIMS[YE],
          error);

      // Compute conservatives
      auto CONS =
          C2P_Hydro<EOS_Tabulated>::compute_conservatives(PRIMS, METRIC );
      PRIMSN[YE] = PRIMS[YE];
      PRIMSN[TEMP] = 0.0;
      ////////    Invert conservatives

#ifdef TAUNEGATIVE
     CONS[TAUENERGY] = -1.e10;
#endif

      typename Margherita_C2P_conventions::error_bits_t error_bits;
      auto residual =Margherita_Con2Prim<C2P_Hydro, EOS_Tabulated>(CONS, PRIMSN, METRIC,
                                                  error_bits);

      auto CONSN =
          C2P_Hydro<EOS_Tabulated>::compute_conservatives(PRIMSN, METRIC );
      CONS = C2P_Hydro<EOS_Tabulated>::compute_conservatives(PRIMS, METRIC);

#ifdef TAUNEGATIVE
      PRIMSN[EPS] = PRIMS[EPS];
      PRIMSN[CS2] = PRIMS[CS2];
      PRIMSN[ENTROPY] = PRIMS[ENTROPY];
      PRIMSN[TEMP] = PRIMS[TEMP];
      PRIMSN[PRESSURE] = PRIMS[PRESSURE];
#endif
      auto errorL1 = L1<>(PRIMS, PRIMSN);

#ifdef OUTPUT
      const std::array<std::string, NUM_PRIMS> names{
          {"RHOB", "PRESSURE", "ZVECX", "ZVECY", "ZVECZ", "EPS", "CS2",
           "ENTROPY", "TEMP", "YE"}};

      std::cout << std::endl;
      // Output old and new primitives
      for (int i = 0; i < NUM_PRIMS; ++i) {
        std::cout << names[i] << ": " << PRIMS[i] << " , " << PRIMSN[i]
                  << std::endl;
      }

#endif

      std::cout << PRIMS[RHOB] << "\t" << PRIMS[TEMP] << "\t" << errorL1 
                << "\t" << residual << std::endl;
    }
  return 0;
}
