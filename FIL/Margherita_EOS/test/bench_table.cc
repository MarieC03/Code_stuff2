//
//
// Quick and dirty test to benchmark the root finding
//

//
//  Newton = false
//  16.1825802 secs on average for max_points = 1000000;

//  Newton = true
//  30.677783 secs on average for max_points = 1000000;

//

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Let's have some debug output
#define STANDALONE
// Table debugging output
//#define DEBUG

#include "../src/Margherita_EOS.h"
#include "../src/margherita.hh"
#include "../src/readtable.cc"

#include "../src/c2p/margherita_c2p.cc"
#include "../src/margherita_c2p.hh"
#include "../src/tabulated_implementation.hh"

// This is from
// https://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c
#include <sys/time.h>
typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
  struct timeval now;
  gettimeofday(&now, NULL);
  return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int main() {
  // how many points should be generated?
  constexpr size_t max_points = 200;
  constexpr size_t it_max = 1;
  //
  std::cout << std::setprecision(15);

  // EOS
  std::string table_name =
      std::string("/home/astro/most/"
                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  constexpr bool use_energy_shift = true;

  EOS_Tabulated::readtable(table_name.c_str(), use_energy_shift);

  std::cout << "Reading table..." << std::endl;

  EOS_Tabulated::error_type error; // TODO How do we handle this?

  // Create random input data

  // pseudo random number generator
  std::random_device rd;  // seed
  std::mt19937 gen(rd()); // generator

  // rho uniform distrib
  std::uniform_real_distribution<> dist_lrho(
      std::log(EOS_Tabulated::eos_rhomin), std::log(EOS_Tabulated::eos_rhomax));
  // ye uniform distrib
  std::uniform_real_distribution<> dist_ye(EOS_Tabulated::eos_yemin,
                                           EOS_Tabulated::eos_yemax);
  // temp uniform distrib
  std::uniform_real_distribution<> dist_ltemp(
      std::log(EOS_Tabulated::eos_tempmin),
      std::log(EOS_Tabulated::eos_tempmax));

  std::cout << "rho max: " << EOS_Tabulated::eos_rhomin
            << " , max: " << EOS_Tabulated::eos_rhomax << std::endl;
  std::cout << "temp max: " << EOS_Tabulated::eos_tempmin
            << " , max: " << EOS_Tabulated::eos_tempmax << std::endl;
  std::cout << "ye max: " << EOS_Tabulated::eos_yemin
            << " , max: " << EOS_Tabulated::eos_yemax << std::endl;

  // total seconds
  double tot_secs = 0;

  for (auto it = 0; it < it_max; ++it) {
    // vector of random data
    std::vector<std::array<double, 4>> data(max_points);

    std::cout << "Creating random data of " << max_points << " length..."
              << std::endl;

    for (auto i = 0; i < max_points; ++i) {
      data[i][0] = std::exp(dist_lrho(gen));
      data[i][1] = dist_ye(gen);
      data[i][2] = std::exp(dist_ltemp(gen));

      // get eps in valid range for rho, ye
      auto eps_range =
          EOS_Tabulated::eps_range__rho_ye(data[i][0], data[i][1], error);
      // apply shift to sample in log
      eps_range[0] = eps_range[0] + EOS_Tabulated::energy_shift;
      eps_range[1] = eps_range[1] + EOS_Tabulated::energy_shift;

      std::uniform_real_distribution<> dist_leps(std::log(eps_range[0]),
                                                 std::log(eps_range[1]));
      //      std::cout << "Got a shifted eps range of: " << eps_range[0] << "
      //      -> " << eps_range[1] << std::endl;
      //      std::cout << "Got an leps range of: " << std::log(eps_range[0]) <<
      //      " -> " << std::log(eps_range[1]) << std::endl;
      data[i][3] = std::exp(dist_leps(gen)) - EOS_Tabulated::energy_shift;
      double leps = dist_leps(gen);
      data[i][3] = std::exp(leps) - EOS_Tabulated::energy_shift;
      //      std::cout << "Picking eps to: " << data[i][3] << " (leps = " <<
      //      leps << ")" << std::endl;
      std::cout << "eps_min: " << eps_range[0] << std::endl;
      std::cout << "eps_max: " << eps_range[1] << std::endl;
      std::cout << "eps_sample: " << data[i][3] << std::endl;
      std::cout << "energy_sample: " << data[i][0] * (1. + data[i][3])
                << std::endl;
    }

    std::cout << "Root finding..." << std::endl;

    // Start timing
    timestamp_t t0 = get_timestamp();
    ////////    Invert conservatives

    // compute press and temp from random input
    std::vector<std::array<double, 2>> results(max_points);
    for (auto i = 0; i < max_points; ++i) {
      double rho_i = data[i][0];
      double ye_i = data[i][1];
      double temp_i = data[i][2];
      double eps_i = data[i][3];

      results[i][1] = data[i][2];

      std::cout << "rho[" << i << "]: " << data[i][0] << std::endl;
      std::cout << "ye[" << i << "]: " << data[i][1] << std::endl;
      std::cout << "temp_guess[" << i << "]: " << data[i][2] << std::endl;
      std::cout << "eps[" << i << "]: " << data[i][3] << std::endl;
      std::cout << std::endl;

      results[i][0] = EOS_Tabulated::press_temp__eps_rho_ye(
          results[i][1], data[i][3], data[i][0], data[i][1], error);
      results[i][3] = EOS_Tabulated::eps__temp_rho_ye(results[i][1], data[i][0],
                                                      data[i][1], error);
      std::cout << std::endl;
      std::cout << "press_result[" << i << "]: " << results[i][0] << std::endl;
      std::cout << "temp_result[" << i << "]: " << results[i][1] << std::endl;
      std::cout << "eps_result[" << i << "]: " << results[i][3] << std::endl;
      std::cout << std::endl;

      //      std::cout << "dRho/Rho:" << std::abs(rho_i - data[i][0])/rho_i <<
      //      std::endl;
      //      std::cout << "dYe/Ye:" << std::abs(ye_i - data[i][1])/ye_i <<
      //      std::endl;
      //      std::cout << "dT/T: " << std::abs(temp_i - results[i][1])/temp_i
      //      << std::endl;
      //      std::cout << "dEps/Eps:" << std::abs(eps_i - data[i][3])/eps_i <<
      //      std::endl;
      //      std::cout << "Press: " << results[i][0] << std::endl;
    }

    ////////
    // Stop timing
    timestamp_t t1 = get_timestamp();

    std::cout << "Finished..." << std::endl;

    double secs = (t1 - t0) / 1000000.0L;
    tot_secs += secs;
    std::cout << "Root finding took " << secs << " secs." << std::endl;
  }

  std::cout << "It took " << tot_secs / it_max << " seconds on average."
            << std::endl;
  /*

    // Add to global timer
    secs += (t1 - t0) / 1000000.0L;

    secs /= MAXIT;

  */
  return 0;
}
