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

#include "../src/3D_Table/tabulated.hh"
#include "../src/3D_Table/tabulated_implementation.hh"
#include "../src/margherita.hh"

// This is from
// https://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c
#include <sys/time.h>
typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
  struct timeval now;
  gettimeofday(&now, NULL);
  return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

constexpr double SQ(const double &x) { return x * x; }

int main() {
  // how many points should be generated?
  constexpr size_t max_points = 10;
  constexpr size_t it_max = 1;
  //
  std::cout << std::setprecision(15);

  // EOS
  std::string table_name =
      std::string("/Users/emost/Downloads/"
                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
  constexpr bool use_energy_shift = true;

  EOS_Tabulated::readtable(table_name.c_str(), use_energy_shift, false);

  std::cout << "Reading table..." << std::endl;

  EOS_Tabulated::error_type error; // TODO How do we handle this?

  // Create random input data

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

  // total seconds
  double tot_secs = 0;

  for (auto it = 0; it < it_max; ++it) {
    // vector of random data
    std::vector<std::array<double, 5>> data(max_points);

    std::cout << "Creating random data of " << max_points << " length..."
              << std::endl;

    for (auto i = 0; i < max_points; ++i) {
      data[i][0] = std::exp(dist_rho(gen));
      data[i][1] = dist_ye(gen);
      data[i][2] = std::exp(dist_temp(gen));

      // get eps in valid range for rho, ye
      // auto eps_range = EOS_Tabulated::eps_range__rho_ye(data[i][0],
      // data[i][1],error);
      data[i][3] = EOS_Tabulated::eps__temp_rho_ye(data[i][2], data[i][0],
                                                   data[i][1], error);
      data[i][4] = EOS_Tabulated::press__temp_rho_ye(data[i][2], data[i][0],
                                                     data[i][1], error);
    }

    std::cout << "Root finding..." << std::endl;

    // Start timing
    timestamp_t t0 = get_timestamp();
    ////////    Invert conservatives
    //
    double eps_LS = 0.;

    // compute press and temp from random input
    std::vector<std::array<double, 2>> results(max_points);
    for (auto i = 0; i < max_points; ++i) {
      results[i][1] = data[i][3];
      std::uniform_real_distribution<> dist_temp2(0.8 * data[i][2],
                                                  1.2 * data[i][2]);
      double temp = dist_temp2(gen);

      results[i][0] = EOS_Tabulated::eps__press_temp_rho_ye(
          data[i][4], temp, data[i][0], data[i][1], error);
      eps_LS += SQ(results[i][1] - results[i][0]);
      std::cout << results[i][0] << " , " << results[i][1] << std::endl;
    }

    ////////
    // Stop timing
    timestamp_t t1 = get_timestamp();

    std::cout << "Finished..." << std::endl;

    double secs = (t1 - t0) / 1000000.0L;
    tot_secs += secs;
    std::cout << "Root finding took " << secs << " secs." << std::endl;
    std::cout << "Eps error is " << sqrt(eps_LS / max_points) << std::endl;
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
