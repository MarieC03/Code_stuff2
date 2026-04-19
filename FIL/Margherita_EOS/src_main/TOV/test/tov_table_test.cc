//
//
// Quick and dirty test to benchmark the inversion process
// and check it for consistency. This file is anything but clean!
//
//  Copyright (C) 2017, Elias Roland Most
//                      <emost@th.physik.uni-frankfurt.de>
//
//
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

//User parameters
constexpr double temp = 0.01;
constexpr int num_points = 10000;
constexpr int downsample_rho = 100;
constexpr bool h_larger_one  = true;
std::string table_name("/Users/carlomusolino/Desktop/tntyst/tntyst_new.h5");


// Let's have some debug output
#define STANDALONE
// Table debugging output
//#define DEBUG

#define COLDTABLE_SETUP

#include "../../margherita.hh"
#include "../../3D_Table/readtable.cc"
#include "../../Margherita_EOS.h"
#include "../../Cold/hot_slice.hh"
#include "../tov.hh"


int main() {

  using namespace Margherita_constants;

  // EOS
  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;

  EOS_Tabulated::readtable_compose(table_name.c_str());
  //Compute Hot slice
  HotSlice_beta_eq(downsample_rho, num_points,temp);

  MargheritaTOV<Hot_Slice> tov(1.e-3);

  tov.solve();

  std::cout << tov <<std::endl;

  return 0;
}
