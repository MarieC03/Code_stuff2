//
//
// Quick Test for the TOV solver. Results should match values
// given in the documentation of WhiskyTOV!
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


// Let's have some debug output
#define STANDALONE
// Table debugging output
//#define DEBUG

#define COLDTABLE_SETUP
#define PWPOLY_SETUP

#include "../../margherita.hh"
#include "../../Margherita_EOS.h"
#include "../../Cold/cold_table_implementation.hh"
#include "../../Cold/setup_cold_table.cc"
#include "../tov.hh"

using this_EOS_t = Cold_EOS_Wrapper<Cold_Table>;

auto const num_models = 200;

int main() {

  using namespace Margherita_constants;

  std::cout << std::setprecision(15);

  std::string name {"./SLy4_betaeq_zeroT.lorene"};

  setup_Cold_Table(name,1000);

  auto tov_seq = std::array<std::unique_ptr<MargheritaTOV<this_EOS_t>>,num_models>();

  for(auto &tov : tov_seq){
	tov = std::make_unique<MargheritaTOV<this_EOS_t>>(1.e-3);
  }
  

// Models:
  auto const rho_min = -3.5;
  auto const rho_max = -1.5;


  auto const delta_rho = (rho_max-rho_min)/(num_models-1);



#pragma omp parallel for
  for(int nn = 0; nn<num_models; ++nn){

    auto rhoL = std::pow(10.,rho_min + nn*delta_rho);
    //Need pressure at center
    auto error = typename this_EOS_t::error_t {};
    double eps;
    auto const pressL = this_EOS_t::EOS_t::press_cold_eps_cold__rho(eps,rhoL,error);

    // EOS
    tov_seq[nn]->press_c= pressL;

    tov_seq[nn]->compact_output=true;
    tov_seq[nn]->adaptive_step = true;

    tov_seq[nn]->solve();

  }

  for(auto &tov : tov_seq){
    std::cout << *tov;
  }

  return 0;
}
