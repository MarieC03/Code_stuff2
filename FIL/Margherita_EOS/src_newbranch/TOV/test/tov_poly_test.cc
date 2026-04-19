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
#include "../../Cold/cold_pwpoly_implementation.hh"
#include "../tov.hh"

using this_EOS_t = Cold_EOS_Wrapper<Cold_PWPoly>;


int main() {

  using namespace Margherita_constants;

  std::cout << std::setprecision(15);

  // EOS
  constexpr int num_pieces = 1;
  constexpr double entropy_min = 1.0e-9;
  constexpr double gamma_th = 2.0;

  // First create Arrays..
  Cold_PWPoly::num_pieces = num_pieces;
  Cold_PWPoly::rhomax = 1.;
  Cold_PWPoly::rhomin = 1.e-20;

  // Set EOS
  // Ideal fluid without cold part
  Cold_PWPoly::k_tab[0] = 100.0;
  Cold_PWPoly::gamma_tab[0] = 2.0;
  Cold_PWPoly::rho_tab[0] = 0.0;
  Cold_PWPoly::eps_tab[0] = 0.0;
  Cold_PWPoly::P_tab[0] = 0.0;

  double rho_c = 1.28e-3;
  auto error = typename this_EOS_t::error_t {};
  double eps;
  auto const pressL = this_EOS_t::EOS_t::press_cold_eps_cold__rho(eps,rho_c,error);
  // EOS
  MargheritaTOV<this_EOS_t> tov(pressL);

  tov.solve();

  std::cout << tov <<std::endl;
  std::cout << "Calculation took " << tov.get_iterations() << " Iterations." <<std::endl;


  std::cout << "Switching to adaptive radius stepping" << std::endl;

  tov.reset();
  tov.adaptive_step = true;

  tov.solve();

  std::cout << tov <<std::endl;
  std::cout << "Calculation took " << tov.get_iterations() << " Iterations." <<std::endl;

  return 0;
}
