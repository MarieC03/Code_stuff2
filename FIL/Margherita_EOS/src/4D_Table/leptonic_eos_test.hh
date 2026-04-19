//
//  Copyright (C) 2021, Harry Ho-Yin Ng
//

#include "leptonic_eos.hh"
#include "leptonic_eos_readtable_leptonic.hh"
#include "leptonic_eos_readtable_scollapse.hh"
#include "leptonic_eos_readtable_compose.hh"
//#include "tabulated.hh"
#include <fstream>
#include <iostream>
#include <iomanip>

#define H5_USE_16_API 1
#include <hdf5.h>

#ifndef EOS_LEPTONIC_GENERATE_COLD_TABLE_HH
#define EOS_LEPTONIC_GENERATE_COLD_TABLE_HH

void EOS_Leptonic::test_leptonic_table(double logrho[], int &nrho){
  using namespace Margherita_constants;

  //code test Passed!  when T < T_thr, the muonic beta eqm condition has bad convergence for rho > 1e-3
  // T_thr = 0.7 MeV for SFHo
  // T_thr = 1.5 MeV for DD2, but DD2 tends to have lower density, those unphysical regions not relevant
  typename EOS_Leptonic::error_type error;
  //// passed very high density region, but low density region, there is no muon -->
  double temp_pt = 0.1;
  double rho_pt  = 1e-14;
  double yle_pt, ymu_pt, e_pt, yp_pt, h_pt, nb_pt;
  double press_pt, eps_pt;
  double baryonic_mass;


  std::ofstream myfile ("Kc2p_relative_diff_fixedYeYmu.log");
  //myfile << nrho << "\n";

    for (int i = 0; i < nrho; ++i) {
      rho_pt = exp(logrho[i]);
    }
  myfile.close();
  CCTK_ERROR("Testing for Kc2p with muons is finished");

};

#endif
