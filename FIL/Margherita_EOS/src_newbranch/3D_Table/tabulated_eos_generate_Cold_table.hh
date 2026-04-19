//
//  Copyright (C) 2021, Harry Ho-Yin Ng
//

#include "tabulated.hh"
#include "tabulated_readtable_scollapse.hh"
#include "tabulated_readtable_compose.hh"
#include <fstream>
#include <iostream>
#include <iomanip>

#define H5_USE_16_API 1
#include <hdf5.h>

#ifndef EOS_TABULATED_GENERATE_COLD_TABLE_HH
#define EOS_TABULATED_GENERATE_COLD_TABLE_HH

void EOS_Tabulated::generate_cold_table(double logrho[], int &nrho){
  using namespace Margherita_constants;

  //code test Passed!  when T < T_thr, the muonic beta eqm condition has bad convergence for rho > 1e-3
  // T_thr = 0.7 MeV for SFHo
  // T_thr = 1.5 MeV for DD2, but DD2 tends to have lower density, those unphysical regions not relevant
  typename EOS_Tabulated::error_type error;
  //// passed very high density region, but low density region, there is no muon -->
  double temp_pt = 0.5;
  double rho_pt  = 1e-14;
  double yle_pt, ymu_pt, e_pt, yp_pt, h_pt, nb_pt;
  double press_pt, eps_pt;
  double baryonic_mass;

  //baryonic_mass = m_neutron_MeV * MeV_to_erg / c2_cgs;
  baryonic_mass = EOS_Tabulated::baryon_mass;
  std::cout <<"\n"<< baryonic_mass <<"_"<<"baryonic mass"<<"\n";

  std::ofstream myfile ("npe_DD2_T0-5mev.rns");
  myfile << nrho << "\n";

    for (int i = 0; i < nrho; ++i) {
      e_pt = 0.0;
      press_pt = 0.0;
      nb_pt = 0.0;
      h_pt = 0.0;
      eps_pt = 0.0;
      
      rho_pt = exp(logrho[i]);
         press_pt = EOS_Tabulated::press_eps_yle_ymu__beta_eq__rho_temp(eps_pt, yp_pt, ymu_pt, 
                                                                  rho_pt, temp_pt, error);

      // all should in CGS
      press_pt = press_pt / PRESSGF;
      e_pt = rho_pt * (1.0 + eps_pt)/RHOGF;
      nb_pt = rho_pt/RHOGF/ baryonic_mass;
      h_pt = c2_cgs * log( (e_pt+press_pt/c2_cgs)/( rho_pt * mnuc_cgs ) );
      // use baryonic mass for h_pt --> -ve
      //h_pt = c2_cgs * log( (e_pt+press_pt/c2_cgs)/(nb_pt * baryonic_mass) );

      myfile<<std::scientific<< std::setprecision (15)<< e_pt <<"  "<< press_pt<<"  "<< h_pt<<"  "<<nb_pt << "\n";

  
    }
  myfile.close();
  std::cout <<"Cold EOS Table for RNS is finished"<<"\n";
  MPI_Abort(MPI_COMM_WORLD, 911);
};

#endif
