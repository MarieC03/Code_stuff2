
//
//  Copyright (C) 2017, Elias Roland Most
//  			<emost@th.physik.uni-frankfurt.de>
//  			Ludwig Jens Papenfort
//                      <papenfort@th.physik.uni-frankfurt.de>
//
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tabulated.hh"
#include "tabulated_implementation.hh"

#include "../4D_Table/leptonic_eos.hh"
//#include "../4D_Table/leptonic_eos_readtable_leptonic.hh"
//#include "../4D_Table/leptonic_eos_readtable_scollapse.hh"
//#include "../4D_Table/leptonic_eos_readtable_compose.hh"
#include "../4D_Table/leptonic_eos_implementation.hh"

#include <iostream>

#ifndef STANDALONE

extern "C" void Margherita_readtable_schedule(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

//code test
      //EOS_Leptonic::readtable_compose(nuceos_table_name);

  // Adding leptonic table first, later may need to subract electron contributions in baryonic table

  // Whether adds electrons contributions to eos
  if (add_ele_contributions_to_table) {
     if (CCTK_EQUALS(eos_type, "tabulated")){
        std::cout << "\nAdding electron contribution to eos (3D)\n" << std::endl;
        EOS_Tabulated::add_ele_contribution = true;
        EOS_Tabulated::readtable_leptonic(leptoniceos_table_name);
     } else {
        std::cout << "\nAdding electron contribution to eos (4D)\n" << std::endl;
        EOS_Leptonic::add_ele_contribution = true;
     }
  }

  // Leptonic eos for 4D eos
  if (CCTK_EQUALS(eos_type, "tabulated_leptonic")){
    if (generate_cold_table) {EOS_Leptonic::generate_cold_table_local = true;}
    if (use_muonic_table) {
      EOS_Leptonic::use_muonic_eos = true;
      std::cout << "\nUse Leptonic table (muonic table)\n" << std::endl;
      EOS_Leptonic::readtable_leptonic(leptoniceos_table_name);
      std::cout << "\nFinished reading leptonic table\n" << std::endl;
    } else {
        CCTK_ERROR("\n Why you need tabulated_leptonic without turning on use_muonic_table??\n") ;
         //EOS_Leptonic::use_muonic_eos = false;
         //std::cout << "\n Using EOS_Leptonic but not reading leptonic table\n" << std::endl;
     }
  } 

  // If your EOS contains electrons contribution, but you want to use leptonic eos contribution instead, 
  // you need to subtract baryon first to avoid double counting electrons
  if (baryon_eos_contains_e_contri && add_ele_contributions_to_table) { 
      if (CCTK_EQUALS(eos_type, "tabulated_leptonic")){
       	  EOS_Leptonic::presubtract_e_contri = true;
      } else {
       	  EOS_Tabulated::presubtract_e_contri = true;
      }
  }

  if(CCTK_EQUALS(nuceos_table_type, "stellarcollapse"))
    if( CCTK_EQUALS(eos_type, "tabulated_leptonic")){
      EOS_Leptonic::readtable_scollapse(nuceos_table_name,do_energy_shift,recompute_mu_nu) ; 
    } else{
      if (generate_cold_table) {EOS_Tabulated::generate_cold_table_local = true;}
      EOS_Tabulated::readtable_scollapse(nuceos_table_name, do_energy_shift, recompute_mu_nu);
    }
  if(CCTK_EQUALS(nuceos_table_type, "compose"))
    if( CCTK_EQUALS(eos_type, "tabulated_leptonic")){
      EOS_Leptonic::readtable_compose(nuceos_table_name);
    } else {
      if (generate_cold_table) {EOS_Tabulated::generate_cold_table_local = true;}
      EOS_Tabulated::readtable_compose(nuceos_table_name);
    }
  if(extend_table){
    EOS_Tabulated::extend_table_high = true;
    EOS_Leptonic::extend_table_high = true;
  }


  if( CCTK_EQUALS(eos_type, "tabulated_leptonic")){
     // Hard code, FIXME: put it to para file?
     // 0.5 or 1 MeV to be chosen since leptonic table beta eqm can only converge with this value of Temp
     EOS_Leptonic::temp_ID = std::max(20., EOS_Leptonic::eos_tempmin);
     //EOS_Leptonic::temp_ID = std::max(0.5, EOS_Leptonic::eos_tempmin);
     //EOS_Leptonic::temp_ID = std::max(0.1, EOS_Leptonic::eos_tempmin);
     std::cout << "Got rhomax,rhomin " << EOS_Leptonic::eos_rhomax << ","
               << EOS_Leptonic::eos_rhomin << std::endl;
     std::cout << "EOS_Leptonic::temp_ID " << EOS_Leptonic::temp_ID << "\n";
  } else {
     // Hard code, FIXME: put it to para file?
     EOS_Tabulated::temp_ID = std::max(20., EOS_Tabulated::eos_tempmin);
     //EOS_Tabulated::temp_ID = std::max(0.5, EOS_Tabulated::eos_tempmin);
     //EOS_Tabulated::temp_ID = std::max(0.1, EOS_Tabulated::eos_tempmin);
     EOS_Tabulated::eos_ylemin = EOS_Tabulated::eos_yemin;
     EOS_Tabulated::eos_ylemax = EOS_Tabulated::eos_yemax;
     std::cout << "Got rhomax,rhomin " << EOS_Tabulated::eos_rhomax << ","
               << EOS_Tabulated::eos_rhomin << std::endl;
     std::cout << "EOS_Tabulated::temp_ID " << EOS_Tabulated::temp_ID << "\n";
  }

  if (use_5_spec_m1 && !use_muonic_table) {
     CCTK_ERROR("When you use_5_spec_m1, you need to turn on use_muonic_table!");
  }
  if (use_muonic_table && slice_is_isentropic) {
     CCTK_ERROR("When you use_muonic_table, you cant use slice_is_isentropic, pls ask Harry to fix");
  }

}
#endif
