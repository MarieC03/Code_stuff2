
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

#include "Weakhub.hh"
#include "Weakhub_implementation.hh"

#include <iostream>

#ifndef STANDALONE

extern "C" void Margherita_setup_Weakhub(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( CCTK_Equals( m1_which_eas, "Weakhub" ) ) {
     M1_Weakhub::M1_Weakhub::use_Weakhub_eas = true;

     std::cout << "\n Use Weakhub neutrino microphysics\n" << std::endl;
     M1_Weakhub::readtable_Weakhub(Weakhub_table_name);
     std::cout << "\n Finished reading Weakhub table\n" << std::endl;
  }

}
#endif
