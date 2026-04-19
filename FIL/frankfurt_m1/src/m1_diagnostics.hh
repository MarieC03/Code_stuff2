//
//  Copyright (C) 2021, Harry Ho-Yin Ng
//  Based on routines by Elias Roland Most
//  			 Ludwig Jens Papenfort
//
// TODO: Check table bounds!
//        When I need to extend the table like EOS_tabulated::

//#include "cctk.h"
#include <array>
#include <bitset>
#include <iostream>

#ifndef M1_DIAGNOSTICS_HH
#define M1_DIAGNOSTICS_HH

class M1_Diagnostics {
 public:
 //

 private:

 public:
  static double munu__rho_temp_ynu(const double &rho, const double &T, const double &ynu);

};


#endif
