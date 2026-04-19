////
//// This file is part of Margherita, the light-weight EOS framework
////
////  Copyright (C) 2017, Elias Roland Most
////                      <emost@th.physik.uni-frankfurt.de>
////
////  This program is free software: you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation, either version 3 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
////  You should have received a copy of the GNU General Public License
////  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _HH_MARGHERITA_C2P
#define _HH_MARGHERITA_C2P

#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include "Margherita_EOS.h"
#include "margherita.hh"
#include "utils/metric.hh"

namespace Margherita_C2P_Registration {
// Atmosphere
static constexpr int RHOB_ATM = 0;
static constexpr int TEMP_ATM = 1;
static constexpr int YE_ATM = 2;
// Con2Prim
static constexpr int EPS_MAX = 3;
static constexpr int TAUENERGY_MIN = 4;
static constexpr int LORENTZ_MAX = 5;
static constexpr int PSI6_BH = 6;

// So in total ALWAYS SET THIS CORRECTLY!
static constexpr int MAX_NUM_LIMITS = 7;
}  // namespace Margherita_C2P_Registration

// Some conventions
namespace Margherita_C2P_conventions {

static constexpr int WVT = 3;
static constexpr int WVX = 0;
static constexpr int WVY = 1;
static constexpr int WVZ = 2;

static constexpr int GXX = 0;
static constexpr int GXY = 1;
static constexpr int GXZ = 2;
static constexpr int GYY = 3;
static constexpr int GYZ = 4;
static constexpr int GZZ = 5;

// IMPORTANT: Convention for the conservatives and primitives:
// Conservatives
static constexpr int RHOSTAR = 0, STILDEX = 1, STILDEY = 2, STILDEZ = 3,
                     TAUENERGY = 4, YESTAR = 5, BX_CENTER = 6, BY_CENTER = 7,
                     BZ_CENTER = 8, ENTROPYSTAR = 9, EX_CENTER = 10,
                     EY_CENTER = 11, EZ_CENTER = 12;
// Primitives
static constexpr int RHOB = 0, PRESSURE = 1, ZVECX = 2, ZVECY = 3, ZVECZ = 4,
                     EPS = 5, CS2 = 6, ENTROPY = 7, TEMP = 8, YE = 9;

static constexpr int NUM_PRIMS = 10;

// This needs to go somewhere else
constexpr double SQ(const double &x) { return x * x; };

enum c2p_errors {
  RHOSTAR_ADJUSTED = 0,
  EPS_ADJUSTED,
  V_MAX_EXCEEDED,
  ENTROPY_FIX,  // This is only used by the entropy c2p
  C2P_FAILED,
  TAU_FIX,
  STILDE_FIX,
  STILDE_FIX_AH,
  YE_ADJUSTED,
  NUM_ERRORS  // This will automatically count the number of errors
};

typedef std::bitset<c2p_errors::NUM_ERRORS> error_bits_t;

};  // namespace Margherita_C2P_conventions

class Margherita_C2P_Limits {
 public:
  static double tau_min;
  static double lorentz_max;
  static double z_max;
  static double eps_max;
  static double psi6_bh;
};

// Include physics modules

#include "c2p/margherita_c2p_hydro.hh"
#include "c2p/margherita_c2p_mhd.hh"

#include "c2p/margherita_c2p.cc"

#endif
