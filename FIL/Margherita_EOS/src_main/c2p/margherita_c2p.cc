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

#include "../margherita_c2p.hh"
#include "../Margherita_EOS.h"

#ifndef _CC_MARGHERITA_C2P
#define _CC_MARGHERITA_C2P
template <template <typename> class physics_t, typename eos>
static double Margherita_Con2Prim(
    std::array<double, physics_t<eos>::numcons> &CONS,
    std::array<double, physics_t<eos>::numprims> &PRIMS, const metric_c &METRIC,
    std::bitset<Margherita_C2P_conventions::c2p_errors::NUM_ERRORS>
        &error_bits) {
  using namespace Margherita_C2P_conventions;

  // Step 0. Make sure that you code has already caught CONS[RHOSTAR] <0  cases
  // BEFORE calling this routine!

  // Step 1. Rescale everything according to the conserved density
  // rho_star(=dens in Whisky based codes)

  CONS[STILDEX] /= CONS[RHOSTAR];
  CONS[STILDEY] /= CONS[RHOSTAR];
  CONS[STILDEZ] /= CONS[RHOSTAR];
  CONS[TAUENERGY] /= CONS[RHOSTAR];
  // Set YE
  PRIMS[YE] = CONS[YESTAR] / CONS[RHOSTAR];
  // Set entropy if present
  if (ENTROPYSTAR == CONS.size() - 1)
    PRIMS[ENTROPY] = CONS[ENTROPYSTAR] / CONS[RHOSTAR];

  // We need to check whether YE is within the table bounds
  if (PRIMS[YE] < eos::eos_yemin) {
    error_bits[c2p_errors::YE_ADJUSTED] = true;
    PRIMS[YE] = eos::eos_yemin;
  } else {
    if (PRIMS[YE] > eos::eos_yemax) {
      error_bits[c2p_errors::YE_ADJUSTED] = true;
      PRIMS[YE] = eos::eos_yemax;
    }
  }


  
  // Remove the sqrt(det) from rho_star
  CONS[RHOSTAR] /= METRIC.sqrtdet;

  // Step 1. Compute STILDE_SQ
  auto stilde_sq =
      METRIC.compute_inv_square_3_gen<STILDEX, physics_t<eos>::numcons>(CONS);
  
  // Enforce limits on conservatives
  // If you would like to know how this works, read the appendix of Etienne et
  // al. 2012
  // https://arxiv.org/pdf/1112.0568.pdf
  //
  if ( physics_t<eos>::fix_conservatives ) {
    // This is Etienne et al. 2012 (A.46)
    const auto stilde_sq_max = physics_t<eos>::limit_tau_and_return_stilde_sq_max(CONS, PRIMS, METRIC, error_bits);
        
#ifdef DEBUG
    std::cout << "TAUENERGY: " << CONS[TAUENERGY] << std::endl;
    std::cout << "Stilde: " << sqrt(stilde_sq)
	      << ",  MAX: " << sqrt(stilde_sq_max) << std::endl;
#endif
    
    // Limit stilde if necessary
    // See Etienne et al. 2012 (A.51)
    if (stilde_sq > stilde_sq_max && !(ENTROPYSTAR == CONS.size() - 1)) {
      const auto convfac = sqrt(stilde_sq_max / stilde_sq);
      
      CONS[STILDEX] *= convfac;
      CONS[STILDEY] *= convfac;
      CONS[STILDEZ] *= convfac;
      
      stilde_sq = stilde_sq_max;
      
      error_bits[c2p_errors::STILDE_FIX] = true;
    }

  } // if conservative fix  

  // Note: Additional limiting e.g. inside black holes highly depends on the
  // underlying physics,
  // e.g. hydro or mhd and is therefore done in the constructor of the physics
  // module

  // Instantiate physics module
  physics_t<eos> physics(&CONS, &PRIMS, &METRIC, &error_bits);

  // Set stilde and compute k = stilde/(tau +1)
  physics.set_stilde(stilde_sq);

  double residual;
  auto z_final = physics.invert(residual);

  // Update primitive variables
  physics.update_primitives(z_final);

  // Return residual error
  return residual;
};

#endif
