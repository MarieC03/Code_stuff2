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

template <typename eos, template <typename> class formulation>
class C2P_MHD_t {
  friend eos;

 public:
  static constexpr bool fix_conservatives = formulation<eos>::fix_conservatives ;
  static constexpr int numcons = formulation<eos>::numcons;
  //static constexpr int numprims = 10;
  static constexpr int numprims = 11;

  typedef const metric_c metric;
  typedef std::array<double, numcons> cons;
  typedef std::array<double, numprims> prims;

  cons *__restrict CONS;
  prims *__restrict PRIMS;
  const metric *__restrict METRIC;
  Margherita_C2P_conventions::error_bits_t *__restrict error_bits;

  static constexpr double tol = formulation<eos>::tol;

  double stilde_sq;
  double B2tilde;
  double SdotBtilde;
  double gamma_lorentz;

  bool inside_BH;

  double aux_vars[formulation<eos>::num_aux];

 private:
  typedef Margherita_C2P_Limits limits;

  inline void enforce_MHD_BH_limits() {
    using namespace Margherita_C2P_conventions;
    if (inside_BH) {
      // Etienne 2012 et al. (A52)
      const double Wm = sqrt(SQ(SdotBtilde) + eos::c2p_h_min);
      // Etienne 2012 et al. (A53)
      const double Sm2 =
          (SQ(Wm) * stilde_sq + SQ(SdotBtilde) * (B2tilde + 2.0 * Wm)) /
          SQ(Wm + B2tilde);
      // Etienne 2012 et al. (A54)
      const double Wmin = sqrt(Sm2 + eos::c2p_h_min);

      // Etienne 2012 et al. (A40)
      const double tau_fluid_min =
          (*CONS)[TAUENERGY] - 0.5 * B2tilde -
          (B2tilde * stilde_sq - SQ(SdotBtilde)) * 0.5 / (SQ(Wmin + B2tilde));

      typename eos::error_type error;  // TODO How do we handle this?
      constexpr double safetyfac_lorentz = 1.2;
      double rhoL =
          (*CONS)[RHOSTAR] / (safetyfac_lorentz * limits::lorentz_max);

         const auto epsrange = eos::eps_range__rho_yle_ymu(rhoL, (*PRIMS)[YE], (*PRIMS)[YMU], error);

         // Etienne 2012 et al. (A55)
         if (tau_fluid_min < epsrange[0]) {
           (*CONS)[TAUENERGY] =
               epsrange[0] + 0.5 * B2tilde +
               (B2tilde * stilde_sq - SQ(SdotBtilde)) * 0.5 / (SQ(Wmin + B2tilde));
           (*error_bits)[c2p_errors::TAU_FIX] = true;
         }

      // Check more appropriate limits

      // Etienne 2012 et al. (A56)
      const double safetyfac = 0.999999;
      const double stilde_sq_max = safetyfac * SQ((*CONS)[TAUENERGY] + 1.);

      if (stilde_sq > stilde_sq_max) {
        const double convfac = sqrt(stilde_sq_max / stilde_sq);

        (*CONS)[STILDEX] *= convfac;
        (*CONS)[STILDEY] *= convfac;
        (*CONS)[STILDEZ] *= convfac;

        stilde_sq = stilde_sq_max;
        SdotBtilde *= convfac;

        (*error_bits)[c2p_errors::STILDE_FIX_AH] = true;
      }
    }
  }

 public:
  static inline double limit_tau_and_return_stilde_sq_max(
      cons &CONS, prims &PRIMS, metric &METRIC,
      typename Margherita_C2P_conventions::error_bits_t &error_bits) {
    using namespace Margherita_C2P_conventions;
    const double B2L = METRIC.compute_square_3_gen<BX_CENTER, numcons>(CONS);
    // We assume that the conservatives have already been rescaled!
    typename eos::error_type error;  // TODO How do we handle this?
    constexpr double safetyfac = 1.;
    double rhoL = CONS[RHOSTAR];  // / (safetyfac * limits::lorentz_max);

       const auto epsrange = eos::eps_range__rho_yle_ymu(rhoL, PRIMS[YE], PRIMS[YMU], error);
       // This is (A43) in Etienne et al 2012 https://arxiv.org/pdf/1112.0568.pdf
       if (CONS[TAUENERGY] - 0.5 * B2L / CONS[RHOSTAR] <
           std::min(0., eos::c2p_eps_min)) {
#ifdef DEBUG
         std::cout << "Limit tau" << std::endl;
#endif
         CONS[TAUENERGY] = epsrange[0] + 0.5 * B2L / CONS[RHOSTAR];
         error_bits[c2p_errors::TAU_FIX] = true;
       }

    // For a derivation see Etienne et al 2012:
    // https://arxiv.org/pdf/1112.0568.pdf
    // FIXME :We might want to check this...
    return SQ(CONS[TAUENERGY] + 1.);
  }

  template <typename ArrayType>
  static inline double limit_tau_and_return_stilde_sq_max(
      cons &CONS, prims &PRIMS, ArrayType &list,
      typename Margherita_C2P_conventions::error_bits_t &error_bits) {
    using namespace Margherita_C2P_conventions;

    // We assume that the conservatives have already been rescaled!
    typename eos::error_type error;  // TODO How do we handle this?
    constexpr double safetyfac = 1.2;
    double rhoL = CONS[RHOSTAR] / (safetyfac * limits::lorentz_max);

       const auto epsrange = eos::eps_range__rho_yle_ymu(rhoL, PRIMS[YE], PRIMS[YMU], error);
       // This is (A43) in Etienne et al 2012 https://arxiv.org/pdf/1112.0568.pdf
       if (CONS[TAUENERGY] - 0.5 * list[1] < epsrange[0]) {
#ifdef DEBUG
         std::cout << "Limit tau" << std::endl;
#endif
         CONS[TAUENERGY] = epsrange[0] + 0.5 * list[1];
         error_bits[c2p_errors::TAU_FIX] = true;
       }

    // For a derivation see Etienne et al 2012:
    // https://arxiv.org/pdf/1112.0568.pdf
    // FIXME :We might want to check this...
    return SQ(CONS[TAUENERGY] + 1.);
  }

  // This is necessary since for MHD we need to store stilde_sq but for HD we
  // just need stilde
  // This is (C2) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf

  inline void set_stilde(const double &stilde_sqL) { stilde_sq = stilde_sqL; }

  C2P_MHD_t(cons *__restrict CONS, prims *__restrict PRIMS,
            metric *__restrict METRIC,
            std::bitset<Margherita_C2P_conventions::c2p_errors::NUM_ERRORS>
                *__restrict error_bits)
      : CONS(CONS), PRIMS(PRIMS), METRIC(METRIC), error_bits(error_bits) {
    using namespace Margherita_C2P_conventions;

    B2tilde = METRIC->compute_square_3_gen<BX_CENTER, numcons>(*CONS) /
              (*CONS)[RHOSTAR];

    const double sqrt_rhostar = sqrt((*CONS)[RHOSTAR]);
    (*CONS)[BX_CENTER] /= sqrt_rhostar;
    (*CONS)[BY_CENTER] /= sqrt_rhostar;
    (*CONS)[BZ_CENTER] /= sqrt_rhostar;

    // Need to also compute SdotB
    SdotBtilde = ((*CONS)[BX_CENTER] * (*CONS)[STILDEX] +
                  (*CONS)[BY_CENTER] * (*CONS)[STILDEY] +
                  (*CONS)[BZ_CENTER] * (*CONS)[STILDEZ]);

    inside_BH = (METRIC->sqrtdet > limits::psi6_bh);
    enforce_MHD_BH_limits();
  }

  template <typename ArrayType>
  C2P_MHD_t(cons *__restrict CONS, prims *__restrict PRIMS, ArrayType &list,
            std::bitset<Margherita_C2P_conventions::c2p_errors::NUM_ERRORS>
                *__restrict error_bits)
      : CONS(CONS), PRIMS(PRIMS), METRIC(nullptr), error_bits(error_bits) {
    using namespace Margherita_C2P_conventions;

    B2tilde = list[1];

    const double sqrt_rhostar = sqrt((*CONS)[RHOSTAR]);
    (*CONS)[BX_CENTER] /= sqrt_rhostar;
    (*CONS)[BY_CENTER] /= sqrt_rhostar;
    (*CONS)[BZ_CENTER] /= sqrt_rhostar;

    // Need to also compute SdotB
    SdotBtilde = list[2];

    inside_BH = false;

    enforce_MHD_BH_limits();
  }

  static inline cons compute_conservatives(
      const std::array<double, numprims> &PRIMSL, const metric_c &METRICL,
      const std::array<double, 3> &Bvec) {
    using namespace Margherita_C2P_conventions;
    // Need to compute Lorentz factor and zvec_low
    std::array<double, 3> zvec_lowL =
        METRICL.lower_index<ZVECX, numprims>(PRIMSL);
    const double z2L = PRIMSL[ZVECX] * zvec_lowL[0] +
                       PRIMSL[ZVECY] * zvec_lowL[1] +
                       PRIMSL[ZVECZ] * zvec_lowL[2];
    const double lorentzL = sqrt(1. + z2L);
#ifdef DEBUG
    std::cout << "Lorentz fac: " << lorentzL << std::endl;
#endif
    const double hL = 1. + PRIMSL[EPS] + PRIMSL[PRESSURE] / PRIMSL[RHOB];

    const auto Bvec_low = METRICL.lower_index<0, 3>(Bvec);

    const double B2L =
        Bvec[0] * Bvec_low[0] + Bvec[1] * Bvec_low[1] + Bvec[2] * Bvec_low[2];

    const double Bdotv = METRICL.scalar_Product(Bvec, zvec_lowL) / lorentzL;

    cons CONSL;
    CONSL[RHOSTAR] = METRICL.sqrtdet * PRIMSL[RHOB] * lorentzL;
#pragma unroll
    for (int i = 0; i < 3; ++i) {
      CONSL[STILDEX + i] =
          (CONSL[RHOSTAR] * hL + METRICL.sqrtdet * B2L / lorentzL) *
              zvec_lowL[i] -
          METRICL.sqrtdet * Bdotv * Bvec_low[i];
    }

    CONSL[TAUENERGY] =
        CONSL[RHOSTAR] *
            ((lorentzL - 1.) +
             lorentzL * (PRIMSL[EPS] + PRIMSL[PRESSURE] / PRIMSL[RHOB])) -
        METRICL.sqrtdet * PRIMSL[PRESSURE]  // Hydro
        + METRICL.sqrtdet *
              (B2L - 0.5 * (B2L / SQ(lorentzL) + SQ(Bdotv)));  // MHD

    CONSL[YMUSTAR] = CONSL[RHOSTAR] * PRIMSL[YMU];
    CONSL[YESTAR] = CONSL[RHOSTAR] * PRIMSL[YE];

    CONSL[BX_CENTER] = Bvec[0];
    CONSL[BY_CENTER] = Bvec[1];
    CONSL[BZ_CENTER] = Bvec[2];

    return CONSL;
  }

  inline void compute_conservatives() {
    *CONS = compute_conservatives(*PRIMS, *METRIC);
  }

  // The root finding routine requires us to provide a functor for the
  // residual.. So
  // let´s provide that
  typename formulation<eos>::type operator()(
      typename formulation<eos>::type &X) {
    // This should be an inlined call to the formuation object!
    return formulation<eos>::evaluate(X, CONS, PRIMS, stilde_sq, B2tilde,
                                      SdotBtilde, aux_vars);
  }

  // Update all primitives
  inline void update_primitives(typename formulation<eos>::type &X) {
    using namespace Margherita_helpers;
    using namespace Margherita_C2P_conventions;

    double lorentz;
    formulation<eos>::update_main_primitives(X, CONS, PRIMS, stilde_sq, B2tilde,
                                             SdotBtilde, aux_vars, lorentz);
    double h;

    typename eos::error_type error;  // TODO How do we handle this?
    // We can't just straight away call the pressure interface since
    // eps might be outside of the table range.
    // Now the table routines will automatically enforce those limits, but
    // then TAUENERGY is inconsitent, hence we need to flag this point for
    // recomputation of the conservatives

       auto epsrange = eos::eps_range__rho_yle_ymu((*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);
       // Also limit in the case of exceeding eps_max
       epsrange[1] = min(epsrange[1], limits::eps_max);
       // This is (C27) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf
       if ((*PRIMS)[EPS] < epsrange[0]) {
         (*error_bits)[c2p_errors::EPS_ADJUSTED] = true;
         (*PRIMS)[EPS] = 1.0001 * epsrange[0];
       }
       if ((*PRIMS)[EPS] > epsrange[1]) {
         (*error_bits)[c2p_errors::EPS_ADJUSTED] = true;
         (*PRIMS)[EPS] = 0.9999 * epsrange[1];
       }
       // 4. Compute a by making a pressure call
       // Update all vars here, cs2, temp, etc..
       (*PRIMS)[PRESSURE] = eos::press_h_csnd2_temp_entropy__eps_rho_yle_ymu(
           h, (*PRIMS)[CS2], (*PRIMS)[TEMP],
           (*PRIMS)[ENTROPY],  // out all but temp (inout)
           (*PRIMS)[EPS], (*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);  // in

    // Raise Stilde
    const auto Stilde_up = METRIC->raise_index<STILDEX, numcons>(*CONS);
    const auto hW = h * lorentz;

    //   This is (31) in Noble et al. 2006
    //   https://arxiv.org/pdf/astro-ph/0512420.pdf
    (*PRIMS)[ZVECX] =
        lorentz * (Stilde_up[WVX] + (*CONS)[BX_CENTER] * (SdotBtilde / hW)) /
        (hW + B2tilde);
    (*PRIMS)[ZVECY] =
        lorentz * (Stilde_up[WVY] + (*CONS)[BY_CENTER] * (SdotBtilde / hW)) /
        (hW + B2tilde);
    (*PRIMS)[ZVECZ] =
        lorentz * (Stilde_up[WVZ] + (*CONS)[BZ_CENTER] * (SdotBtilde / hW)) /
        (hW + B2tilde);

    // Do we need to limit the Lorentz factor?
    if (lorentz > limits::lorentz_max) {
      const auto zL = sqrt(SQ(lorentz) - 1.);

      (*error_bits)[c2p_errors::V_MAX_EXCEEDED] = true;

      (*PRIMS)[ZVECX] *= limits::z_max / zL;
      (*PRIMS)[ZVECY] *= limits::z_max / zL;
      (*PRIMS)[ZVECZ] *= limits::z_max / zL;

      // Important we keep RHOSTAR constant so
      // this changes RHOB
      // Why can't we just do this further up when we compute
      // RHOB for the first time? The answer is that we need
      // a selfconsistent solution to obtain the correct velocities
      // from Stilde^2 since we only limit at the very end.
      (*PRIMS)[RHOB] = (*CONS)[RHOSTAR] / limits::lorentz_max;

      // 4. Compute a by making a pressure call
      // Update all vars here, cs2, temp, etc..
         (*PRIMS)[PRESSURE] = eos::press_h_csnd2_temp_entropy__eps_rho_yle_ymu(
             h, (*PRIMS)[CS2], (*PRIMS)[TEMP],
             (*PRIMS)[ENTROPY],  // out all but temp (inout)
             (*PRIMS)[EPS], (*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);  // in

    // code test
      //EOS_Leptonic::c2p_roofinding_call_iteration_masfunconly = EOS_Leptonic::c2p_roofinding_call_iteration_masfunconly + 1;
    }
  }

  inline double error_estimate(double const &hW, double const &lorentz) {
    using namespace Margherita_helpers;
    using namespace Margherita_C2P_conventions;

    const auto z_B2_sq = SQ(hW + B2tilde);
    const auto tau = (*CONS)[TAUENERGY] - (hW - 1.) - 0.5 * B2tilde +
                     (*PRIMS)[PRESSURE] / (*CONS)[RHOSTAR] -
                     0.5 * (B2tilde * stilde_sq - SQ(SdotBtilde)) / z_B2_sq;

    auto const residual = std::fabs(tau / (hW - 1.));

    if (std::isfinite(residual)) {
      return residual;
    } else {
      return 1.e99;
    }
  }

  inline typename formulation<eos>::type invert(double &residual) {
    typename formulation<eos>::type A, B;
    formulation<eos>::compute_braketing_interval(A, B, CONS, PRIMS, stilde_sq,
                                                 B2tilde, SdotBtilde, aux_vars);

    auto res = formulation<eos>::template find_root<>(A, B, *this);
    residual = std::fabs(formulation<eos>::template error<>(res, *this));

    return res;
  }
};
