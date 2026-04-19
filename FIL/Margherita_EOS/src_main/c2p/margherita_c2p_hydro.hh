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

#include "../utils/brent.hh"

template <typename eos>
class C2P_Hydro {
  friend eos;

 private:
  double stilde;
  typedef Margherita_C2P_Limits limits;

  // Compute eps as a function of z
  // This is (C4,C16) in Galeazzi et al 2013.
  // https://arxiv.org/pdf/1306.4953.pdf
  inline double compute_eps_from_Z(const double &Z, const double &lorentz) {
    using namespace Margherita_C2P_conventions;
    return lorentz * (*CONS)[TAUENERGY] - Z * (stilde - Z / (1. + lorentz));
  }

  // Compute the specific enthalpy h= 1 + eps+ press/rho
  // This is (C5) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf
  inline double compute_h(const double &a) {
    using namespace Margherita_C2P_conventions;
    return (1. + (*PRIMS)[EPS]) * (1. + a);
  }

  // Compute the specific enthalpy a= press/(rho*( 1 + eps))
  // This is (C1) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf
  inline double compute_a() {
    using namespace Margherita_helpers;
    using namespace Margherita_C2P_conventions;
    // This is (C1) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf
    const double a =
        (*PRIMS)[PRESSURE] / ((*PRIMS)[RHOB] * (1.0 + (*PRIMS)[EPS]));

    // Enforce the dominant energy condition:
    // This is (C6) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf
    // return min(a, 1.0);
    return a;  // Kastaun actually doesn't enforce the dominant energy
               //  condition at this stage, although the paper says so...
  }

 public:
  static constexpr bool fix_conservatives = true ;
  static constexpr int numcons = 6;
  static constexpr int numprims = 10;

  static double k_max;
  double gamma_lorentz;

  typedef const metric_c metric;
  typedef std::array<double, numcons> cons;
  typedef std::array<double, numprims> prims;

  typename Margherita_C2P_conventions::error_bits_t *error_bits;

  cons *CONS;
  prims *PRIMS;
  const metric_c *METRIC;

  //  static constexpr double kmax = 0.9999975198504966; // lorentz factor 15.
  static constexpr double tol = 1.e-15;  // tolerance on z

  C2P_Hydro(cons *CONS, prims *PRIMS, metric *METRIC,
            typename Margherita_C2P_conventions::error_bits_t *error_bits)
      : CONS(CONS), PRIMS(PRIMS), METRIC(METRIC), error_bits(error_bits){};

  template <typename ArrayType>
  C2P_Hydro(cons *CONS, prims *PRIMS, ArrayType &list,
            typename Margherita_C2P_conventions::error_bits_t *error_bits)
      : CONS(CONS), PRIMS(PRIMS), METRIC(nullptr), error_bits(error_bits){};

  static inline double limit_tau_and_return_stilde_sq_max(
      cons &CONS, prims &PRIMS, metric &METRIC,
      typename Margherita_C2P_conventions::error_bits_t &error_bits) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    // We assume that the conservatives have already been rescaled!
    typename eos::error_type error;  // TODO How do we handle this?
    constexpr double safetyfac = 1.0;
    double rhoL = CONS[RHOSTAR] / (safetyfac * limits::lorentz_max);
    const auto epsrange = eos::eps_range__rho_ye(rhoL, PRIMS[YE], error);

    // This is (A43) in Etienne et al 2012 https://arxiv.org/pdf/1112.0568.pdf
    if (CONS[TAUENERGY] < epsrange[0]) {
      //      if(CONS[TAUENERGY] < 0 ){
      CONS[TAUENERGY] = epsrange[0];
      error_bits[c2p_errors::TAU_FIX] = true;
    }
    // This is (1) in  Deaton et al. 2013: https://arxiv.org/abs/1304.3384
    // For a derivation see Etienne et al 2012:
    // https://arxiv.org/pdf/1112.0568.pdf
    // FIXME :We might want to check this...
    //
    // return CONS[TAUENERGY] * (CONS[TAUENERGY] + 2.) +
    //       (1. - SQ(eos::h_min)) / SQ(METRIC.sqrtdet * CONS[RHOSTAR]);
    return SQ((CONS[TAUENERGY] + 1.) * k_max);
  }

  template <typename ArrayType>
  static inline double limit_tau_and_return_stilde_sq_max(
      cons &CONS, prims &PRIMS, ArrayType &list,
      typename Margherita_C2P_conventions::error_bits_t &error_bits) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    // We assume that the conservatives have already been rescaled!
    typename eos::error_type error;  // TODO How do we handle this?
    constexpr double safetyfac = 1.0;
    double rhoL = CONS[RHOSTAR] / (safetyfac * limits::lorentz_max);
    const auto epsrange = eos::eps_range__rho_ye(rhoL, PRIMS[YE], error);

    // This is (A43) in Etienne et al 2012 https://arxiv.org/pdf/1112.0568.pdf
    if (CONS[TAUENERGY] < epsrange[0]) {
      CONS[TAUENERGY] = epsrange[0];
      error_bits[c2p_errors::TAU_FIX] = true;
    }
    // This is (1) in  Deaton et al. 2013: https://arxiv.org/abs/1304.3384
    // For a derivation see Etienne et al 2012:
    // https://arxiv.org/pdf/1112.0568.pdf
    // FIXME :We might want to check this...
    //
    // return CONS[TAUENERGY] * (CONS[TAUENERGY] + 2.) +
    //       (1. - SQ(eos::h_min)) / SQ(METRIC.sqrtdet * CONS[RHOSTAR]);
    return SQ((CONS[TAUENERGY] + 1.) * k_max);
  }

  // This is necessary since for MHD we need to store stilde_sq but for HD we
  // just need stilde
  // This is (C2) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf
  inline void set_stilde(const double &stilde_sq) { stilde = sqrt(stilde_sq); }

  // Compute the braketing interval of the root
  // This is (C23) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf
  static inline void compute_braketing_interval(double &A, double &B,
                                                const double &k) {
    using namespace Margherita_C2P_conventions;

    A = 0.5 * k / sqrt(1. - 0.25 * SQ(k));
    B = 1.e-6 + k / sqrt(1. - SQ(k));
  }

  // Update all primitives
  inline void update_primitives(const double &z) {
    using namespace Margherita_helpers;
    using namespace Margherita_C2P_conventions;

    bool conservatives_inconsistent = false;

    const double lorentz = sqrt(1.0 + SQ(z));

    // Update rho and eps
    // 2. Update rho
    (*PRIMS)[RHOB] = (*CONS)[RHOSTAR] / lorentz;
    // 3. Compute eps and h
    (*PRIMS)[EPS] = compute_eps_from_Z(z, lorentz);

    typename eos::error_type error;  // TODO How do we handle this?
    // We can't just straight away call the pressure interface since
    // eps might be outside of the table range.
    // Now the table routines will automatically enforce those limits, but
    // then TAUENERGY is inconsitent, hence we need to flag this point for
    // recomputation of the conservatives
    auto epsrange = eos::eps_range__rho_ye((*PRIMS)[RHOB], (*PRIMS)[YE], error);
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
    double h;
    // 4. Compute a by making a pressure call
    // Update all vars here, cs2, temp, etc..
    (*PRIMS)[PRESSURE] = eos::press_h_csnd2_temp_entropy__eps_rho_ye(
        h, (*PRIMS)[CS2], (*PRIMS)[TEMP],
        (*PRIMS)[ENTROPY],  // out all but temp (inout)
        (*PRIMS)[EPS], (*PRIMS)[RHOB], (*PRIMS)[YE], error);  // in

#ifdef DEBUG
    std::cout << "\n";
    std::cout << "PRIMS[PRESSURE] :" << (*PRIMS)[PRESSURE] << std::endl;
    std::cout << "PRIMS[EPS] :" << (*PRIMS)[EPS] << std::endl;
    std::cout << "PRIMS[CS2] :" << (*PRIMS)[CS2] << std::endl;
    std::cout << "PRIMS[ENTROPY] :" << (*PRIMS)[ENTROPY] << std::endl;
    std::cout << "h :" << h << std::endl;
    std::cout << "lorentz final :" << lorentz << std::endl;
#endif

    // Compute h the complicated way again to enforce the dominant energy
    // condition
    h = compute_h(compute_a());

    const auto stilde_up = METRIC->raise_index<STILDEX, numcons>(*CONS);
    (*PRIMS)[ZVECX] = stilde_up[0] / h;
    (*PRIMS)[ZVECY] = stilde_up[1] / h;
    (*PRIMS)[ZVECZ] = stilde_up[2] / h;

    // Do we need to limit the Lorentz factor?
    if (lorentz > limits::lorentz_max) {
      (*error_bits)[c2p_errors::V_MAX_EXCEEDED] = true;

      (*PRIMS)[ZVECX] *= limits::z_max / z;
      (*PRIMS)[ZVECY] *= limits::z_max / z;
      (*PRIMS)[ZVECZ] *= limits::z_max / z;

      // Important we keep RHOSTAR constant so
      // this changes RHOB
      // Why can't we just do this further up when we compute
      // RHOB for the first time? The answer is that we need
      // a selfconsistent solution to obtain the correct velocities
      // from Stilde^2 since we only limit at the very end.
      (*PRIMS)[RHOB] = (*CONS)[RHOSTAR] / limits::lorentz_max;

      // 4. Compute a by making a pressure call
      // Update all vars here, cs2, temp, etc..
      (*PRIMS)[PRESSURE] = eos::press_h_csnd2_temp_entropy__eps_rho_ye(
          h, (*PRIMS)[CS2], (*PRIMS)[TEMP],
          (*PRIMS)[ENTROPY],  // out all but temp (inout)
          (*PRIMS)[EPS], (*PRIMS)[RHOB], (*PRIMS)[YE], error);  // in
    }
  }

  static inline cons compute_conservatives(
      const std::array<double, numprims> &PRIMSL, const metric_c &METRICL) {
    using namespace Margherita_C2P_conventions;
    // Need to compute Lorentz factor and zvec_low
    const auto zvec_lowL = METRICL.lower_index<ZVECX, numprims>(PRIMSL);
    const double z2L = PRIMSL[ZVECX] * zvec_lowL[0] +
                       PRIMSL[ZVECY] * zvec_lowL[1] +
                       PRIMSL[ZVECZ] * zvec_lowL[2];
    const double lorentzL = sqrt(1. + z2L);
    const double hL = 1. + PRIMSL[EPS] + PRIMSL[PRESSURE] / PRIMSL[RHOB];

    cons CONSL;
    CONSL[RHOSTAR] = METRICL.sqrtdet * PRIMSL[RHOB] * lorentzL;
    CONSL[STILDEX] = CONSL[RHOSTAR] * hL * zvec_lowL[0];
    CONSL[STILDEY] = CONSL[RHOSTAR] * hL * zvec_lowL[1];
    CONSL[STILDEZ] = CONSL[RHOSTAR] * hL * zvec_lowL[2];
    CONSL[TAUENERGY] =
        CONSL[RHOSTAR] *
            ((lorentzL - 1.) +
             lorentzL * (PRIMSL[EPS] + PRIMSL[PRESSURE] / PRIMSL[RHOB])) -
        METRICL.sqrtdet * PRIMSL[PRESSURE];
    CONSL[YESTAR] = CONSL[RHOSTAR] * PRIMSL[YE];
    return CONSL;
  }

  inline void compute_conservatives() {
    *CONS = compute_conservatives(*PRIMS, *METRIC);
  }

  // The root finding routine requires us to provide a functor for the
  // residual.. So
  // let´s provide that
  double operator()(const double &X) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    // This and the next line are (C15) in Galeazzi et al 2013.
    // https://arxiv.org/pdf/1306.4953.pdf
    // 1. Compute the Lorentz factor
    const double lorentz = sqrt(1. + SQ(X));
    gamma_lorentz = lorentz;
    // 2. Update rho
    (*PRIMS)[RHOB] = (*CONS)[RHOSTAR] / lorentz;

    typename eos::error_type error;
    const auto epsrange =
        eos::eps_range__rho_ye((*PRIMS)[RHOB], (*PRIMS)[YE], error);
    // 3. Compute eps and h
    (*PRIMS)[EPS] = compute_eps_from_Z(X, lorentz);
    (*PRIMS)[EPS] = min(epsrange[1], max(epsrange[0], (*PRIMS)[EPS]));

    // 4. Compute a by making a pressure call
    (*PRIMS)[PRESSURE] = eos::press_temp__eps_rho_ye(
        (*PRIMS)[TEMP],                                       // inout
        (*PRIMS)[EPS], (*PRIMS)[RHOB], (*PRIMS)[YE], error);  // in
    // TODO Need error handling!

    const double a = compute_a();

#ifdef DEBUG
    std::cout << "\n";
    std::cout << "PRIMS[PRESSURE] :" << (*PRIMS)[PRESSURE] << std::endl;
    std::cout << "PRIMS[EPS] :" << (*PRIMS)[EPS] << std::endl;
    std::cout << "lorentz fac : " << lorentz << std::endl;
    std::cout << "residual :" << X - stilde / compute_h(a) << std::endl;
#endif

    // This is (C22) in Galeazzi et al 2013: https://arxiv.org/pdf/1306.4953.pdf
    // Finally: Compute residual
    return X - stilde / compute_h(a);
  }

  inline double invert(double &residual) {
    using namespace Margherita_C2P_conventions;

    // This is (C2) in Galeazzi et al 2013. https://arxiv.org/pdf/1306.4953.pdf
    // I am not sure this will work for MHD, but will try
    // FIXME: First try without the limiting
    // const double k= min(sqrt(stilde_sq)/(CONS[TAUENERGY]+1.), kmax);
    double k = stilde / ((*CONS)[TAUENERGY] + 1.);
#ifdef DEBUG
    std::cout << "k, larger" << k << " , " << (k > k_max) << std::endl;
#endif
    // Limit k
    k = std::min(k_max, k);

    // We need a braketing
    double A, B;
    // This interface might have to change for MHD
    compute_braketing_interval(A, B, k);

    // Step 3: the stage is set, let's try to find a solution (This should
    // always
    // converge to some possibly unphysical value
    const auto res = zero_brent<>(A, B, tol, *this);

    // Potentially dangerous!!
    // const auto norm_fac = (std::fabs(res) > 1.e-40) ? 1./res : 1.;
    // residual = std::fabs((*this)(res))*norm_fac;
    residual = std::fabs((*this)(res));

    return res;
  };
};
