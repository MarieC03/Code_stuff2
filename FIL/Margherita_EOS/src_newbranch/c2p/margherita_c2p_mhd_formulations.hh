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

#include "../Cold/cold_pwpoly.hh"
#include "../Cold/cold_pwpoly_implementation.hh"
//#include "../utils/ND-newton-raphson.cc"
#include "../utils/NewtonRaphson.hh"
#include "../utils/brent.hh"

template <typename eos>
class C2P_MHD_Palenzuela_f {
  friend eos;

 public:
  static constexpr bool fix_conservatives = true ;
  //static constexpr int numcons = 9;
  //static constexpr int numprims = 10;
  static constexpr int numcons = 10;
  static constexpr int numprims = 11;
  static constexpr int num_aux = 1;

  typedef const metric_c metric;
  typedef std::array<double, numcons> cons;
  typedef std::array<double, numprims> prims;

  typedef double type;

  static constexpr double tol = 1.e-15;  // tolerance on hW

  static inline void update_main_primitives(
      const type &hW, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars, double &lorentz) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    // 1. Compute the Lorentz factor
    // This is step 1. in the appendix of Palenzuela et al. 2015

    double lorentz_inv_sqL =
        1. - (stilde_sq + (2. * hW + B2tilde) * SQ(SdotBtilde / hW)) /
                 (SQ(hW + B2tilde));
    // Enforce limits on the Lorentz factor to ensure that we never get
    // superluminal
    lorentz_inv_sqL =
        min(1.0, max(1. / SQ(1.2 * Margherita_C2P_Limits::lorentz_max),
                     lorentz_inv_sqL));
    const double lorentz_inv = sqrt(lorentz_inv_sqL);

    // 2. Update rho
    (*PRIMS)[RHOB] = (*CONS)[RHOSTAR] * lorentz_inv;

    lorentz = 1. / lorentz_inv;

    // 3. Compute eps
    // Compute eps as a function of z
    // This is step 3. in the appendix of Palenzuela et al. 2015
    // https://arxiv.org/pdf/1505.01607.pdf
    const double epsL =
        lorentz - 1. + hW * lorentz_inv * (1. - SQ(lorentz)) +
        lorentz * ((*CONS)[TAUENERGY] - B2tilde + 0.5 * SQ(SdotBtilde / hW) +
                   0.5 * B2tilde * lorentz_inv_sqL);

    typename eos::error_type error;  // TODO How do we handle this?
    // We can't just straight away call the pressure interface since
    // eps might be outside of the table range.
    // Now the table routines will automatically enforce those limits, but
    // then TAUENERGY is inconsitent, hence we need to flag this point for
    // recomputation of the conservatives
       const auto epsrange =
           eos::eps_range__rho_yle_ymu((*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);
       (*PRIMS)[EPS] = min(epsrange[1], max(epsrange[0], epsL));
  }

  // This is the main iteration of Palenzuela et al. 2015
  // https://arxiv.org/pdf/1505.01607.pdf

  static inline type evaluate(const type &hW, cons *__restrict CONS,
                              prims *__restrict PRIMS, const double &stilde_sq,
                              const double &B2tilde, const double &SdotBtilde,
                              double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    // This is the algorithm given in appendix A of Palenzuela et al. 2015
    double lorentz;

    // Compute epsilon
    update_main_primitives(hW, CONS, PRIMS, stilde_sq, B2tilde, SdotBtilde,
                           aux_vars, lorentz);

    // 4. Compute a by making a pressure call
    typename eos::error_type error;

       (*PRIMS)[PRESSURE] = eos::press_temp__eps_rho_yle_ymu(    // out
           (*PRIMS)[TEMP],                                       // inout
           (*PRIMS)[EPS], (*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);  // in
       // TODO Need error handling!

    // Compute the enthalpy
    const double a =
        (*PRIMS)[PRESSURE] / ((*PRIMS)[RHOB] * (1.0 + (*PRIMS)[EPS]));
    // We need something like the dominant energy condition for MHD
    // See here: Pimentel et al. 2016
    // https://arxiv.org/pdf/1612.03299.pdf

    const double h = (1. + (*PRIMS)[EPS]) * (1. + a);

#ifdef DEBUG
    std::cout << "\n";
    std::cout << "PRIMS[PRESSURE] :" << (*PRIMS)[PRESSURE] << std::endl;
    std::cout << "PRIMS[EPS] :" << (*PRIMS)[EPS] << std::endl;
    std::cout << "residual :" << hW - h * lorentz << std::endl;
    ;
    std::cout << "hW :" << hW << std::endl;
    std::cout << "W :" << lorentz << std::endl;
#endif
    // Return the residual
    return hW - h * lorentz;
  }

  // Compute the braketing interval of the root
  // This is (A8) in the appendix of Palenzuela et al. 2015
  // https://arxiv.org/pdf/1505.01607.pdf

  static inline void compute_braketing_interval(
      double &A, double &B, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;

    // Palenzuela says:
    typename eos::error_type error;
    constexpr double safety_fac = 1.2;
    double rho_min =
        (*CONS)[RHOSTAR] / (safety_fac * Margherita_C2P_Limits::lorentz_max);
    double rho_max = (*CONS)[RHOSTAR];
    B = 2. + 2. * (*CONS)[TAUENERGY] - B2tilde;
    const double Stilde_bound =
        (B2tilde > 1.e-40) ? SdotBtilde / sqrt(B2tilde) : 0.;
    A = max(1. + (*CONS)[TAUENERGY] - B2tilde, Stilde_bound);

#ifdef DEBUG
    std::cout << "A, B: " << A << " , " << B << std::endl;
#endif
  }

  template <typename F_t>
  static inline double find_root(const double &A, const double &B, F_t &F) {
    return zero_brent<>(A, B, tol, F);
  }

  template <typename F_t>
  static inline double error(double &Z, F_t &F) {
    using namespace Margherita_C2P_conventions;
    // const auto res= std::fabs(F(Z))/Z;
    const auto h =
        (1. + (*F.PRIMS)[EPS] + (*F.PRIMS)[PRESSURE] / (*F.PRIMS)[RHOB]);
    F.gamma_lorentz = Z / h;
#ifdef DEBUG
    std::cout << "Error estimate: " << F.error_estimate(Z, F.gamma_lorentz)
              << std::endl;
#endif
    const auto res = F.error_estimate(Z, F.gamma_lorentz);
    return res;
  }
};

template <typename eos>
class C2P_MHD_Entropy_Fix_f {
  friend eos;

 public:
  static constexpr bool fix_conservatives = true ;
  //static constexpr int numcons = 10;
  //static constexpr int numprims = 10;
  static constexpr int numcons = 11;
  static constexpr int numprims = 11;
  static constexpr int num_aux = 1;

  typedef const metric_c metric;
  typedef std::array<double, numcons> cons;
  typedef std::array<double, numprims> prims;

  typedef double type;

  static constexpr double tol = 1.e-15;  // tolerance on hW

 private:
  // This is step 1. in the appendix of Palenzuela et al. 2015
  // https://arxiv.org/pdf/1505.01607.pdf

  static inline double compute_lorentz_factor(const double &stilde_fluid_sq,
                                              const double &h_cold) {
    using namespace Margherita_C2P_conventions;
    return sqrt(1. + stilde_fluid_sq / SQ(h_cold));
  }

 public:
  // Compute the braketing interval of the root
  // This is (A59) in the appendix of Etienne et al. 2012
  // https://arxiv.org/pdf/1112.0568.pdf

  static inline void compute_braketing_interval(
      double &A, double &B, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;

    //    const double lorentz = compute_lorentz_factor(stilde_sq, aux_vars[0]);

    A = (*CONS)[RHOSTAR] / (1.2 * Margherita_C2P_Limits::lorentz_max);
    B = (*CONS)[RHOSTAR];

    aux_vars[0] = eos::c2p_h_min;

#ifdef DEBUG
    std::cout << "A, B, h_min: " << A << " , " << B << " , " << aux_vars[0]
              << std::endl;
#endif
  }

  static inline void update_main_primitives(
      const type &X, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars, double &lorentz) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;

    // 0. Update rho
    (*PRIMS)[RHOB] = X;
    // 1. Compute h_cold
    // 4. Compute a by making a pressure call
    typename eos::error_type error;
    // auto entropyrange = eos::entropy_range__rho_ye((*PRIMS)[RHOB],
    // (*PRIMS)[YE],error);
    // 3. Compute eps and h
    // double entL =
    // min(entropyrange[1],max(entropyrange[0],(*PRIMS)[ENTROPY]));
       double entL = (*PRIMS)[ENTROPY];
       (*PRIMS)[PRESSURE] = eos::press_h_csnd2_temp_eps__entropy_rho_yle_ymu(
           aux_vars[0], (*PRIMS)[CS2], (*PRIMS)[TEMP], (*PRIMS)[EPS], entL,
           (*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);  // in
       // TODO Need error handling!

    double A = (B2tilde > 1.e-40) ? SQ(SdotBtilde) / B2tilde : 0.;
    double B = stilde_sq;

#ifdef DEBUG
    std::cout << "Sfluid: A, B: " << A << " , " << B << std::endl;
    std::cout << "Entropy: " << entL << std::endl;

#endif

    const auto func = [&](const double sl) {
      // This is equation (A60) in the appendix of Etienne et al. 2012
      // Compute W (not the lorentz factor!)
      // https://arxiv.org/pdf/1112.0568.pdf
      const auto W = sqrt(sl + SQ(aux_vars[0]));

      // Compute Stilde_fluid
      // This is equation (A61) in the appendix of Etienne et al. 2012
      // https://arxiv.org/pdf/1112.0568.pdf
      const auto s_fl2 =
          (SQ(W) * stilde_sq + SQ(SdotBtilde) * (B2tilde + 2. * W)) /
          SQ(W + B2tilde);
#ifdef DEBUG
      std::cout << "\n";
      std::cout << "sfluid : " << sl << std::endl;
      std::cout << "sfluid residual :" << s_fl2 - sl << std::endl;
#endif
      return s_fl2 - sl;
    };  // END lambda function

    const auto stilde_fl_sq = zero_brent<>(A, B, 1.e-14, func);

    // 1. Compute the Lorentz factor
    lorentz = compute_lorentz_factor(stilde_fl_sq, aux_vars[0]);
  }

  // The root finding routine requires us to provide a functor for the
  // residual.. So
  // let´s provide that

  static inline type evaluate(const type &X, cons *__restrict CONS,
                              prims *__restrict PRIMS, const double &stilde_sq,
                              const double &B2tilde, const double &SdotBtilde,
                              double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;

    double lorentz;
    update_main_primitives(X, CONS, PRIMS, stilde_sq, B2tilde, SdotBtilde,
                           aux_vars, lorentz);

#ifdef DEBUG
    std::cout << "\n";
    std::cout << "residual :" << X - (*CONS)[RHOSTAR] / lorentz << std::endl;
    ;
#endif
    // Return the residual
    return X - (*CONS)[RHOSTAR] / lorentz;
  }

  template <typename F_t>
  static inline double find_root(const double &A, const double &B, F_t &F) {
    return zero_brent<>(A, B, tol, F);
  }

  template <typename F_t>
  static inline double error(double &Z, F_t &F) {
    using namespace Margherita_C2P_conventions;
    const auto res= std::fabs(F(Z))/Z;
    return res;
  }
};



template <typename eos>
class C2P_MHD_Entropy_Fix_New_f {
  friend eos;

 public:
  static constexpr bool fix_conservatives = true ;
  //static constexpr int numcons = 10;
  //static constexpr int numprims = 10;
  static constexpr int numcons = 11;
  static constexpr int numprims = 11;
  static constexpr int num_aux = 1;

  typedef const metric_c metric;
  typedef std::array<double, numcons> cons;
  typedef std::array<double, numprims> prims;

  typedef double type;

  static constexpr double tol = 1.e-15;  // tolerance on hW

 public:
  // Compute the braketing interval of the root
  // This is (A59) in the appendix of Etienne et al. 2012
  // https://arxiv.org/pdf/1112.0568.pdf

  static inline void compute_braketing_interval(
      double &A, double &B, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;

    //    const double lorentz = compute_lorentz_factor(stilde_sq, aux_vars[0]);

    A = 1.;
    B = (1.2 * Margherita_C2P_Limits::lorentz_max);

#ifdef DEBUG
    std::cout << "A, B" << A << " , " << B << std::endl;
#endif
  }

  static inline void update_main_primitives(
      const type &X, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars, double &lorentz) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;

    lorentz = X;
    // 0. Update rho
    (*PRIMS)[RHOB] = (*CONS)[RHOSTAR] / lorentz;
    // 1. Compute h_cold
    // 4. Compute a by making a pressure call
    typename eos::error_type error;
    // auto entropyrange = eos::entropy_range__rho_ye((*PRIMS)[RHOB],
    // (*PRIMS)[YE],error);
    // 3. Compute eps and h
    // double entL =
    // min(entropyrange[1],max(entropyrange[0],(*PRIMS)[ENTROPY]));
       double entL = (*PRIMS)[ENTROPY];
       (*PRIMS)[PRESSURE] = eos::press_h_csnd2_temp_eps__entropy_rho_yle_ymu(
           aux_vars[0], (*PRIMS)[CS2], (*PRIMS)[TEMP], (*PRIMS)[EPS], entL,
           (*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);  // in


  }

  // The root finding routine requires us to provide a functor for the
  // residual.. So
  // let´s provide that

  static inline type evaluate(const type &X, cons *__restrict CONS,
                              prims *__restrict PRIMS, const double &stilde_sq,
                              const double &B2tilde, const double &SdotBtilde,
                              double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;

    double lorentz;
    update_main_primitives(X, CONS, PRIMS, stilde_sq, B2tilde, SdotBtilde,
                           aux_vars, lorentz);

    auto const hW = aux_vars[0]*lorentz;
    auto const lorentz_sq = lorentz*lorentz;
    auto const residual = -lorentz_sq*stilde_sq 
      +(lorentz_sq-1.)*SQ(B2tilde + hW) - lorentz_sq*(B2tilde + 2.*hW)*SQ(SdotBtilde/hW);

    // Enforce limits on the Lorentz factor to ensure that we never get
    // superluminal

#ifdef DEBUG
    std::cout << "\n";
    std::cout << "lorentz :" << lorentz << std::endl;
    std::cout << "residual :" << residual << std::endl;
    ;
#endif
    // Return the residual
    return residual;
  }

  template <typename F_t>
  static inline double find_root(const double &A, const double &B, F_t &F) {
    return zero_brent<>(A, B, tol, F);
  }

  template <typename F_t>
  static inline double error(double &Z, F_t &F) {
    using namespace Margherita_C2P_conventions;
    const auto h =
        (1. + (*F.PRIMS)[EPS] + (*F.PRIMS)[PRESSURE] / (*F.PRIMS)[RHOB]);
    const auto res= std::fabs(F(Z))/SQ(F.B2tilde + h*Z);
    return res;
  }
};

/////////////////////////////////////////////////////////////////////////////////
//
//  Algorithm by W.Newman and N.Hamlin
//
//  "PRIMITIVE VARIABLE DETERMINATION IN CONSERVATIVE RELATIVISTIC
//   MAGNETOHYDRODYNAMIC SIMULATIONS"
//   SIAM J. SCI. COMPUT. Vol. 36, No. 4, pp. B661–B683
//
//   This method uses a fixed point iteration on the pressure.
//   Once the pressure is fixed the inversion is reduced to finding the
//   root of a cubic polynomial for which the analytic solution is known.
//   The algorithm automatically selects the correct root and converges
//   for almost all initial pressure guesses. To speed up the computation
//   and Aitken delta-squared process is used.
//
//   !Be aware that this algorithm can be slow!
//
/////////////////////////////////////////////////////////////////////////////////

template <typename eos>
class C2P_MHD_Newman_f {
  friend eos;

 public:
  static constexpr bool fix_conservatives = true ;
  //static constexpr int numcons = 9;
  //static constexpr int numprims = 10;
  static constexpr int numcons = 10;
  static constexpr int numprims = 11;
  static constexpr int num_aux = 5;

  typedef const metric_c metric;
  typedef std::array<double, numcons> cons;
  typedef std::array<double, numprims> prims;

  typedef double type;

  static constexpr double tol = 1.e-14;  // tolerance on hW

  static inline void update_main_primitives(
      const type &hW, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars, double &lorentz) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    // 1. Compute the Lorentz factor

    auto const vp_sq = (stilde_sq + (2. * hW + B2tilde) * SQ(SdotBtilde / hW));

    auto const Z_sq = std::fabs(vp_sq / ((SQ(hW + B2tilde) - vp_sq)));

    // Enforce limits on the Lorentz factor to ensure that we never get
    // superluminal

    auto const lorentz_sq =
        min(1. + Z_sq, SQ(1.2 * Margherita_C2P_Limits::lorentz_max));

    lorentz = sqrt(lorentz_sq);
    auto lorentz_inv = 1. / lorentz;

    // 2. Update rho
    (*PRIMS)[RHOB] = (*CONS)[RHOSTAR] * lorentz_inv;

    // 3. Compute eps from hW
    const double epsL =
        (hW - lorentz) * lorentz_inv - (*PRIMS)[PRESSURE] / (*PRIMS)[RHOB];
    //    const double epsL =
    //        lorentz - 1. + hW * lorentz_inv * (1. - SQ(lorentz)) +
    //        lorentz * ((*CONS)[TAUENERGY] - B2tilde + 0.5 * SQ(SdotBtilde /
    //        hW) +
    //                   0.5 * B2tilde * lorentz_inv_sqL);

    typename eos::error_type error;  // TODO How do we handle this?
       const auto epsrange =
           eos::eps_range__rho_yle_ymu((*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);
       (*PRIMS)[EPS] = min(epsrange[1], max(epsrange[0], epsL));
  }

  static inline type evaluate(type &hW, cons *__restrict CONS,
                              prims *__restrict PRIMS, const double &stilde_sq,
                              const double &B2tilde, const double &SdotBtilde,
                              double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    double lorentz;

    // This is Newman&Hamlin's variable a which contains a term p/D, so this
    // needs  to be corrected as a new pressure guess is available
    aux_vars[0] += (aux_vars[4] - aux_vars[3]) / (*CONS)[RHOSTAR];

    // phi is a root of the cubic polynomial (E1)^3 - a * (E1)^2 + d  = 0
    // Step 2. in section 8.
    double phi = acos(sqrt(27. / 4. * aux_vars[1] / aux_vars[0]) / aux_vars[0]);
    double E1 = aux_vars[0] * 1. / 3. -
                2. / 3. * aux_vars[0] * cos(2. / 3. * phi + 2. / 3. * M_PI);

    // Note that E1 = hW + B^2/D
    hW = E1 - B2tilde;

    // Compute epsilon
    update_main_primitives(hW, CONS, PRIMS, stilde_sq, B2tilde, SdotBtilde,
                           aux_vars, lorentz);

    // Rotate
    aux_vars[2] = aux_vars[3];
    aux_vars[3] = aux_vars[4];

    // 4. Compute a by making a pressure call
    typename eos::error_type error;
       (*PRIMS)[PRESSURE] = eos::press_temp__eps_rho_yle_ymu(    // out
           (*PRIMS)[TEMP],                                       // inout
           (*PRIMS)[EPS], (*PRIMS)[RHOB], (*PRIMS)[YE], (*PRIMS)[YMU], error);  // in

    //(5.9) in Newman & Hamlin 2014
    // We need to ensure a minimum pressure to guarantee the existence of a
    // root!
    double P_over_Dmin = std::cbrt(27. / 4. * aux_vars[1]) - 0.5 * B2tilde -
                         (*CONS)[TAUENERGY] - 1.;
    (*PRIMS)[PRESSURE] =
        max((*PRIMS)[PRESSURE], (*CONS)[RHOSTAR] * P_over_Dmin);

    aux_vars[4] = (*PRIMS)[PRESSURE];

#ifdef DEBUG
    std::cout << "\n";
    std::cout << "PRIMS[PRESSURE] :" << (*PRIMS)[PRESSURE] << std::endl;
    std::cout << "PRIMS[EPS] :" << (*PRIMS)[EPS] << std::endl;
    std::cout << "residual :"
              << std::fabs(aux_vars[4] - aux_vars[3]) / std::fabs(aux_vars[3])
              << std::endl;
    ;
    std::cout << "hW :" << hW << std::endl;
    std::cout << "W :" << lorentz << std::endl;
#endif

    // Return the residual
    return aux_vars[4] - aux_vars[3];
  }

  // Compute the braketing interval of the root
  // This is (A8) in the appendix of Palenzuela et al. 2015
  // https://arxiv.org/pdf/1505.01607.pdf

  static inline void compute_braketing_interval(
      double &A, double &B, cons *__restrict CONS, prims *__restrict PRIMS,
      const double &stilde_sq, const double &B2tilde, const double &SdotBtilde,
      double *__restrict aux_vars) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;

    aux_vars[1] = 0.5 * (stilde_sq * B2tilde - SQ(SdotBtilde));

    double P_over_Dmin = std::cbrt(27. / 4. * aux_vars[1]) - 0.5 * B2tilde -
                         (*CONS)[TAUENERGY] - 1.;
    (*PRIMS)[PRESSURE] =
        max((*PRIMS)[PRESSURE], (*CONS)[RHOSTAR] * P_over_Dmin);

    aux_vars[4] = (*PRIMS)[PRESSURE];
    aux_vars[3] = (*PRIMS)[PRESSURE];
    aux_vars[2] = (*PRIMS)[PRESSURE];

    aux_vars[0] = (*CONS)[TAUENERGY] + 1. + 0.5 * B2tilde +
                  aux_vars[4] / (*CONS)[RHOSTAR];
  }

  template <typename F_t>
  static inline double find_root(const double &A, const double &B, F_t &F) {
    using namespace Margherita_C2P_conventions;
    using namespace Margherita_helpers;
    constexpr size_t iterations = 300;
    constexpr size_t eps = 0;

    double hW;

    bool fail = false;
    size_t iter = 0;
    double error = 1.e99;
    while(std::fabs(error) > tol*F.aux_vars[3] && iter < iterations){
     F(hW);
     error = F(hW);
     if(iter > iterations-2) fail = true;

     double lipshitz = (F.aux_vars[4] - F.aux_vars[3])/(F.aux_vars[3]-F.aux_vars[2]);

#ifdef DEBUG
     std::cout<< "Fixed point error: " <<error/F.aux_vars[3] <<std::endl;
     std::cout<< "Lipshitz: " <<lipshitz <<std::endl;
#endif
/*
     if(std::fabs(lipshitz)<1. && std::fabs(error) < eps*F.aux_vars[3]){
	auto tmp = F.aux_vars[4];
        F.aux_vars[4]=F.aux_vars[3] + (F.aux_vars[4] - F.aux_vars[3])/(1.- lipshitz);

	(*F.PRIMS)[PRESSURE] = F.aux_vars[4];
        F.aux_vars[2]=F.aux_vars[3];
        F.aux_vars[3]=tmp;

//	error = F.aux_vars[4] - F.aux_vars[3];
#ifdef DEBUG
     std::cout<< "Aitken error: " <<std::fabs(F.aux_vars[4] - F.aux_vars[3])/F.aux_vars[3] <<std::endl;
#endif
 
     }
*/
	 iter++;
      }
#ifdef DEBUG
    std::cout << "Inversion took " << iter << " iterations." << std::endl;
#endif

    if (fail) {
      F.aux_vars[2] = 1.e80;
      F.aux_vars[3] = 1.e85;
      F.aux_vars[4] = 1.e90;
    }

    F(hW);
    return hW;
  }

  template <typename F_t>
  static inline double error(double &Z, F_t &F) {
    using namespace Margherita_C2P_conventions;
    const auto h =
        (1. + (*F.PRIMS)[EPS] + (*F.PRIMS)[PRESSURE] / (*F.PRIMS)[RHOB]);
    F.gamma_lorentz = Z / h;
#ifdef DEBUG
    std::cout << "Error estimate: " << F.error_estimate(Z, F.gamma_lorentz)
              << std::endl;
#endif

    //    return std::fabs(F.aux_vars[4] -
    //    F.aux_vars[3])/std::fabs(F.aux_vars[3]);
    return F.error_estimate(Z, F.gamma_lorentz);
  }
};
