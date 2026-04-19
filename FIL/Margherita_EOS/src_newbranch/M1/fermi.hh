/*
 * =====================================================================================
 *
 *       Filename:  fermi.hh
 *
 *    Description:  Approximate formulae for Fermi Integrals
 *    		    See Takahashi et al. 1978
 *
 *        Version:  1.0
 *        Created:  13/05/2017 19:31:13
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#ifndef __HH__FERMI
#define __HH__FERMI

#include "helpers.hh"

constexpr double eta_min = 1.e-3;

template <int order>
class Fermi_Dirac {};

template <>
class Fermi_Dirac<0> {
 public:
  static inline double get(const double& eta) {
    if (eta > eta_min)
      return eta + log(1. + exp(-eta));
    else
      return log(1. + exp(eta));
  }
};
template <>
class Fermi_Dirac<1> {
 public:
  static inline double get(const double& eta) {
    using namespace m1_helpers ;
    if (eta > eta_min)
      return horner(eta, 1.6449, 0., 0.5) / (1. + exp(-1.6855 * eta));
    else
      return exp(eta) / (1. + 0.2159 * exp(0.8857 * eta));
  }

  // not to be used by the user!
  static inline double no_exp_inv(const double& eta) {
    return (1. + 0.2159 * exp(0.8857 * eta));
  }
};
template <>
class Fermi_Dirac<2> {
 public:
  static inline double get(const double& eta) {
    using namespace m1_helpers ;
    if (eta > eta_min)
      return horner(eta, 0., 3.2899, 0., 1. / 3.) / (1. - exp(-1.8246 * eta));
    else
      return 2. * exp(eta) / (1. + 0.1092 * exp(0.8908 * eta));
  }
  // not to be used by the user!
  static inline double no_exp_inv(const double& eta) {
    return 0.5 * (1. + 0.1092 * exp(0.8908 * eta));
  }
};
template <>
class Fermi_Dirac<3> {
 public:
  static inline double get(const double& eta) {
    using namespace m1_helpers ;
    if (eta > eta_min)
      return horner(eta, 11.3644, 0., 4.9348, 0., 0.25) /
             (1. + exp(-1.9039 * eta));
    else
      return 6. * exp(eta) / (1. + 0.0559 * exp(0.9069 * eta));
  }
  // not to be used by the user!
  static inline double no_exp_inv(const double& eta) {
    return 1. / 6. * (1. + 0.0559 * exp(0.9069 * eta));
  }
};

template <>
class Fermi_Dirac<4> {
 public:
  static inline double get(const double& eta) {
    using namespace m1_helpers ;
    if (eta > eta_min)
      return horner(eta, 0., 45.4576, 0., 6.5797, 0., 0.2) /
             (1. - exp(-1.9484 * eta));
    else
      return 24. * exp(eta) / (1. + 0.0287 * exp(0.9257 * eta));
  }
  // not to be used by the user!
  static inline double no_exp_inv(const double& eta) {
    return 1. / 24. * (1. + 0.0287 * exp(0.9257 * eta));
  }
};

template <>
class Fermi_Dirac<5> {
 public:
  static inline double get(const double& eta) {
    using namespace m1_helpers ;
    if (eta > eta_min)
      return horner(eta, 236.5323, 0., 113.6439, 0., 8.2247, 0., 1. / 6.) /
             (1. + exp(-1.9727 * eta));
    else
      return 120. * exp(eta) / (1. + 0.0147 * exp(0.9431 * eta));
  }
  // not to be used by the user!
  static inline double no_exp_inv(const double& eta) {
    return 1. / 120. * (1. + 0.0147 * exp(0.9431 * eta));
  }
};

template <int N, int M>
class Fermi_Dirac_Ratio {
 public:
  static inline double get(const double& eta) {
    if (eta > eta_min)
      return Fermi_Dirac<N>::get(eta) / Fermi_Dirac<M>::get(eta);
    else
      return Fermi_Dirac<M>::no_exp_inv(eta) / Fermi_Dirac<N>::no_exp_inv(eta);
  }
};

#endif
