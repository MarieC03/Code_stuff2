/*
 * =====================================================================================
 *
 *       Filename:  leakage.hh
 *
 *    Description:  Computes neutrino emission in a leakage scheme
 *    		    based on the GR1D Leakage Scheme by Evan O'Connor
 *    		    and the (closed source) WeakRates codes by Filippo Galeazzi
 *    		    that both implement a simplified Neutrino Leakage scheme by
 * Ruffert et al 1998 and Rosswog, Liebendörfer 2006 et al.
 *
 *        Version:  1.0
 *        Created:  14/05/2017 01:03:33
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#include "../margherita.hh"
#include "fermi.hh"
#include "fugacities.hh"
#include "helpers.hh"

#ifndef __HH_LEAKAGE
#define __HH_LEAKAGE

class Leakage {
 private:
  enum { NUE = 0, NUA, NUX };

  template <int N>
  using FD = Fermi_Dirac<N>;
  template <int N, int M>
  using FDR = Fermi_Dirac_Ratio<N, M>;

  using FG = Fugacities;

  const std::array<double, 3> g_nu{{1., 1., 4.}};

 public:
  bool brems = true;

  Fugacities F;

  std::array<double, 3> Q{{0, 0, 0}};           // MeV/s/cm^3
  std::array<double, 3> Q_diff_inv{{0, 0, 0}};  // MeV/s/cm^3
  std::array<double, 3> Q_eff{{0, 0, 0}};       // MeV/s/cm^3
  std::array<double, 3> R{{0, 0, 0}};           // number/s/cm^3
  std::array<double, 3> R_diff_inv{{0, 0, 0}};  // number/s/cm^3
  std::array<double, 3> R_eff{{0, 0, 0}};       // number/s/cm^3

  std::array<double, 3> n_nu{{0, 0, 0}};
  std::array<double, 3> en_nu{{0, 0, 0}};

  std::array<double, 3> kappa_n{{0, 0, 0}};  // 1/cm
  std::array<double, 3> kappa_e{{0, 0, 0}};  // 1/cm

  std::array<double, 3> kappa_n_abs{{0, 0, 0}};  // 1/cm
  std::array<double, 3> kappa_e_abs{{0, 0, 0}};  // 1/cm

  inline void calc_emission() {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    auto block_factor = std::array<double, 3>{{1, 1, 1}};

    // BETA DECAY !

    // Ruffert et al. (B3)
    block_factor[NUE] =
        1. + exp(F.eta[NUE] - FDR<5, 4>::get(F.eta[FG::ELECTRON]));

    // Ruffert et al. (B4)
    block_factor[NUA] =
        1. + exp(F.eta[NUA] - FDR<5, 4>::get(-F.eta[FG::ELECTRON]));

    // Ruffert et al. (B1)
    R[NUE] =
        beta * F.eta_pn * int_pow<5>(F.temp) * FD<4>::get(F.eta[FG::ELECTRON]);
    // Ruffert et al. (B2)
    R[NUA] =
        beta * F.eta_np * int_pow<5>(F.temp) * FD<4>::get(-F.eta[FG::ELECTRON]);

    // Ruffert et al. (B15)
    Q[NUE] =
        beta * F.eta_pn * int_pow<6>(F.temp) * FD<5>::get(F.eta[FG::ELECTRON]);
    Q[NUA] =
        beta * F.eta_np * int_pow<6>(F.temp) * FD<5>::get(-F.eta[FG::ELECTRON]);

// This is not taken into account by Evan O'Connor
#pragma unroll
    for (int i = 0; i < R.size(); ++i) {
      R[i] /= block_factor[i];
      Q[i] /= block_factor[i];

      if (!std::isfinite(Q[i])) Q[i] = 0;
      if (!std::isfinite(R[i])) R[i] = 0;
    }

#ifdef DEBUG

    std::cout << F.eta_np << std::endl;

    std::cout << "Beta decay" << std::endl;
    std::cout << "Q : ";
    for (auto &q : Q) std::cout << q << " , ";
    std::cout << std::endl;
    std::cout << "R : ";
    for (auto &q : R) std::cout << q << " , ";
    std::cout << std::endl;
#endif

// ELECTRON POSITRON PAIR ANNIHILATION

// Blocking factors
// Ruffert et al. (B9)
#pragma unroll
    for (int i = NUE; i <= NUX; ++i) {
      block_factor[i] =
          1. + exp(F.eta[i] - 0.5 * (FDR<4, 3>::get(F.eta[FG::ELECTRON]) +
                                     FDR<4, 3>::get(-F.eta[FG::ELECTRON])));
    }

    constexpr auto eps_const = 8. * pi / int_pow<3>(hc_mevcm);

    // Ruffert et al. (B5)
    const auto eps_m =
        eps_const * int_pow<4>(F.temp) * FD<3>::get(F.eta[FG::ELECTRON]);
    const auto eps_p =
        eps_const * int_pow<4>(F.temp) * FD<3>::get(-F.eta[FG::ELECTRON]);

    // Ruffert et al. (B6)
    const auto eps_tilde_m =
        eps_const * int_pow<5>(F.temp) * FD<4>::get(F.eta[FG::ELECTRON]);
    const auto eps_tilde_p =
        eps_const * int_pow<5>(F.temp) * FD<4>::get(-F.eta[FG::ELECTRON]);

    //	 const auto eps_fraction = 0.5*(eps_tilde_m *eps_p + eps_tilde_p
    //*eps_m)/(eps_p*eps_m);
    const auto eps_fraction = 0.5 * F.temp *
                              (FDR<4, 3>::get(F.eta[FG::ELECTRON]) +
                               FDR<4, 3>::get(-F.eta[FG::ELECTRON]));

    // Pair constant in Ruffert (B8)
    const auto pair_const =
        sigma_0 * clite / int_pow<2>(me_mev) * eps_m * eps_p;
    const auto R_pair = pair_const *
                        (int_pow<2>(Cv - Ca) + int_pow<2>(Cv + Ca)) /
                        (36. * block_factor[NUE] * block_factor[NUA]);

    if (std::isfinite(R_pair)) {
      R[NUE] += R_pair;
      R[NUA] += R_pair;

      // Ruffert et al. (B16)
      Q[NUE] += R_pair * eps_fraction;
      Q[NUA] += R_pair * eps_fraction;
    }

    // Ruffert et al. (B10)
    const auto R_pair_x = 1. / 9. * pair_const *
                          (int_pow<2>(Cv - Ca) + int_pow<2>(Cv + Ca - 2.)) /
                          int_pow<2>(block_factor[NUX]);

    if (std::isfinite(R_pair_x)) {
      R[NUX] += R_pair_x;
      // Ruffert et al. (B16)
      Q[NUX] += R_pair_x * eps_fraction;
    }

#ifdef DEBUG
    std::cout << "Electrion Positron pair annihilation" << std::endl;
    std::cout << "Q : ";

    std::cout << R_pair * eps_fraction << " , ";
    std::cout << R_pair * eps_fraction << " , ";
    std::cout << R_pair_x * eps_fraction << " , ";
    std::cout << std::endl;

    std::cout << "R : ";
    std::cout << R_pair << " , ";
    std::cout << R_pair << " , ";
    std::cout << R_pair_x << " , ";

    std::cout << std::endl;
#endif

    // PLASMON DECAY!
    // Ruffert et al. see comment after (B12)
    const auto gamma =
        gamma_0 *
        sqrt(1. / 3. *
             (int_pow<2>(pi) + 3.0 * int_pow<2>(F.eta[FG::ELECTRON])));

    // Ruffert et al. (B11)
    const auto gamma_const =
        int_pow<3>(pi) * sigma_0 * clite /
        (int_pow<2>(me_mev) * 3.0 * fsc * int_pow<6>(hc_mevcm)) *
        int_pow<6>(gamma) * exp(-gamma) * (1.0 + gamma) * int_pow<8>(F.temp);

// Ruffert et al. (B13)
#pragma unroll
    for (int i = NUE; i <= NUX; ++i) {
      block_factor[i] =
          1. + exp(F.eta[i] - (1. + 0.5 * int_pow<2>(gamma) / (1. + gamma)));
    }

    // Ruffert et al. (B11)
    const auto R_gamma =
        int_pow<2>(Cv) * gamma_const / (block_factor[NUE] * block_factor[NUA]);
    const auto Q_gamma = F.temp * 0.5 * (2. + int_pow<2>(gamma) / (1. + gamma));

    if (std::isfinite(R_gamma)) {
      R[NUE] += R_gamma;
      R[NUA] += R_gamma;

      // Ruffert et al. (B17)
      Q[NUE] += Q_gamma * R_gamma;
      Q[NUA] += Q_gamma * R_gamma;
    }

    // Ruffert et al. (B12)
    const auto R_gamma_x =
        4. * int_pow<2>(Cv - 1.) * gamma_const / int_pow<2>(block_factor[NUX]);
    if (std::isfinite(R_gamma_x)) {
      R[NUX] += R_gamma_x;
      // Ruffert et al. (B17)
      Q[NUX] += Q_gamma * R_gamma_x;
    }

#ifdef DEBUG
    std::cout << "Plasmon decay" << std::endl;
    std::cout << "Q : ";

    std::cout << R_gamma * Q_gamma << " , ";
    std::cout << R_gamma * Q_gamma << " , ";
    std::cout << R_gamma_x * Q_gamma << " , ";
    std::cout << std::endl;

    std::cout << "R : ";
    std::cout << R_gamma << " , ";
    std::cout << R_gamma << " , ";
    std::cout << R_gamma_x << " , ";

    std::cout << std::endl;
#endif
    if (brems) {
      // Bremsstrahlung fitting formula described in
      // A. Burrows et al. Nuclear Physics A 777 (2006) 356-394
      // This seems to differ significantly from Sekiguchi 2011!!
      // ERM: I have checked the formula again with the paper and believe it to
      // be ok...
      const auto R_brems =
          0.231 * (2.0778e2 * erg_to_mev) * 0.5 *
          (int_pow<2>(F.Ym[FG::NEUTRON]) + int_pow<2>(F.Ym[FG::PROTON]) +
           28.0 / 3.0 * F.Ym[FG::NEUTRON] * F.Ym[FG::PROTON]) *
          int_pow<2>(F.rho) * int_pow<4>(F.temp) * sqrt(F.temp);

      const auto Q_brems = R_brems * F.temp / 0.231 * 0.504;

#pragma unroll
      for (int i = NUE; i <= NUX; ++i) {
        Q[i] += g_nu[i] * Q_brems;
        R[i] += g_nu[i] * R_brems;
      }

#ifdef DEBUG
      std::cout << "Nucleon-Nuclean Bremsstrahlung" << std::endl;
      std::cout << "Q : ";

      std::cout << Q_brems << " , ";
      std::cout << Q_brems << " , ";
      std::cout << Q_brems << " , ";
      std::cout << std::endl;

      std::cout << "R : ";
      std::cout << R_brems << " , ";
      std::cout << R_brems << " , ";
      std::cout << R_brems << " , ";

      std::cout << std::endl;

      std::cout << "Total" << std::endl;

      std::cout << "Q : ";

      std::cout << Q[0] << " , ";
      std::cout << Q[1] << " , ";
      std::cout << Q[2] << " , ";
      std::cout << std::endl;

      std::cout << "R : ";
      std::cout << R[0] << " , ";
      std::cout << R[1] << " , ";
      std::cout << R[2] << " , ";

      std::cout << std::endl;

#endif
    }
  }

  inline void calc_number_densities() {
    using namespace Leak_Constants;
    using namespace Margherita_constants;

#pragma unroll
    for (int i = NUE; i <= NUX; ++i) {
      // Ruffert et al. (B18)
      n_nu[i] = g_nu[i] * 4. * pi / int_pow<3>(hc_mevcm) * int_pow<3>(F.temp) *
                FD<2>::get(F.eta[i]);

      // Ruffert et al. (B19)
      en_nu[i] = g_nu[i] * 4. * pi / int_pow<3>(hc_mevcm) * int_pow<4>(F.temp) *
                 FD<3>::get(F.eta[i]);
    }
  }

  inline void calc_opacities() {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    using namespace Margherita_helpers;

    // The notation in the papers is utterly confusing, let's summarise the main
    // conventions  in the following:
    // Ultimately we are interested in the optical depth \tau, so really we
    // would like \tau = \int \kappa d s, but \kappa = \frac{1}{mean free path}
    // the mean free path \lambda = \frac{1}{n \sigma}: so really we would like
    // \kappa = n \sigma.
    // There is an additional simplification, we assue that we can factor out an
    // energy dependence E^2, so what we compute is in fact
    // \zeta = \frac{1}{\lambda E^2}.
    // We compute the separate \zetas from absorption and scattering, sum them
    // up and reconstruct \lambda. The big problem now is, that we don't know
    // E^2, so we take some very approximate frequency averaged value of E^2 =
    // F_4/F_2 T^2 (11) in Rosswog 2002. Although we will not make us of this,
    // another notation commonly found is \int \kappa ds = \xi E^2.

    // Scattering
    // This follows Ruffert appendix A
    // Ruffert et al (A1)
    constexpr auto Cs_n = (1. + 5. * int_pow<2>(alpha)) / 24. * sigma_0;
    constexpr auto Cs_p =
        (4. * int_pow<2>(Cv - 1.) + 5. * int_pow<2>(alpha)) / 24. * sigma_0;

    // The following equations are different from Ruffert in the sense, that we
    // do  include the presence of heavy nucleons, so that Y_n + Y_p !=1

    // Ruffert (A8)
    // Pauli blocking for degenerate matter
    auto ymp = F.Ym[FG::PROTON] / (1. + 2. / 3. * max(0., F.eta[FG::PROTON]));
    auto ymn = F.Ym[FG::NEUTRON] / (1. + 2. / 3. * max(0., F.eta[FG::NEUTRON]));

    if (!std::isnormal(ymp)) ymp = 0;
    if (!std::isnormal(ymn)) ymn = 0;

    // Ruffert (A6)
    // Note that kappa_n stores zeta for now until we multiply with E^2
    // Also we are missing a factor of (T/me c^2) but since all terms are
    // multiplied by  it, we add it in at the very end
    kappa_n[NUE] += F.nb * Cs_n * ymn;  // F.Ym[FG::NEUTRON];
    kappa_n[NUA] += F.nb * Cs_n * ymn;  // F.Ym[FG::NEUTRON];
    kappa_n[NUX] += F.nb * Cs_n * ymn;  // F.Ym[FG::NEUTRON];

    kappa_n[NUE] += F.nb * Cs_p * ymp;  // F.Ym[FG::PROTON];
    kappa_n[NUA] += F.nb * Cs_p * ymp;  // F.Ym[FG::PROTON];
    kappa_n[NUX] += F.nb * Cs_p * ymp;  // F.Ym[FG::PROTON];

    // Also include coherent neutrino nucleus scattering
    // as used by Rosswog et al 2003 (A18)
    auto C_nucleus = 1. / 16. * sigma_0 * F.Ym[FG::ABAR] *
                     int_pow<2>(1. - F.Ym[FG::ZBAR] / F.Ym[FG::ABAR]);
    if (!std::isnormal(C_nucleus)) C_nucleus = 0;
    kappa_n[NUE] += C_nucleus * F.nb * F.Ym[FG::HEAVY];
    kappa_n[NUA] += C_nucleus * F.nb * F.Ym[FG::HEAVY];
    kappa_n[NUX] += C_nucleus * F.nb * F.Ym[FG::HEAVY];

    // ABSORPTION
    constexpr auto abs_const = 0.25 * (1. + 3. * int_pow<2>(alpha)) * sigma_0;
    const auto block_factor_n =
        1. + exp(F.eta[FG::ELECTRON] - FDR<5, 4>::get(F.eta[NUE]));

    // neutrons can only absorb nue since n + nu_e -> p + e-
    // due to lepton number conservation
    const auto absorp_e = F.eta_np * abs_const / block_factor_n;
    if (std::isfinite(absorp_e)) {
      kappa_n_abs[NUE] += absorp_e;
    }

    const auto block_factor_p =
        1. + exp(-F.eta[FG::ELECTRON] - FDR<5, 4>::get(F.eta[NUA]));

    // protons can only absorb nua since p + nu_a + e- -> n
    // due to lepton number conservation
    const auto absorp_a = F.eta_pn * abs_const / block_factor_p;

    if (std::isfinite(absorp_a)) {
      kappa_n_abs[NUA] += absorp_a;
    }

    // kappa_n and kappa_e differ only by different fermi functions

    for(int i = NUE; i <= NUX; ++i){
      kappa_n[i]+= kappa_n_abs[i];
      kappa_e[i] = kappa_n[i];  // Copy zetas
      kappa_e_abs[i] = kappa_n_abs[i];
    }

// Go from \zetas to \kappas
// Also include the missing (T/me c^2)^2 factor
#pragma unroll
    for (int i = NUE; i <= NUX; ++i) {
      auto nfac = int_pow<2>(F.temp / me_mev) *
                    FDR<4, 2>::get(F.eta[i]);  // This is now in 1/cm
      kappa_n[i] *= nfac;
      kappa_n_abs[i] *= nfac;

      auto efac = int_pow<2>(F.temp / me_mev) *
                    FDR<5, 3>::get(F.eta[i]);  // This is now in 1/cm
      kappa_e[i] *= efac;
      kappa_e_abs[i] *= efac;

      if (!std::isfinite(kappa_n[i]) || kappa_n[i] < 0.) kappa_n[i] = 0;
      if (!std::isfinite(kappa_e[i]) || kappa_e[i] < 0.) kappa_e[i] = 0;

      if (!std::isfinite(kappa_n_abs[i]) || kappa_n_abs[i] < 0.) kappa_n_abs[i] = 0;
      if (!std::isfinite(kappa_e_abs[i]) || kappa_e_abs[i] < 0.) kappa_e_abs[i] = 0;
    }
#ifdef DEBUG

    std::cout << F.eta_np << std::endl;

    std::cout << "Absorption" << std::endl;
    std::cout << "kappa_n : ";
    for (auto &q : kappa_n) std::cout << q << " , ";
    std::cout << std::endl;
    std::cout << "kappa_e : ";
    for (auto &q : kappa_e) std::cout << q << " , ";
    std::cout << std::endl;
#endif
  }

  template <typename tau_type>
  inline void calc_diffusion_rates_Ruffert(tau_type &__tau_n,
                                           tau_type &__tau_e) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;

    // Diffusion rates according to Ruffert 1996 (as used in THC_LeakageBase)

    calc_number_densities();

    for (int i = NUE; i <= NUX; ++i) {
      // Ruffert et al. 1996 (B20) and (B22)
      //	R_diff[i] =  g_nu[i] * clite * n_nu[i] *
      //kappa_n[i]/int_pow<2>(__tau_n[i])/(a_diff);
      R_diff_inv[i] = a_diff /
                      (g_nu[i] * clite * normfact * n_nu[i] * kappa_n[i]) *
                      int_pow<2>(__tau_n[i]);

      // Ruffert et al. 1996 (B21) and (B23)
      // Q_diff[i] =  g_nu[i] * clite * en_nu[i] *
      // kappa_e[i]/int_pow<2>(__tau_e[i])/(a_diff);
      Q_diff_inv[i] = a_diff /
                      (g_nu[i] * normfact * clite * en_nu[i] * kappa_e[i]) *
                      int_pow<2>(__tau_e[i]);
    }
  }

  template <typename tau_type>
  inline void calc_diffusion_rates(tau_type &__tau_n, tau_type &__tau_e) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;

    constexpr auto factor =
        4. / a_diff * pi * clite / int_pow<3>(hc_mevcm) * normfact;

    // Diffusion rates according to Rosswog et al 2002.
    // See also Sekiguchi 2012
    for (int i = NUE; i <= NUX; ++i) {
      const auto E_sq = int_pow<2>(F.temp) * FDR<4, 2>::get(F.eta[i]);  // In
                                                                        // MeV

      //	R_diff[i] = factor * g_nu[i] *
      //(kappa_n[i]/int_pow<2>(__tau_n[i]))
      //	            * F.temp * FD<0>::get(F.eta[i]) * E_sq;
      R_diff_inv[i] =
          int_pow<2>(__tau_n[i]) / (factor * g_nu[i] * kappa_n[i] * F.temp *
                                    FD<0>::get(F.eta[i]) * E_sq);

      const auto E_sq_e = int_pow<2>(F.temp) * FDR<5, 3>::get(F.eta[i]);

      // Q_diff[i] = factor * g_nu[i] * (kappa_e[i]/int_pow<2>(__tau_e[i]))
      //  	    * int_pow<2>(F.temp) * FD<1>::get(F.eta[i]) * E_sq_e;
      Q_diff_inv[i] = int_pow<2>(__tau_e[i]) /
                      (factor * g_nu[i] * kappa_e[i] * int_pow<2>(F.temp) *
                       FD<1>::get(F.eta[i]) * E_sq_e);
    }
  }

  inline void calc_diffusion_rates() {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    calc_diffusion_rates<>(F.tau_n, F.tau_e);
  }

  inline void calc_diffusion_rates_Ruffert() {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    calc_diffusion_rates_Ruffert<>(F.tau_n, F.tau_e);
  }

  inline void calc_effective_rates() {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    for (int i = NUE; i <= NUX; ++i) {
      // Catch floating point overflow or division by zero
      if (std::isfinite(R[i] * R_diff_inv[i]) *
          std::isfinite(10. * R_diff_inv[i])) {
        R_eff[i] = R[i] / (1. + normfact * R[i] * R_diff_inv[i]);
        // Catch floating point underflow
        if (!std::isnormal(R_diff_inv[i])) {
          R_eff[i] = R[i];
        }
      } else {
        R_eff[i] = 0;
      }

      // Catch floating point overflow or division by zero
      if (std::isfinite(Q[i] * Q_diff_inv[i]) *
          std::isfinite(10. * Q_diff_inv[i])) {
        Q_eff[i] = Q[i] / (1. + normfact * Q[i] * Q_diff_inv[i]);

        // Catch floating point underflow
        if (!std::isnormal(Q_diff_inv[i])) {
          Q_eff[i] = Q[i];
        }
      } else {
        Q_eff[i] = 0;
      }
    }
  }

  inline void calc_leakage() {
    calc_emission();
    calc_opacities();
    calc_diffusion_rates();
    calc_effective_rates();
  }

  inline void calc_leakage_Ruffert() {
    calc_emission();
    calc_opacities();
    calc_diffusion_rates_Ruffert();
    calc_effective_rates();
  }

  // CONSTRUCTORS //

  Leakage(Fugacities &&_F) : F(_F){};

  Leakage(const double &rho, const double &ye, const double &temp)
      : F(Fugacities(rho, temp, ye)){};

  // FIXME: These constructors don't compile yet...
  template <typename F>
  Leakage(const double &rho, const double &ye, const double &temp, F &&tau_n)
      : F(Fugacities(rho, temp, ye, std::forward<F>(tau_n))){};

  template <typename F>
  Leakage(const double &rho, const double &ye, const double &temp, F &&tau_n,
          F &&tau_e)
      : F(Fugacities(rho, temp, ye, std::forward<F>(tau_n),
                     std::forward<F>(tau_e))){};

  // Access operators

  // ENERGY RATES //
  inline double Q_cgs(int i) const {
    using namespace Leak_Constants;
    return Q_eff[i] * mev_to_erg;
  }

  inline double Q_cactus(int i) const {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return Q_eff[i] * mev_to_erg * RHOGF * EPSGF / TIMEGF;
  }

  inline double Q_emission_cgs(int i) const {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return Q[i] * mev_to_erg;
  }

  inline double Q_emission_cactus(int i) const {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return Q[i] * mev_to_erg * RHOGF * EPSGF / TIMEGF;
  }

  inline double Q_diff_cgs(int i) const {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return 1. / (normfact * Q_diff_inv[i]) * mev_to_erg;
  }

  inline double Q_diff_cactus(int i) const {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return 1. / (normfact * Q_diff_inv[i]) * mev_to_erg * RHOGF * EPSGF /
           TIMEGF;
  }

  // NUMBER RATES //
  inline double R_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return R_eff[i] * massn_cgs;
  }

  inline double R_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return R_eff[i] * massn_cgs * RHOGF / TIMEGF;
  }

  inline double R_emission_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return R[i] * massn_cgs;
  }

  inline double R_emission_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return R[i] * massn_cgs * RHOGF / TIMEGF;
  }

  inline double R_diff_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return 1. / (normfact * R_diff_inv[i]) * massn_cgs;
  }

  inline double R_diff_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return 1. / (normfact * R_diff_inv[i]) * massn_cgs * RHOGF / TIMEGF;
  }

  // OPACITES //
  inline double kappa_n_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return kappa_n[i];
  }

  inline double kappa_n_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return kappa_n[i] / LENGTHGF;
  }

  inline double kappa_e_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return kappa_e[i];
  }

  inline double kappa_e_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return kappa_e[i] / LENGTHGF;
  }

  inline double kappa_n_abs_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return kappa_n_abs[i];
  }

  inline double kappa_n_abs_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return kappa_n_abs[i] / LENGTHGF;
  }

  inline double kappa_e_abs_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return kappa_e_abs[i];
  }

  inline double kappa_e_abs_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return kappa_e_abs[i] / LENGTHGF;
  }

  // NUMBER DENSITIES //
  inline double n_nu_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return n_nu[i];
  }

  inline double n_nu_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return n_nu[i] / (int_pow<3>(LENGTHGF));
  }

  inline double en_nu_cgs(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return en_nu[i] * mev_to_erg;
  }

  inline double en_nu_cactus(int i) {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    return en_nu[i] * mev_to_erg * RHOGF * EPSGF;
  }
};

inline double compute_analytic_opacity(const double rho) {
  using namespace Margherita_constants;

  // This computes an empirical value of the opacity for cold neutron stars
  // as given in Deaton et al. 2013. (10)

  const auto rho_cgs = rho * INVRHOGF;
  const auto log_tau_approx = 0.96 * (log(rho_cgs) / log(10.) - 11.7);
  return exp(log(10.) * log_tau_approx);
}

#endif
