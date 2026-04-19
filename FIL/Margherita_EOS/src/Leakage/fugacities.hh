/*
 * =====================================================================================
 *
 *       Filename:  fugacities.hh
 *
 *    Description:  Class to store neutrino fugacities
 *
 *        Version:  1.0
 *        Created:  13/05/2017 21:14:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#include <cmath>
#include <limits>
#include "../margherita.hh"
#include "leak_constants.hh"

#ifndef __HH_FUGACITY
#define __HH_FUGACITY

class Fugacities {
 public:
  enum { NUE = 0, NUA, NUX, ELECTRON, PROTON, NEUTRON, HAT, NUM_INDICES };
  enum {
    ALPHA = 0,
    HEAVY,
    ABAR,
    ZBAR,
    MNUM_INDICES = ZBAR + 1 + 2
  };  // ALSO account for PROTON and NEUTRON
  std::array<double, NUM_INDICES> eta;
  std::array<double, MNUM_INDICES> Ym;

  double eta_np;
  double eta_pn;
  double nb;

  double rho;
  double temp;
  double ye;

  std::array<double, 3> tau_n{{std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::quiet_NaN()}};
  std::array<double, 3> tau_e;

  inline void compute_fugacities() {
    using namespace Leak_Constants;
    using namespace Margherita_constants;
    using namespace Margherita_helpers;

    // Making this call will automatically enforce table bounds
    typename EOS_Tabulated::error_type error;
    eta[ELECTRON] =
        EOS_Tabulated::mue_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_ye(
            eta[PROTON], eta[NEUTRON], Ym[ALPHA], Ym[HEAVY], Ym[NEUTRON],
            Ym[PROTON], Ym[ABAR], Ym[ZBAR], temp, rho, ye, error);

    // We need to convert the density to cgs now
    rho *= INVRHOGF * (mnuc_cgs / EOS_Tabulated::baryon_mass);

    // Eta still stores chemical potentials, we now fill the remaining gaps
    // and then normalize by T
    eta[HAT] = eta[NEUTRON] - eta[PROTON] - Qnp;
    eta[NUE] = eta[ELECTRON] + eta[PROTON] - eta[NEUTRON];
    eta[NUA] = -eta[NUE];
    eta[NUX] = 0;

    eta[PROTON] += Qnp;

    // Interpolate between optically thin and thick limit

    const auto eta_e_fac = (1. - exp(-tau_n[NUE]));
    if (std::isfinite(eta_e_fac)) {
      eta[NUE] *= eta_e_fac;
    }
    const auto eta_a_fac = (1. - exp(-tau_n[NUA]));
    if (std::isfinite(eta_a_fac)) {
      eta[NUA] *= eta_a_fac;
    }

    // Now we obtain the fugacities
    for (auto &etaL : eta) etaL /= temp;

    // Enforce positivity of mass fractions
    for (auto &YL : Ym) YL = max(0, YL);

    nb = rho * avogadro;  // Neutrino number density

    // nucleon final state blocking, see Rosswog
    // Ruffert et al. (A13)
    eta_np = nb * (Ym[PROTON] - Ym[NEUTRON]) / (exp(-eta[HAT]) - 1.);
    eta_pn = nb * (Ym[NEUTRON] - Ym[PROTON]) / (exp(eta[HAT]) - 1.);
    // eta_np = exp(eta[HAT])*eta_pn;

    // Bruenn (ApJSS 58 1985) (3.1): non degenerate matter limit.
    if (rho < 2.e11) {
      eta_pn = nb * Ym[PROTON];
      eta_np = nb * Ym[NEUTRON];
    }

    if (!std::isnormal(eta_np)) eta_np = 0.;
    if (!std::isnormal(eta_pn)) eta_pn = 0.;

    // Rosswog (A9)
    eta_pn = std::max(eta_pn, 0.);
    eta_np = std::max(eta_np, 0.);
  }

  Fugacities(const double &rho, const double &temp, const double &ye)
      : rho(rho), temp(temp), ye(ye) {
    compute_fugacities();
  }

  template <typename F>
  Fugacities(const double &rho, const double &temp, const double &ye,
             F &&__tau_n)
      : rho(rho), temp(temp), ye(ye), tau_n(std::forward<F>(__tau_n)) {
    compute_fugacities();
  }

  template <typename F>
  Fugacities(const double &rho, const double &temp, const double &ye, F &&tau_n,
             F &&tau_e)
      : rho(rho),
        temp(temp),
        ye(ye),
        tau_n(std::forward<F>(tau_n)),
        tau_e(std::forward<F>(tau_e)) {
    compute_fugacities();
  }
};
#endif
