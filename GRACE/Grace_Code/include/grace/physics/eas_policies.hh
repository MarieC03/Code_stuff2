/**
 * @file eas_policies.hh
 * @author Carlo Musolino (carlo.musolino@aei.mpg.de)
 * @brief EAS source-term operators (tests, photons, neutrinos)
 * @date 2024-05-13
 *
 * @copyright This file is part of of the General Relativistic Astrophysics
 * Code for Exascale.
 * GRACE is an evolution framework that uses Finite Volume
 * methods to simulate relativistic spacetimes and plasmas
 * Copyright (C) 2023 Carlo Musolino
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef GRACE_PHYSICS_EAS_POLICIES_HH
#define GRACE_PHYSICS_EAS_POLICIES_HH

#include <grace_config.h>

#include <grace/physics/m1_helpers.hh>
#include <grace/physics/m1.hh>

#include <grace/utils/device.h>
#include <grace/utils/inline.h>

#include <grace/physics/eos/eos_base.hh>
#include <grace/physics/eos/eos_storage.hh>
#include <grace/physics/eos/physical_constants.hh>

#include <grace/physics/eas_neutrino_rates_analytic.hh>

#include <grace/config/config_parser.hh>
#include <grace/errors/assert.hh>

#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

#include <string>
#include <type_traits>

namespace grace {

inline std::string resolve_m1_eas_kind_host()
{
    std::string kind = "test";
    try {
        kind = grace::get_param<std::string>("m1", "eas", "kind");
    } catch(...) {}

    bool use_analytic = false;
    bool use_weakhub = false;
    try { use_analytic = grace::get_param<bool>("m1", "eas", "use_analytic"); } catch(...) {}
    try { use_weakhub = grace::get_param<bool>("m1", "eas", "use_weakhub"); } catch(...) {}

    if (use_analytic) return "neutrino_analytic";
    if (use_weakhub) return "neutrino_weakhub";
    return kind;
}

template <typename ViewT>
KOKKOS_INLINE_FUNCTION inline void set_m1_species_rates(
    ViewT const& view,
    int const ispec,
    double const kappa_a,
    double const kappa_s,
    double const eta_e,
    double const eta_n,
    double const kappa_n)
{
    view(m1_aux_index(KAPPAA_, ispec)) = kappa_a;
    view(m1_aux_index(KAPPAS_, ispec)) = kappa_s;
    view(m1_aux_index(ETA_, ispec)) = eta_e;
    view(m1_aux_index(ETAN_, ispec)) = eta_n;
    view(m1_aux_index(KAPPAAN_, ispec)) = kappa_n;
}

template <typename ViewT>
KOKKOS_INLINE_FUNCTION inline void fill_m1_species_rates(
    ViewT const& view,
    double const kappa_a,
    double const kappa_s,
    double const eta_e,
    double const eta_n,
    double const kappa_n)
{
    for (int ispec = 0; ispec < m1_num_species(); ++ispec) {
        set_m1_species_rates(view, ispec, kappa_a, kappa_s, eta_e, eta_n, kappa_n);
    }
}

GRACE_HOST_DEVICE inline constexpr int active_m1_species_to_physical_species(int ispec)
{
#ifdef M1_NU_FIVESPECIES
    return ispec;
#elif defined(M1_NU_THREESPECIES)
    return (ispec == 0) ? NUE : ((ispec == 1) ? NUEBAR : NUX);
#else
    return NUE;
#endif
}

//------------------------------------------------------------------------------
// Test EAS operator
//------------------------------------------------------------------------------
struct test_eas_op {
    enum test_t  {
        ZERO_EAS=0,
        LARGE_KS,
        EMITTING_SPHERE,
        SHADOW_CAST,
        COUPLING_TEST
    } ;
    test_eas_op(
        grace::var_array_t _aux
    ) : aux(_aux)
    {
        auto _which_test = grace::get_param<std::string>(
            "m1", "id_type"
        ) ;
        if (_which_test == "straight_beam" or
            _which_test == "curved_beam" )
        {
            which_test = ZERO_EAS ;
        } else if (
            _which_test == "scattering"
            or _which_test == "moving_scattering"
        ) {
            which_test = LARGE_KS ;
            _ks_value = grace::get_param<double>("m1","scattering_test","k_s") ;
        } else if (
            _which_test == "shadow"
        ) {
            which_test = SHADOW_CAST;
        } else if (
            _which_test == "emitting_sphere"
        ) {
            which_test = EMITTING_SPHERE ;
            _emitting_sphere_temperature = grace::get_param<double>("m1","emitting_sphere_test","temperature") ;
            _emitting_sphere_cross_section = grace::get_param<double>("m1","emitting_sphere_test","cross_section") ;
        } else if ( _which_test == "coupling_test") {
            which_test = COUPLING_TEST ;
        } else {
            ERROR("Unknown m1 test") ;
        }
    }

    void KOKKOS_INLINE_FUNCTION
    operator() (
        VEC(const int i, const int j, const int k), int64_t q
        , double* xyz
    ) const
    {
        auto u = Kokkos::subview(aux,VEC(i,j,k),Kokkos::ALL(),q) ;
        double r=0;
        switch (which_test) {
            case ZERO_EAS:
            fill_m1_species_rates(u, 0.0, 0.0, 0.0, 0.0, 0.0);
            break ;
            case LARGE_KS:
            fill_m1_species_rates(u, 0.0, _ks_value, 0.0, 0.0, 0.0);
            break ;
            case SHADOW_CAST:
            // we assume pcoords is cartesian
            r = sqrt(
                SQR(xyz[0]) + SQR(xyz[1]) + SQR(xyz[2])
            ) ;
            fill_m1_species_rates(u, 0.0, 0.0, 0.0, 0.0, 0.0);
            if ( r<0.046875 ) {
                fill_m1_species_rates(u, 1e06, 0.0, 0.0, 0.0, 1e06);
            }
            break ;
            case EMITTING_SPHERE:
            // we assume pcoords is cartesian
            r = sqrt(
                SQR(xyz[0]) + SQR(xyz[1]) + SQR(xyz[2])
            ) ;
            fill_m1_species_rates(u, 0.0, 0.0, 0.0, 0.0, 0.0);
            if ( r < 1. ) {
                double T = _emitting_sphere_temperature ;
                double sigma = _emitting_sphere_cross_section ;
                // we set the rates according to LTE,
                // for simplicity the Stefan Boltzmann constant is 1 here
                fill_m1_species_rates(
                    u,
                    sigma,
                    0.0,
                    SQR(T) * SQR(T) * sigma,
                    SQR(T) * T * sigma,
                    sigma);
            }
            break;
            case COUPLING_TEST:
            r = sqrt(
                SQR(xyz[0]) + SQR(xyz[1]) + SQR(xyz[2])
            ) ;
            fill_m1_species_rates(u, 0.0, 0.0, 0.0, 0.0, 0.0);
            if ( r < 1.0 ) {
                if ( r < 0.5 ) {
                    // effectively T = 1
                    fill_m1_species_rates(u, 1.0, 0.0, 0.01, 0.01, 1.0);
                } else {
                    double T = 1. - (r-0.5)/0.5 ;
                    fill_m1_species_rates(u, 1.0, 0.0, 0.01 * T * T * T * T, 0.01 * T * T * T, 1.0);
                }
            }
            break ;
        }
    }
    var_array_t aux ;
    test_t which_test;
    double _ks_value ;
    double _emitting_sphere_cross_section, _emitting_sphere_temperature;
} ;

//------------------------------------------------------------------------------
// Photon EAS operator
//------------------------------------------------------------------------------
struct photon_eas_op {
  explicit photon_eas_op(var_array_t _aux)
      : mass_scale(grace::get_param<double>("coordinate_system", "mass_scale")), aux(_aux) {}

  void KOKKOS_INLINE_FUNCTION
  operator()(VEC(const int i, const int j, const int k), int64_t q, double* xyz)
  const {
      using namespace grace::physical_constants ;
      double mu = 0.5 ; // fully ionized, fixme
      double const rho = aux(VEC(i,j,k),RHO_,q)  * Msun_cgs / (Msun_to_cm*SQR(Msun_to_cm*mass_scale)) ;
      double const T   = k_cgs * aux(VEC(i,j,k),TEMP_,q) / (mu * mp_cgs);// erg
      // Energy rates
      double eta_cgs = 1.4e-27 * sqrt(T/k_cgs) * SQR(rho)/me_cgs/mp_cgs ; // gaunt factor 1
      // Kirchoff law
      // 2 *( k_cgs *T )**4 * np.pi**4 / ( 15 * clight**2 * h_cgs**3)
      double BB = 2 * SQR(M_PI*T)*SQR(M_PI*T) / ( 15. * SQR(clight) * SQR(h_cgs) * h_cgs) ;
      // Planck mean opacity
      double kappa_cgs = eta_cgs / BB ;
      double BBn = 4 * SQR(T)*T / (SQR(clight)*SQR(h_cgs)*h_cgs) * 1.20206 ; // numerical factor is Zeta(3)
      double const eta_e = eta_cgs  / Msun_to_erg * SQR(Msun_to_cm)*Msun_to_cm * Msun_to_s * SQR(mass_scale) * mass_scale ;
      double const kappa_a = kappa_cgs * Msun_to_cm * mass_scale ;
      // Number rates
      // The integral over frequency is IR divergence, we cutoff at the plasma frequency
      // and plop the result (~log(h nu_min/(kT))) into the Gaunt factor
      // now -- we don't really care about eta anyway, eta/kappa is unaffected by this
      // for now we prentend the gaunt factor stays one and sweep this under the rug
      // but FIXME please
      double const eta_n = eta_cgs / T * SQR(Msun_to_cm)*Msun_to_cm * Msun_to_s * SQR(mass_scale) * SQR(mass_scale) / Msun_to_erg ;
      double const kappa_n = eta_cgs / T / BBn * Msun_to_cm * mass_scale ;
      // Scattering
      auto u = Kokkos::subview(aux,VEC(i,j,k),Kokkos::ALL(),q) ;
      fill_m1_species_rates(u, kappa_a, 0.0, eta_e, eta_n, kappa_n);
  }

  var_array_t aux;
  double mass_scale{1.0};
};

//------------------------------------------------------------------------------
// Neutrino EAS operator
//------------------------------------------------------------------------------
template <typename eos_t>
struct neutrinos_eas_op {
  using tau_kind_t = int;
  enum : tau_kind_t {
    TAU_NONE = 0,
    TAU_ANALYTIC_DENSITY = 1,
    TAU_LOCAL_SPHERICAL = 2
  };
  enum : int {
    BETAEQ_OFF = 0,
    BETAEQ_CHEMICAL = 1
  };

  explicit neutrinos_eas_op(var_array_t _aux)
      : eos(eos::get().get_eos<eos_t>()),
        aux(_aux),
        mass_scale(grace::get_param<double>("coordinate_system", "mass_scale")),
        beta_decay(grace::get_param<bool>("m1", "eas", "beta_decay")),
        plasmon_decay(grace::get_param<bool>("m1", "eas", "plasmon_decay")),
        bremsstrahlung(grace::get_param<bool>("m1", "eas", "bremsstrahlung")),
        pair_annihilation(grace::get_param<bool>("m1", "eas", "pair_annihilation")),
        apply_temp_correction(grace::get_param<bool>("m1", "eas", "temperature_correction")),
        use_weakhub(parse_use_weakhub_host()),
        betaeq_mode(parse_betaeq_mode_host()),
        tau_kind(parse_tau_kind_host()),
        // betaeq_mode(BETAEQ_OFF),
        // tau_kind(TAU_NONE),
        weakhub(grace::weakhub::get_device_handle())
  {
    validate_configuration_host();
    // Host-only parsing: strings are not device-friendly.
    // TODO: parser optimize
    //betaeq_mode = parse_betaeq_mode_host();
    //tau_kind = parse_tau_kind_host();
    spherical_tau.r_outer_code = grace::get_param<double>("m1", "eas", "tau_outer_radius_code");
    spherical_tau.seed_with_analytic = true;
  }

  int parse_betaeq_mode_host() const {
    try {
      const auto mode = grace::get_param<std::string>("m1", "eas", "betaeq_policy");
      if (mode == "chemical") return BETAEQ_CHEMICAL;
    } catch(...) {}
    return BETAEQ_OFF;
  }

  tau_kind_t parse_tau_kind_host() const {
    try {
      const auto mode = grace::get_param<std::string>("m1", "eas", "tau_policy");
      if (mode == "none") return TAU_NONE;
      if (mode == "analytic_density") return TAU_ANALYTIC_DENSITY;
      if (mode == "local_spherical") return TAU_LOCAL_SPHERICAL;
    } catch(...) {}
    return TAU_NONE;
  }

  bool parse_use_weakhub_host() const {
    return grace::weakhub::weakhub_enabled_from_params();
  }

  void validate_configuration_host() const {
    const auto eas_kind = resolve_m1_eas_kind_host();
    if (eas_kind != "neutrino_analytic" && eas_kind != "neutrino_weakhub") {
      return;
    }

    const auto species_mode = grace::get_param<std::string>("m1", "eas", "species_mode");
    const bool use_analytic = grace::get_param<bool>("m1", "eas", "use_analytic");
    const bool use_weakhub_cfg = grace::get_param<bool>("m1", "eas", "use_weakhub");
    const auto eos_type = grace::get_param<std::string>("eos", "eos_type");

    ASSERT(eos_type == "tabulated" || eos_type == "leptonic_4d",
           "Microphysical neutrino M1 rates require a Ye-dependent tabulated EOS "
           "or the leptonic_4d EOS.");

    if (eas_kind == "neutrino_analytic") {
      ASSERT(!use_weakhub_cfg,
             "m1.eas.kind is neutrino_analytic, but m1.eas.use_weakhub is true. "
             "Choose one rates treatment.");
    } else if (eas_kind == "neutrino_weakhub") {
      ASSERT(!use_analytic,
             "m1.eas.kind is neutrino_weakhub, but m1.eas.use_analytic is true. "
             "Choose one rates treatment.");
    }

    if (eos_type == "leptonic_4d") {
      const bool use_muonic_eos =
          grace::get_param<bool>("eos", "leptonic_4d", "use_muonic_eos");
#ifdef M1_NU_FIVESPECIES
      ASSERT(use_muonic_eos,
             "The 5-species M1 build requires eos.leptonic_4d.use_muonic_eos = true "
             "so muonic matter couples consistently to nu_mu and anti-nu_mu.");
#else
      (void)use_muonic_eos;
#endif
    }

#ifdef M1_NU_FIVESPECIES
    ASSERT(species_mode == "build_default" || species_mode == "five",
           "This GRACE executable was built with five M1 neutrino species. "
           "Set m1.eas.species_mode to \"five\" or \"build_default\".");
    ASSERT(eas_kind == "neutrino_weakhub",
           "The 5-species muonic transport path is only supported with Weakhub "
           "rates, matching the FIL implementation.");
    ASSERT(eos_type == "leptonic_4d",
           "The 5-species muonic transport path requires eos.eos_type = \"leptonic_4d\".");
#elif defined(M1_NU_THREESPECIES)
    ASSERT(species_mode == "build_default" || species_mode == "three",
           "This GRACE executable was built with three M1 neutrino species. "
           "Set m1.eas.species_mode to \"three\" or \"build_default\".");
#else
    ASSERT(species_mode == "build_default" || species_mode == "one",
           "This GRACE executable was built with one M1 neutrino species. "
           "Set m1.eas.species_mode to \"one\" or \"build_default\".");
#endif
  }

  // --- beta equilibrium: mu_e + mu_p - mu_n - Qnp = 0 (chemical equilibrium) ---
  GRACE_HOST_DEVICE inline double beta_eq_residual(double rho, double T, double Ye) const {
    double mu_p=0.0, mu_n=0.0;
    double Xa=0.0, Xh=0.0, Xn=0.0, Xp=0.0, Abar=1.0, Zbar=1.0;
    typename eos_t::error_type err = 0;
    double ye_loc = Ye, T_loc = T, rho_loc = rho;
    const double mu_e = eos.mue_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_ye(mu_p, mu_n, Xa, Xh, Xn, Xp, Abar, Zbar, T_loc, rho_loc, ye_loc, err);
    if (err.any() || !::isfinite(mu_e) || !::isfinite(mu_p) || !::isfinite(mu_n)) return 0.0;
    return (mu_e + mu_p - mu_n - nu_constants::Qnp);
  }

  GRACE_HOST_DEVICE inline double find_ye_betaeq(double rho, double T, double Ye0) const {
    // Device-friendly bisection. If no bracket -> return Ye0.
    double a = 1.0e-6, b = 0.60;
    double fa = beta_eq_residual(rho, T, a), fb = beta_eq_residual(rho, T, b);
    if (!::isfinite(fa) || !::isfinite(fb) || fa * fb > 0.0) return Ye0;

    double left = a, right = b, fleft = fa, mid = Ye0;
    for (int it = 0; it < 40; ++it) {
      mid = 0.5 * (left + right);
      const double fm = beta_eq_residual(rho, T, mid);
      if (!::isfinite(fm) || fm == 0.0) break;
      if (fleft * fm <= 0.0) right = mid; else { left = mid; fleft = fm; }
      if ((right - left) < 1.0e-8) break;
    }
    return mid;
  }

#ifdef GRACE_ENABLE_LEPTONIC_4D
  template <typename eos_u = eos_t,
            typename std::enable_if_t<std::is_same_v<eos_u, leptonic_eos_4d_t>, int> = 0>
  GRACE_HOST_DEVICE inline double beta_eq_residual_4d(double rho, double T, double Ye, double Ymu) const {
    typename eos_t::error_type err = 0;
    double rho_loc = rho, temp_loc = T, ye_loc = Ye, ymu_loc = Ymu;
    double mu_mu = 0.0, mu_p = 0.0, mu_n = 0.0;
    const double mu_e = eos.mue_mumu_mup_mun__temp_rho_ye_ymu(
        mu_mu, mu_p, mu_n, temp_loc, rho_loc, ye_loc, ymu_loc, err);
    if (err.any() || !::isfinite(mu_e) || !::isfinite(mu_p) || !::isfinite(mu_n)) return 0.0;
    return mu_e + mu_p - mu_n - nu_constants::Qnp;
  }

  template <typename eos_u = eos_t,
            typename std::enable_if_t<std::is_same_v<eos_u, leptonic_eos_4d_t>, int> = 0>
  GRACE_HOST_DEVICE inline double mu_eq_residual_4d(double rho, double T, double Ye, double Ymu) const {
    typename eos_t::error_type err = 0;
    double rho_loc = rho, temp_loc = T, ye_loc = Ye, ymu_loc = Ymu;
    double mu_mu = 0.0, mu_p = 0.0, mu_n = 0.0;
    const double mu_e = eos.mue_mumu_mup_mun__temp_rho_ye_ymu(
        mu_mu, mu_p, mu_n, temp_loc, rho_loc, ye_loc, ymu_loc, err);
    if (err.any() || !::isfinite(mu_e) || !::isfinite(mu_mu)) return 0.0;
    return mu_mu - mu_e;
  }

  template <typename eos_u = eos_t,
            typename std::enable_if_t<std::is_same_v<eos_u, leptonic_eos_4d_t>, int> = 0>
  GRACE_HOST_DEVICE inline double solve_ye_for_ymu_betaeq_4d(
      double rho, double T, double Ymu, double Ye0) const
  {
    const double a = eos.get_c2p_ye_min();
    const double b = eos.get_c2p_ye_max();
    const double fa = beta_eq_residual_4d(rho, T, a, Ymu);
    const double fb = beta_eq_residual_4d(rho, T, b, Ymu);
    if (!::isfinite(fa) || !::isfinite(fb)) return Ye0;
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if (fa * fb > 0.0) return (Kokkos::fabs(fa) < Kokkos::fabs(fb)) ? a : b;

    double left = a, right = b, fleft = fa, mid = Ye0;
    for (int it = 0; it < 48; ++it) {
      mid = 0.5 * (left + right);
      const double fm = beta_eq_residual_4d(rho, T, mid, Ymu);
      if (!::isfinite(fm) || fm == 0.0) break;
      if (fleft * fm <= 0.0) right = mid; else { left = mid; fleft = fm; }
      if ((right - left) < 1.0e-10) break;
    }
    return mid;
  }
#endif

  GRACE_HOST_DEVICE inline void solve_betaeq_state(double rho, double T, double& Ye, double& Ymu) const {
#if defined(GRACE_ENABLE_LEPTONIC_4D) && defined(M1_NU_FIVESPECIES)
    if constexpr (std::is_same_v<eos_t, leptonic_eos_4d_t>) {
      const double a = eos.get_c2p_ymu_min();
      const double b = eos.get_c2p_ymu_max();
      const double fa = mu_eq_residual_4d(rho, T, solve_ye_for_ymu_betaeq_4d(rho, T, a, Ye), a);
      const double fb = mu_eq_residual_4d(rho, T, solve_ye_for_ymu_betaeq_4d(rho, T, b, Ye), b);
      double ymu_mid = Ymu;

      if (!::isfinite(fa) || !::isfinite(fb)) {
        Ye = solve_ye_for_ymu_betaeq_4d(rho, T, Ymu, Ye);
        return;
      }
      if (fa == 0.0) {
        ymu_mid = a;
      } else if (fb == 0.0) {
        ymu_mid = b;
      } else if (fa * fb > 0.0) {
        ymu_mid = (Kokkos::fabs(fa) < Kokkos::fabs(fb)) ? a : b;
      } else {
        double left = a, right = b, fleft = fa;
        for (int it = 0; it < 48; ++it) {
          ymu_mid = 0.5 * (left + right);
          const double ye_mid = solve_ye_for_ymu_betaeq_4d(rho, T, ymu_mid, Ye);
          const double fm = mu_eq_residual_4d(rho, T, ye_mid, ymu_mid);
          if (!::isfinite(fm) || fm == 0.0) break;
          if (fleft * fm <= 0.0) right = ymu_mid; else { left = ymu_mid; fleft = fm; }
          if ((right - left) < 1.0e-10) break;
        }
      }

      Ymu = ymu_mid;
      Ye = solve_ye_for_ymu_betaeq_4d(rho, T, Ymu, Ye);
      return;
    }
#endif
    Ye = find_ye_betaeq(rho, T, Ye);
  }

   // Main Kernel
  void KOKKOS_INLINE_FUNCTION operator()(VEC(const int i, const int j, const int k), int64_t q, double* xyz) const {
    const double rho = aux(VEC(i,j,k),RHO_,q);
    const double T   = aux(VEC(i,j,k),TEMP_,q);
    double Ye        = aux(VEC(i,j,k),YE_,q);
    double Ymu       = 0.0;
#ifdef GRACE_ENABLE_LEPTONIC_4D
    if constexpr (m1_explicit_muon_transport_enabled()) {
      Ymu = aux(VEC(i,j,k),YMU_,q);
    }
#endif
    if (betaeq_mode == BETAEQ_CHEMICAL) solve_betaeq_state(rho, T, Ye, Ymu);

    tau_policy_none tau_none{};
    tau_policy_analytic_density tau_ana{};
    nu_rates_all_out all{};

    if (use_weakhub && weakhub.valid) {
      switch (tau_kind) {
        case TAU_LOCAL_SPHERICAL:
          all = compute_all_species_weakhub(eos, weakhub, rho, T, Ye, Ymu, mass_scale, plasmon_decay, bremsstrahlung, pair_annihilation, xyz, spherical_tau, apply_temp_correction);
          break;
        case TAU_ANALYTIC_DENSITY:
          all = compute_all_species_weakhub(eos, weakhub, rho, T, Ye, Ymu, mass_scale, plasmon_decay, bremsstrahlung, pair_annihilation, xyz, tau_ana, apply_temp_correction);
          break;
        default:
          all = compute_all_species_weakhub(eos, weakhub, rho, T, Ye, Ymu, mass_scale, plasmon_decay, bremsstrahlung, pair_annihilation, xyz, tau_none, apply_temp_correction);
          break;
      }
    } else {
      switch (tau_kind) {
        case TAU_LOCAL_SPHERICAL:
          all = compute_all_species(eos, rho, T, Ye, Ymu, mass_scale, beta_decay, plasmon_decay, bremsstrahlung, pair_annihilation, xyz, spherical_tau, apply_temp_correction);
          break;
        case TAU_ANALYTIC_DENSITY:
          all = compute_all_species(eos, rho, T, Ye, Ymu, mass_scale, beta_decay, plasmon_decay, bremsstrahlung, pair_annihilation, xyz, tau_ana, apply_temp_correction);
          break;
        default:
          all = compute_all_species(eos, rho, T, Ye, Ymu, mass_scale, beta_decay, plasmon_decay, bremsstrahlung, pair_annihilation, xyz, tau_none, apply_temp_correction);
          break;
      }
    }

    auto u = Kokkos::subview(aux,VEC(i,j,k),Kokkos::ALL(),q) ;
    for (int ispec = 0; ispec < m1_num_species(); ++ispec) {
      const nu_rates_out r = all.out[active_m1_species_to_physical_species(ispec)];
      set_m1_species_rates(u, ispec, r.kappa_a, r.kappa_s, r.eta_E, r.eta_N, r.kappa_n);
    }
  }

  eos_t eos;
  var_array_t aux;
  double mass_scale;
  bool beta_decay, plasmon_decay, bremsstrahlung, pair_annihilation;
  bool apply_temp_correction;
  bool use_weakhub;
  int betaeq_mode;
  tau_kind_t tau_kind;
  grace::weakhub::device_handle weakhub;
  grace::tau_policy_local_spherical spherical_tau{};
};


} /* namespace grace */

#endif /*GRACE_PHYSICS_EAS_POLICIES_HH*/
