/**
 * @file c2p_lep4d.hh
 * @brief  Top-level conservative-to-primitive conversion driver for the 4D
 *         leptonic EOS (rho, T, Y_e, Y_mu) in GRACE.
 * @date   2025
 * @copyright This file is part of GRACE. GPL-3 or later.
 */

#ifndef GRACE_PHYSICS_EOS_C2P_LEP4D_HH
#define GRACE_PHYSICS_EOS_C2P_LEP4D_HH

#include <grace_config.h>
#include <grace/physics/eos/c2p.hh>
#include <grace/physics/eos/kastaun_c2p_lep4d.hh>
#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/grmhd_helpers.hh>
#include <grace/utils/device.h>
#include <grace/utils/metric_utils.hh>

#include <Kokkos_Core.hpp>

namespace grace {

template<typename eos_t>
KOKKOS_INLINE_FUNCTION static double
compute_beta_lep4d(lep4d_prims_array_t const& prims,
                   metric_array_t const& metric,
                   eos_t const&)
{
    double const* const betau = metric._beta.data();
    double const* const gdd = metric._g.data();
    double const* const z = &(prims[LEP4D_ZXL]);
    double const* const B = &(prims[LEP4D_BXL]);
    double const alp = metric.alp();

    double W{};
    grmhd_get_W(gdd, z, &W);

    double smallbu[4];
    double smallb2{};
    grmhd_get_smallbu_smallb2(betau, gdd, B, z, W, alp, &smallbu, &smallb2);
    return 2.0 * prims[LEP4D_PRESSL] / Kokkos::fmax(smallb2, 1e-50);
}

template<typename eos_t>
KOKKOS_INLINE_FUNCTION static void
limit_primitives_lep4d(
    lep4d_prims_array_t& prims,
    metric_array_t const& metric,
    eos_t const& eos,
    double max_w,
    double max_sigma,
    c2p_sig_t& sig)
{
    double const z2max = SQR(max_w) - 1.0;
    double const z2 = metric.square_vec(
        {prims[LEP4D_ZXL], prims[LEP4D_ZYL], prims[LEP4D_ZZL]});
    if (z2 >= z2max) {
        double const znorm = 0.99 * Kokkos::sqrt(z2max / z2);
        prims[LEP4D_ZXL] *= znorm;
        prims[LEP4D_ZYL] *= znorm;
        prims[LEP4D_ZZL] *= znorm;
        sig.set(c2p_sig_enum_t::C2P_VEL_TOO_HIGH);
    }

    double const* const betau = metric._beta.data();
    double const* const gdd = metric._g.data();
    double const* const z = &(prims[LEP4D_ZXL]);
    double const* const B = &(prims[LEP4D_BXL]);
    double const alp = metric.alp();

    double W{};
    grmhd_get_W(gdd, z, &W);

    double smallbu[4];
    double smallb2{};
    grmhd_get_smallbu_smallb2(betau, gdd, B, z, W, alp, &smallbu, &smallb2);
    double const sigma = smallb2 / prims[LEP4D_RHOL];
    if (sigma >= max_sigma) {
        double const rhofact = 1.001 * sigma / max_sigma;
        prims[LEP4D_RHOL] *= rhofact;
        eos_err_t eos_err;
        double csnd2{};
        prims[LEP4D_PRESSL] = eos.press_eps_csnd2_entropy__temp_rho_ye_ymu(
            prims[LEP4D_EPSL], csnd2, prims[LEP4D_ENTL], prims[LEP4D_TEMPL],
            prims[LEP4D_RHOL], prims[LEP4D_YEL], prims[LEP4D_YMUL], eos_err);
        sig.set(c2p_sig_enum_t::C2P_SIGMA_TOO_HIGH);
    }
}

template<typename eos_t>
KOKKOS_INLINE_FUNCTION static void
limit_conservatives_lep4d(
    lep4d_cons_array_t& cons,
    metric_array_t const& metric,
    eos_t const& eos,
    c2p_err_t& err)
{
    double rhoL = cons[LEP4D_DENSL];
    double yeL = cons[LEP4D_YESL] / cons[LEP4D_DENSL];
    double ymuL = cons[LEP4D_YMUSL] / cons[LEP4D_DENSL];

    yeL = Kokkos::fmin(eos.get_c2p_ye_max(), Kokkos::fmax(eos.get_c2p_ye_min(), yeL));
    ymuL = Kokkos::fmin(eos.get_c2p_ymu_max(), Kokkos::fmax(eos.get_c2p_ymu_min(), ymuL));

    double epsmin{}, epsmax{};
    eos_err_t eoserr;
    eos.eps_range__rho_ye_ymu(epsmin, epsmax, rhoL, yeL, ymuL, eoserr);

    auto B2L = metric.square_vec(
        {cons[LEP4D_BSXL], cons[LEP4D_BSYL], cons[LEP4D_BSZL]});

    if (cons[LEP4D_TAUL] / cons[LEP4D_DENSL] - 0.5 * B2L / cons[LEP4D_DENSL]
        < Kokkos::fmin(0.0, epsmin))
    {
        cons[LEP4D_TAUL] = epsmin * cons[LEP4D_DENSL] + 0.5 * B2L;
        err.set(c2p_err_enum_t::C2P_RESET_TAU);
    }

    double const st2_max = SQR(cons[LEP4D_TAUL] / cons[LEP4D_DENSL] + 1.0);
    auto const st2L = metric.square_covec(
        {cons[LEP4D_STXL], cons[LEP4D_STYL], cons[LEP4D_STZL]}) / SQR(cons[LEP4D_DENSL]);
    if (st2L > st2_max) {
        double const fac = Kokkos::sqrt(st2_max / st2L);
        cons[LEP4D_STXL] *= fac;
        cons[LEP4D_STYL] *= fac;
        cons[LEP4D_STZL] *= fac;
        err.set(c2p_err_enum_t::C2P_RESET_STILDE);
    }
}

template<typename eos_t>
GRACE_HOST_DEVICE static void
reset_to_atmosphere_lep4d(
    lep4d_cons_array_t& cons,
    lep4d_prims_array_t& prims,
    metric_array_t const& metric,
    atmo_params_t const& atmo,
    eos_t const& eos,
    c2p_err_t& err)
{
    prims[LEP4D_RHOL] = atmo.rho_fl;
    prims[LEP4D_YEL] = atmo.ye_fl;
    prims[LEP4D_YMUL] = eos.get_c2p_ymu_atm();
    prims[LEP4D_TEMPL] = atmo.temp_fl;
    prims[LEP4D_ZXL] = 0.0;
    prims[LEP4D_ZYL] = 0.0;
    prims[LEP4D_ZZL] = 0.0;

    eos_err_t eoserr;
    double csnd2{};
    prims[LEP4D_PRESSL] = eos.press_eps_csnd2_entropy__temp_rho_ye_ymu(
        prims[LEP4D_EPSL], csnd2, prims[LEP4D_ENTL], prims[LEP4D_TEMPL],
        prims[LEP4D_RHOL], prims[LEP4D_YEL], prims[LEP4D_YMUL], eoserr);

    prims_to_conservs_lep4d(prims, cons, metric);
    err.set_all();
}

GRACE_HOST_DEVICE inline void
prims_to_conservs_lep4d(
    lep4d_prims_array_t& prims,
    lep4d_cons_array_t& cons,
    metric_array_t const& metric)
{
    double const* const betau = metric._beta.data();
    double const* const gdd = metric._g.data();
    double const* const z = &prims[LEP4D_ZXL];
    double const* const B = &prims[LEP4D_BXL];
    double const alp = metric.alp();

    double W{};
    grmhd_get_W(gdd, z, &W);

    double smallbu[4], smallb2{};
    grmhd_get_smallbu_smallb2(betau, gdd, B, z, W, alp, &smallbu, &smallb2);

    double D{}, tau{}, sstar{};
    double Stilde[3];
    grmhd_get_conserved(
        W, prims[LEP4D_RHOL], smallbu, smallb2,
        alp, prims[LEP4D_EPSL], prims[LEP4D_PRESSL],
        betau, z, gdd, prims[LEP4D_ENTL],
        &D, &tau, &Stilde, &sstar);

    double const sqrtg = metric.sqrtg();
    cons[LEP4D_DENSL] = sqrtg * D;
    cons[LEP4D_STXL] = sqrtg * Stilde[0];
    cons[LEP4D_STYL] = sqrtg * Stilde[1];
    cons[LEP4D_STZL] = sqrtg * Stilde[2];
    cons[LEP4D_TAUL] = sqrtg * tau;
    cons[LEP4D_ENTSL] = sqrtg * sstar;
    cons[LEP4D_YESL] = cons[LEP4D_DENSL] * prims[LEP4D_YEL];
    cons[LEP4D_YMUSL] = cons[LEP4D_DENSL] * prims[LEP4D_YMUL];
}

template<typename eos_t>
GRACE_HOST_DEVICE void
conservs_to_prims_lep4d(
    lep4d_cons_array_t& cons,
    lep4d_prims_array_t& prims,
    metric_array_t const& metric,
    eos_t const& eos,
    atmo_params_t const& atmo,
    excision_params_t const& excision,
    c2p_params_t const& c2p_pars,
    double* rtp,
    c2p_err_t& c2p_err)
{
    using c2p_main_t = kastaun_c2p_lep4d_t<eos_t>;
    using c2p_backup_t = entropy_fix_c2p_lep4d_t<eos_t>;

    bool c2p_failed{false};
    c2p_sig_t c2p_ret;
    c2p_err.reset();
    c2p_err.set(c2p_err_enum_t::C2P_RESET_ENTROPY);

    for (auto& c : cons) {
        c /= metric.sqrtg();
    }

    prims[LEP4D_BXL] = cons[LEP4D_BSXL];
    prims[LEP4D_BYL] = cons[LEP4D_BSYL];
    prims[LEP4D_BZL] = cons[LEP4D_BSZL];

    if (cons[LEP4D_DENSL] < 0.0) {
        reset_to_atmosphere_lep4d(cons, prims, metric, atmo, eos, c2p_err);
        return;
    }

    prims[LEP4D_YEL] = cons[LEP4D_YESL] / cons[LEP4D_DENSL];
    prims[LEP4D_YMUL] = cons[LEP4D_YMUSL] / cons[LEP4D_DENSL];
    if (prims[LEP4D_YEL] > eos.get_c2p_ye_max()) {
        prims[LEP4D_YEL] = eos.get_c2p_ye_max();
        c2p_err.set(c2p_err_enum_t::C2P_RESET_YE);
    }
    if (prims[LEP4D_YEL] < eos.get_c2p_ye_min()) {
        prims[LEP4D_YEL] = eos.get_c2p_ye_min();
        c2p_err.set(c2p_err_enum_t::C2P_RESET_YE);
    }
    if (prims[LEP4D_YMUL] > eos.get_c2p_ymu_max()) {
        prims[LEP4D_YMUL] = eos.get_c2p_ymu_max();
        c2p_err.set(c2p_err_enum_t::C2P_RESET_YE);
    }
    if (prims[LEP4D_YMUL] < eos.get_c2p_ymu_min()) {
        prims[LEP4D_YMUL] = eos.get_c2p_ymu_min();
        c2p_err.set(c2p_err_enum_t::C2P_RESET_YE);
    }

    bool const c2p_is_lenient = metric.alp() < c2p_pars.alp_bh_thresh;
    limit_conservatives_lep4d(cons, metric, eos, c2p_err);

    if (cons[LEP4D_DENSL] > atmo.rho_fl * (1.0 + atmo.atmo_tol)) {
        c2p_main_t c2p(eos, metric, cons);
        double const residual = c2p.invert(prims, c2p_ret);

        c2p_failed = (Kokkos::fabs(residual) > c2p_pars.tol)
                  || c2p_ret.test(c2p_sig_enum_t::C2P_RHO_TOO_HIGH)
                  || c2p_ret.test(c2p_sig_enum_t::C2P_EPS_TOO_HIGH)
                  || (!Kokkos::isfinite(residual));

        double const beta = compute_beta_lep4d(prims, metric, eos);
        if ((c2p_failed || beta <= c2p_pars.beta_fallback) && c2p_pars.use_ent_backup) {
            c2p_ret.reset();
            c2p_backup_t backup(eos, metric, cons);
            double const backup_residual = backup.invert(prims, c2p_ret);
            c2p_failed = (Kokkos::fabs(backup_residual) > c2p_pars.tol)
                      || (!Kokkos::isfinite(backup_residual))
                      || (prims[LEP4D_EPSL] >= eos.get_c2p_eps_max());
        }

        c2p_handle_signals(c2p_ret, c2p_is_lenient, c2p_err);
    } else {
        c2p_failed = true;
    }

    bool const excise = excision.excise_by_radius
        ? rtp[0] <= excision.r_ex
        : metric.alp() <= excision.alp_ex;

    if ((prims[LEP4D_RHOL] < (1.0 + atmo.atmo_tol) * atmo.rho_fl) || c2p_failed || excise) {
        reset_to_atmosphere_lep4d(cons, prims, metric, atmo, eos, c2p_err);
    } else if (prims[LEP4D_TEMPL] < atmo.temp_fl) {
        prims[LEP4D_TEMPL] = atmo.temp_fl;
        eos_err_t eos_err;
        double csnd2{};
        prims[LEP4D_PRESSL] = eos.press_eps_csnd2_entropy__temp_rho_ye_ymu(
            prims[LEP4D_EPSL], csnd2, prims[LEP4D_ENTL], prims[LEP4D_TEMPL],
            prims[LEP4D_RHOL], prims[LEP4D_YEL], prims[LEP4D_YMUL], eos_err);
        c2p_handle_eos_signals(eos_err, c2p_is_lenient, c2p_err);
        c2p_err.set(c2p_err_enum_t::C2P_RESET_TAU);
        c2p_err.set(c2p_err_enum_t::C2P_RESET_STILDE);
        c2p_err.set(c2p_err_enum_t::C2P_RESET_ENTROPY);
    } else {
        c2p_sig_t c2p_sig;
        limit_primitives_lep4d(
            prims, metric, eos, c2p_pars.max_w, c2p_pars.max_sigma, c2p_sig);
        c2p_handle_signals(c2p_sig, c2p_is_lenient, c2p_err);
    }

    prims_to_conservs_lep4d(prims, cons, metric);
}

#define INSTANTIATE_LEP4D_C2P(EOS) \
extern template \
void GRACE_HOST_DEVICE \
conservs_to_prims_lep4d<EOS>(lep4d_cons_array_t&, \
                             lep4d_prims_array_t&, \
                             metric_array_t const&, \
                             EOS const&, \
                             atmo_params_t const&, \
                             excision_params_t const&, \
                             c2p_params_t const&, \
                             double*, \
                             c2p_err_t&)

INSTANTIATE_LEP4D_C2P(grace::leptonic_eos_4d_t);
#undef INSTANTIATE_LEP4D_C2P

} /* namespace grace */

#endif /* GRACE_PHYSICS_EOS_C2P_LEP4D_HH */
