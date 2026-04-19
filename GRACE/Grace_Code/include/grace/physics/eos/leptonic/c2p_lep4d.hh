/**
 * @file c2p_lep4d.hh
 * @brief  Top-level conservative-to-primitive conversion driver for the 4D
 *         leptonic EOS (rho, T, Y_e, Y_mu) in GRACE.
 *
 *         This file mirrors c2p.hh but operates on lep4d_cons_array_t /
 *         lep4d_prims_array_t, and handles the extra YMUSL conservative
 *         variable throughout.
 *
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifndef GRACE_PHYSICS_EOS_C2P_LEP4D_HH
#define GRACE_PHYSICS_EOS_C2P_LEP4D_HH

#include <grace_config.h>
#include <grace/utils/device.h>
#include <grace/utils/metric_utils.hh>
#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/eos/kastaun_c2p_lep4d.hh>
#include <grace/physics/eos/c2p.hh>   // for c2p_err_t, c2p_sig_t, atmo_params_t, etc.
#include <grace/physics/grmhd_helpers.hh>

#include <Kokkos_Core.hpp>

namespace grace {

// ===========================================================================
//  Helper: limit conservatives for 4D leptonic EOS
// ===========================================================================
template< typename eos_t >
KOKKOS_INLINE_FUNCTION
static void limit_conservatives_lep4d(
    lep4d_cons_array_t&    cons,
    metric_array_t const&  metric,
    eos_t const&           eos,
    c2p_err_t&             err)
{
    double rhoL = cons[LEP4D_DENSL] ;
    double yeL  = cons[LEP4D_YESL]  / cons[LEP4D_DENSL] ;
    double ymuL = cons[LEP4D_YMUSL] / cons[LEP4D_DENSL] ;

    // Clamp fractions
    yeL  = Kokkos::fmin(eos.eos_yemax,  Kokkos::fmax(eos.eos_yemin,  yeL))  ;
    ymuL = Kokkos::fmin(eos.eos_ymumax, Kokkos::fmax(eos.eos_ymumin, ymuL)) ;

    double epsmin, epsmax ;
    eos_err_t eoserr ;
    eos.eps_range__rho_ye_ymu(epsmin, epsmax, rhoL, yeL, ymuL, eoserr) ;

    auto B2L = metric.square_vec({cons[LEP4D_BSXL], cons[LEP4D_BSYL], cons[LEP4D_BSZL]}) ;

    if (cons[LEP4D_TAUL]/cons[LEP4D_DENSL] - 0.5*B2L/cons[LEP4D_DENSL]
            < Kokkos::fmin(0., epsmin))
    {
        cons[LEP4D_TAUL] = epsmin*cons[LEP4D_DENSL] + 0.5*B2L ;
        err.set(c2p_err_enum_t::C2P_RESET_TAU) ;
    }

    double st2_max = SQR(cons[LEP4D_TAUL]/cons[LEP4D_DENSL] + 1.) ;
    auto st2L = metric.square_covec({cons[LEP4D_STXL],
                                     cons[LEP4D_STYL],
                                     cons[LEP4D_STZL]}) / SQR(cons[LEP4D_DENSL]) ;
    if (st2L > st2_max) {
        double fac = std::sqrt(st2_max/st2L) ;
        cons[LEP4D_STXL] *= fac ;
        cons[LEP4D_STYL] *= fac ;
        cons[LEP4D_STZL] *= fac ;
        err.set(c2p_err_enum_t::C2P_RESET_STILDE) ;
    }
}

// ===========================================================================
//  Helper: reset to atmosphere for 4D leptonic EOS
// ===========================================================================
template< typename eos_t >
GRACE_HOST_DEVICE
static void reset_to_atmosphere_lep4d(
    lep4d_cons_array_t&    cons,
    lep4d_prims_array_t&   prims,
    metric_array_t const&  metric,
    atmo_params_t const&   atmo,
    eos_t const&           eos,
    c2p_err_t&             err)
{
    prims[LEP4D_RHOL]   = atmo.rho_fl  ;
    prims[LEP4D_YEL]    = atmo.ye_fl   ;
    prims[LEP4D_YMUL]   = eos.get_c2p_ymu_atm() ;
    prims[LEP4D_TEMPL]  = atmo.temp_fl ;
    prims[LEP4D_ZXL]    = 0. ;
    prims[LEP4D_ZYL]    = 0. ;
    prims[LEP4D_ZZL]    = 0. ;

    eos_err_t eoserr ;
    double yel{prims[LEP4D_YEL]}, ymul{prims[LEP4D_YMUL]} ;
    double rhol{prims[LEP4D_RHOL]}, templ{prims[LEP4D_TEMPL]} ;
    double csnd2 ;
    prims[LEP4D_PRESSL] = eos.press_eps_csnd2__temp_rho_ye_ymu(
        prims[LEP4D_EPSL], csnd2, templ, rhol, yel, ymul, eoserr) ;
    prims[LEP4D_ENTL] = eos.entropy_cold__rho_impl(rhol, eoserr) ;

    // recompute conservs consistently
    prims_to_conservs_lep4d(prims, cons, metric) ;
    err.set_all() ;
}

// ===========================================================================
//  primitives → conservatives  (4D leptonic)
// ===========================================================================
GRACE_HOST_DEVICE
inline void prims_to_conservs_lep4d(
    lep4d_prims_array_t&   prims,
    lep4d_cons_array_t&    cons,
    metric_array_t const&  metric)
{
    double const * const betau = metric._beta.data() ;
    double const * const gdd   = metric._g.data()    ;
    double const * const z     = &prims[LEP4D_ZXL]   ;
    double const * const B     = &prims[LEP4D_BXL]   ;
    double const alp           = metric.alp()         ;

    double W ;
    grmhd_get_W(gdd, z, &W) ;

    double smallbu[4], smallb2 ;
    grmhd_get_smallbu_smallb2(betau,gdd,B,z,W,alp,&smallbu,&smallb2) ;

    double D, tau, sstar ;
    double Stilde[3] ;
    grmhd_get_conserved(
        W, prims[LEP4D_RHOL], smallbu, smallb2,
        alp, prims[LEP4D_EPSL], prims[LEP4D_PRESSL],
        betau, z, gdd, prims[LEP4D_ENTL],
        &D, &tau, &Stilde, &sstar) ;

    double const sqrtg = metric.sqrtg() ;
    cons[LEP4D_DENSL]  = sqrtg * D          ;
    cons[LEP4D_STXL]   = sqrtg * Stilde[0]  ;
    cons[LEP4D_STYL]   = sqrtg * Stilde[1]  ;
    cons[LEP4D_STZL]   = sqrtg * Stilde[2]  ;
    cons[LEP4D_TAUL]   = sqrtg * tau         ;
    cons[LEP4D_ENTSL]  = sqrtg * sstar       ;
    cons[LEP4D_YESL]   = cons[LEP4D_DENSL] * prims[LEP4D_YEL]  ;
    cons[LEP4D_YMUSL]  = cons[LEP4D_DENSL] * prims[LEP4D_YMUL] ;
    // magnetic (face-staggered B already handled upstream)
}

// ===========================================================================
//  Main c2p dispatcher  (mirrors conservs_to_prims in c2p.cpp)
// ===========================================================================
template< typename eos_t >
GRACE_HOST_DEVICE
void conservs_to_prims_lep4d(
    lep4d_cons_array_t&    cons,
    lep4d_prims_array_t&   prims,
    metric_array_t const&  metric,
    eos_t const&           eos,
    atmo_params_t const&   atmo,
    excision_params_t const& excision,
    c2p_params_t const&    c2p_pars,
    double*                rtp,
    c2p_err_t&             c2p_err)
{
    using c2p_main_t   = kastaun_c2p_lep4d_t<eos_t>       ;
    using c2p_backup_t = entropy_fix_c2p_lep4d_t<eos_t>   ;

    c2p_sig_t c2p_ret ;
    c2p_err.reset() ;
    c2p_err.set(c2p_err_enum_t::C2P_RESET_ENTROPY) ;

    // Undensitize
    for (auto& c : cons) c /= metric.sqrtg() ;

    // Set B from conservs (B is undensitized cell-centred)
    prims[LEP4D_BXL] = cons[LEP4D_BSXL] ;
    prims[LEP4D_BYL] = cons[LEP4D_BSYL] ;
    prims[LEP4D_BZL] = cons[LEP4D_BSZL] ;

    // Guard negative density
    if (cons[LEP4D_DENSL] < 0.) {
        reset_to_atmosphere_lep4d(cons, prims, metric, atmo, eos, c2p_err) ;
        return ;
    }

    // Set Y_e, Y_mu and clamp
    prims[LEP4D_YEL]  = cons[LEP4D_YESL]  / cons[LEP4D_DENSL] ;
    prims[LEP4D_YMUL] = cons[LEP4D_YMUSL] / cons[LEP4D_DENSL] ;
    {
        bool adj = false ;
        if (prims[LEP4D_YEL]  > eos.eos_yemax)  { prims[LEP4D_YEL]  = eos.eos_yemax  ; adj = true ; }
        if (prims[LEP4D_YEL]  < eos.eos_yemin)  { prims[LEP4D_YEL]  = eos.eos_yemin  ; adj = true ; }
        if (prims[LEP4D_YMUL] > eos.eos_ymumax) { prims[LEP4D_YMUL] = eos.eos_ymumax ; adj = true ; }
        if (prims[LEP4D_YMUL] < eos.eos_ymumin) { prims[LEP4D_YMUL] = eos.eos_ymumin ; adj = true ; }
        if (adj) c2p_err.set(c2p_err_enum_t::C2P_RESET_YE) ;
    }

    bool c2p_is_lenient = (metric.alp() < c2p_pars.alp_bh_thresh) ;

    limit_conservatives_lep4d(cons, metric, eos, c2p_err) ;

    if (cons[LEP4D_DENSL] > atmo.rho_fl * (1. + atmo.atmo_tol)) {

        // --- primary c2p ---
        c2p_main_t c2p(eos, metric, cons) ;
        double residual = c2p.invert(prims, c2p_ret) ;

        c2p_handle_signals(c2p_ret, c2p_is_lenient, c2p_err) ;
        if (rtp) *rtp = residual ;

        bool c2p_failed = c2p_err.any() ;

        // --- entropy backup ---
        if (c2p_failed && c2p_pars.use_entropy_backup) {
            c2p_sig_t backup_ret ;
            c2p_backup_t c2p_bk(eos, metric, cons) ;
            double res_bk = c2p_bk.invert(prims, backup_ret) ;
            c2p_handle_signals(backup_ret, c2p_is_lenient, c2p_err) ;
            if (rtp) *rtp = res_bk ;
        }

    } else {
        reset_to_atmosphere_lep4d(cons, prims, metric, atmo, eos, c2p_err) ;
        return ;
    }

    // Post-c2p: recompute consistent conservs only for flagged quantities
    prims_to_conservs_lep4d(prims, cons, metric) ;
}

// ===========================================================================
//  Explicit template instantiation declarations  (definition in .cpp)
// ===========================================================================
#define INSTANTIATE_LEP4D_C2P(EOS) \
extern template \
void GRACE_HOST_DEVICE \
conservs_to_prims_lep4d<EOS>( lep4d_cons_array_t& \
                             , lep4d_prims_array_t& \
                             , metric_array_t const& \
                             , EOS const& \
                             , atmo_params_t const& \
                             , excision_params_t const& \
                             , c2p_params_t const& \
                             , double* \
                             , c2p_err_t& )

INSTANTIATE_LEP4D_C2P(grace::leptonic_eos_4d_t) ;
#undef INSTANTIATE_LEP4D_C2P

} /* namespace grace */
#endif /* GRACE_PHYSICS_EOS_C2P_LEP4D_HH */
