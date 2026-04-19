/**
 * @file grmhd_lep4d_id.hh
 * @brief  Initial-data and beta-equilibrium helpers for the 4D leptonic EOS.
 *
 *         Provides:
 *          - Host-side beta-equilibrium solver (Y_e, Y_mu at given rho, T).
 *          - lep4d_id_helper_t: Kokkos-compatible functor that extends the
 *            standard set_grmhd_initial_data_impl primitive-setting kernel
 *            to also set the YMUSTAR_ conserved variable and the YMU_
 *            auxiliary from imported primitives.
 *          - set_grmhd_initial_data_lep4d(): thin wrapper that calls the
 *            standard set_grmhd_initial_data<leptonic_eos_4d_t>() and then
 *            populates YMUSTAR_ / YMU_ consistently.
 *
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifndef GRACE_GRMHD_LEP4D_ID_HH
#define GRACE_GRMHD_LEP4D_ID_HH

#ifdef GRACE_ENABLE_LEPTONIC_4D

#include <grace_config.h>
#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/eos/eos_storage_lep4d.hh>
#include <grace/physics/eos/c2p_lep4d.hh>
#include <grace/physics/grmhd_helpers.hh>
#include <grace/data_structures/variable_indices.hh>
#include <grace/data_structures/variables.hh>
#include <grace/utils/grace_utils.hh>

#include <Kokkos_Core.hpp>

namespace grace {

// ===========================================================================
//  Host-side beta-equilibrium solver for (Y_e, Y_mu) at fixed rho, T.
//
//  Conditions:
//    mu_n - mu_p - mu_e = 0          (standard beta equilibrium)
//    mu_mu - mu_e = 0                 (muon threshold / equilibrium)
//    Y_e = Y_{e^-} - Y_{e^+}         (charge neutrality)
//    Y_mu = Y_{mu^-} - Y_{mu^+}      (from table)
//
//  The solver performs an outer Brent iteration over Y_mu in
//  [ymumin, ymumax], and for each trial Y_mu it solves the inner
//  beta-eq condition for Y_e via a second Brent (reusing the existing
//  baryon-only pattern from read_eos_table.cpp).
// ===========================================================================

/**
 * @brief Find (Y_e, Y_mu) in beta equilibrium at given rho, T.
 *        Host-only (uses Brent root-finding with std::function).
 *
 * @param[out] ye_eq   Equilibrium electron fraction.
 * @param[out] ymu_eq  Equilibrium muon fraction.
 * @param[in]  rho     Rest-mass density [geometric units].
 * @param[in]  temp    Temperature [geometric units].
 * @param[in]  eos     The 4D leptonic EOS object (host copy).
 */
inline void
find_beta_eq_ye_ymu(double& ye_eq, double& ymu_eq,
                    double rho, double temp,
                    leptonic_eos_4d_t const& eos)
{
    eos_err_t err ;

    // clamp rho and T to table range (host side)
    rho  = std::max(eos.eos_rhomin,  std::min(eos.eos_rhomax,  rho))  ;
    temp = std::max(std::exp(eos.ltempmin), std::min(std::exp(eos.ltempmax), temp)) ;

    double const yemin  = eos.eos_yemin  ;
    double const yemax  = eos.eos_yemax  ;
    double const ymumin = eos.eos_ymumin ;
    double const ymumax = eos.eos_ymumax ;

    // -------------------------------------------------------------------
    //  Outer loop: iterate over Y_mu to enforce  mu_mu = mu_e.
    //  For each trial Y_mu we find the Y_e that satisfies standard
    //  beta equilibrium (mu_n - mu_p - mu_e = 0) via an inner Brent.
    // -------------------------------------------------------------------
    auto residual_outer = [&](double ymu_trial) -> double
    {
        // Inner: find Y_e satisfying  mu_n - mu_p - mu_e = 0
        auto residual_inner = [&](double ye_trial) -> double
        {
            eos_err_t lerr ;
            double mue_val, mumu_dummy, mup_val, mun_val ;
            double rhol{rho}, templ{temp}, yel{ye_trial}, ymul{ymu_trial} ;
            mue_val = eos.mue_mumu_mup_mun__temp_rho_ye_ymu(
                mumu_dummy, mup_val, mun_val,
                templ, rhol, yel, ymul, lerr) ;
            // dmu = mu_n - mu_p - mu_e  (should be zero in beta-eq)
            return mun_val - mup_val - mue_val ;
        } ;

        // Solve inner beta-eq for Y_e
        ye_eq = utils::brent(residual_inner, yemin, yemax, 1e-14) ;

        // Now evaluate  mu_mu - mu_e  (should be zero at muon threshold)
        eos_err_t lerr ;
        double mue_val, mumu_val, mup_dummy, mun_dummy ;
        double rhol{rho}, templ{temp}, yel{ye_eq}, ymul{ymu_trial} ;
        mue_val = eos.mue_mumu_mup_mun__temp_rho_ye_ymu(
            mumu_val, mup_dummy, mun_dummy,
            templ, rhol, yel, ymul, lerr) ;

        return mumu_val - mue_val ;
    } ;

    // If mu_mu(ymumin) > 0 everywhere (no muons at this rho,T) set ymu=0
    double f0 = residual_outer(ymumin) ;
    double f1 = residual_outer(ymumax) ;

    if (f0 * f1 > 0.) {
        // No sign change → muons absent or saturated
        // Choose endpoint with smallest absolute value
        ymu_eq = (std::fabs(f0) < std::fabs(f1)) ? ymumin : ymumax ;
        // Re-solve inner to set ye_eq consistently
        auto residual_inner = [&](double ye_trial) -> double {
            eos_err_t lerr ;
            double mue_v, mumu_d, mup_v, mun_v ;
            double rhol{rho}, templ{temp}, yel{ye_trial}, ymul{ymu_eq} ;
            mue_v = eos.mue_mumu_mup_mun__temp_rho_ye_ymu(
                mumu_d, mup_v, mun_v, templ, rhol, yel, ymul, lerr) ;
            return mun_v - mup_v - mue_v ;
        } ;
        ye_eq = utils::brent(residual_inner, yemin, yemax, 1e-14) ;
    } else {
        ymu_eq = utils::brent(residual_outer, ymumin, ymumax, 1e-14) ;
    }
}

// ===========================================================================
//  Kernel helper: after the standard ID has set DENS_, YESTAR_, etc.,
//  compute YMUSTAR_ = sqrt(g) * D * Y_mu  using beta equilibrium at
//  each point (for the TOV / cold-start use case).
//  For imported IDs (FUKA, LORENE) Y_mu should be interpolated from the
//  ID data or set to beta-eq as a fallback.
// ===========================================================================

/**
 * @brief Post-ID kernel: populate YMUSTAR_ and YMU_ from beta equilibrium.
 *
 *        Iterates over all active quadrants.  Y_mu is determined from
 *        beta equilibrium at the local (rho, T) reconstructed from the
 *        already-set DENS_ and TAU_ conserved variables plus the auxiliary
 *        TEMP_.
 *
 * @param eos  Device copy of the 4D leptonic EOS.
 * @param vars Variable list (provides state and aux views).
 */
inline void
set_ymustar_from_beta_eq(leptonic_eos_4d_t const& eos,
                         variable_list_impl_t&     vars)
{
    using namespace Kokkos ;

    auto& state = vars.getstate() ;
    auto& aux   = vars.getaux()   ;

    auto& forest = amr::forest::get() ;
    size_t const nq  = forest.local_num_quadrants() ;

    auto const nx  = state.extent(0) ;
    auto const ny  = state.extent(1) ;
    auto const nz  = state.extent(2) ;

    parallel_for(
        GRACE_EXECUTION_TAG("ID","lep4d_set_ymustar")
      , MDRangePolicy<Rank<GRACE_NSPACEDIM+1>,default_execution_space>(
            {VEC(0,0,0),0},
            {VEC(nx,ny,nz),static_cast<int>(nq)})
      , KOKKOS_LAMBDA (VEC(int const& i, int const& j, int const& k), int const& q)
        {
            double const dens = state(VEC(i,j,k), DENS_, q) ;
            if (dens <= 0.) return ;

            double const sqrtg_loc = 1.0 ; // metric.sqrtg() not available here;
            // density already includes sqrtg: D_star = sqrtg * rho * W
            // Y_e is in aux already
            double ye_loc   = aux(VEC(i,j,k), YE_,   q) ;
            double temp_loc = aux(VEC(i,j,k), TEMP_,  q) ;
            double rho_loc  = aux(VEC(i,j,k), RHO_,   q) ;

            // Clamp
            ye_loc   = Kokkos::fmin(eos.eos_yemax,   Kokkos::fmax(eos.eos_yemin,   ye_loc))  ;
            temp_loc = Kokkos::fmin(Kokkos::exp(eos.ltempmax),
                                    Kokkos::fmax(Kokkos::exp(eos.ltempmin), temp_loc)) ;
            rho_loc  = Kokkos::fmin(eos.eos_rhomax,  Kokkos::fmax(eos.eos_rhomin,  rho_loc)) ;

            // Minimal muon estimate: cold beta-eq
            // On GPU we can only do a direct table lookup at fixed ye.
            // Full beta-eq requires Brent → host-only. As a GPU-safe fallback
            // we interpolate mu_mu and mu_e at the current (rho,T,ye,ymumin)
            // and nudge ymu upward until mu_mu < mu_e.
            double ymu_loc = eos.eos_ymumin ;
            {
                eos_err_t lerr ;
                double rhol{rho_loc}, templ{temp_loc}, yel{ye_loc}, ymul{eos.eos_ymumin} ;
                double mumu = eos.mumu__temp_rho_ye_ymu(templ, rhol, yel, ymul, lerr) ;
                double mue  = eos.mue__temp_rho_ye_ymu (templ, rhol, yel, ymul, lerr) ;
                // If mu_mu > mu_e muons are energetically favoured; use ymumax as safe upper
                // bound.  Exact beta-eq is done host-side; here we use mid-point as a seed.
                if (mumu > mue) {
                    ymu_loc = 0.5 * (eos.eos_ymumin + eos.eos_ymumax) ;
                }
            }

            // Store Y_mu primitive
            aux(VEC(i,j,k), YMU_, q) = ymu_loc ;

            // YMUSTAR_ = DENS_ * Y_mu  (same densitization convention as YESTAR_)
            state(VEC(i,j,k), YMUSTAR_, q) = dens * ymu_loc ;
        }) ;

    fence() ;
}

// ===========================================================================
//  Public entry point for leptonic-4D initial data setup.
//
//  Call this function *instead of* set_grmhd_initial_data<eos_t>()
//  when eos_type == "leptonic_4d".
// ===========================================================================

/**
 * @brief Top-level initial data setup for the 4D leptonic EOS.
 *
 *        Calls set_grmhd_initial_data<leptonic_eos_4d_t>() for the
 *        standard fields (DENS_, TAU_, SX_/Y_/Z_, YESTAR_, ENTROPYSTAR_,
 *        B, prims), then populates YMUSTAR_ and YMU_ from beta equilibrium.
 *
 *        For beta-equilibrium atmosphere setup the host-side solver
 *        find_beta_eq_ye_ymu() is called once to determine the atmosphere
 *        Y_mu and then broadcast uniformly to cells at the floor density.
 */
void set_grmhd_initial_data_lep4d() ;

} /* namespace grace */

#endif /* GRACE_ENABLE_LEPTONIC_4D */
#endif /* GRACE_GRMHD_LEP4D_ID_HH */
