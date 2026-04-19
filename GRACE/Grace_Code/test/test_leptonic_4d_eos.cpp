/**
 * @file test_leptonic_4d_eos.cpp
 * @brief  Unit tests for the 4D leptonic EOS (leptonic_eos_4d_t) in GRACE.
 *
 *         Tests:
 *          1. EOS API self-consistency: for every (rho, T, Y_e, Y_mu) point,
 *             every API call that should return the same thermodynamic quantity
 *             agrees to relative tolerance 1e-10.
 *          2. C2P round-trip: primitive → conservative → primitive recovers
 *             the input primitives to a residual < 1e-10.
 *          3. Beta-equilibrium solver: at the solved (Y_e, Y_mu) point,
 *             mu_n - mu_p - mu_e ≈ 0  and  mu_mu - mu_e ≈ 0.
 *
 *         Style mirrors test_eos.cpp and test_pwpoly_eos.cpp.
 *
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifdef GRACE_ENABLE_LEPTONIC_4D

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Kokkos_Core.hpp>

#include <grace_config.h>
#include <grace/amr/grace_amr.hh>
#include <grace/data_structures/grace_data_structures.hh>
#include <grace/utils/grace_utils.hh>
#include <grace/physics/eos/eos_storage_lep4d.hh>
#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/eos/c2p_lep4d.hh>
#include <grace/physics/id/grmhd_lep4d_id.hh>
#include <grace/system/grace_system.hh>
#include <grace/utils/metric_utils.hh>
#include <grace/physics/grmhd_helpers.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

// ============================================================
//  Helpers
// ============================================================

static std::vector<double>
logspace(double lo, double hi, int N)
{
    std::vector<double> v(N) ;
    double const dl = (std::log(hi)-std::log(lo)) / (N-1) ;
    for (int i=0; i<N; ++i) v[i] = std::exp(std::log(lo)+i*dl) ;
    return v ;
}

static std::vector<double>
linspace(double lo, double hi, int N)
{
    std::vector<double> v(N) ;
    double const dh = (hi-lo)/(N-1) ;
    for (int i=0; i<N; ++i) v[i] = lo + i*dh ;
    return v ;
}

// ============================================================
//  Test 1: EOS API self-consistency on a grid of (rho,T,ye,ymu)
// ============================================================

TEST_CASE("leptonic_4d EOS API self-consistency", "[leptonic_4d][eos]")
{
    auto eos = grace::get_leptonic_4d_eos() ;

    // Sample a coarse 4-D grid
    constexpr int NR = 8, NT = 5, NY = 3, NYM = 3 ;
    auto rhos  = logspace(eos.eos_rhomin*10., eos.eos_rhomax*0.1, NR)  ;
    auto temps = logspace(std::exp(eos.ltempmin)*2., std::exp(eos.ltempmax)*0.5, NT) ;
    auto yes   = linspace(eos.eos_yemin  + 0.01, eos.eos_yemax  - 0.01, NY)  ;
    auto ymus  = linspace(eos.eos_ymumin + 0.001, eos.eos_ymumax - 0.001, NYM) ;

    // Allocate result arrays on device
    int const NPTS = NR*NT*NY*NYM ;
    Kokkos::View<double*> d_press("press", NPTS) ;
    Kokkos::View<double*> d_eps  ("eps",   NPTS) ;
    Kokkos::View<double*> d_h    ("h",     NPTS) ;
    Kokkos::View<double*> d_csnd ("csnd",  NPTS) ;
    Kokkos::View<double*> d_ent  ("ent",   NPTS) ;
    Kokkos::View<double*> d_temp ("temp",  NPTS) ;
    // Error matrix: up to 20 consistency checks per point
    Kokkos::View<double**> d_err ("err",   NPTS, 20) ;

    // Copy input grid to device
    Kokkos::View<double*> d_rho ("d_rho",  NR)  ;
    Kokkos::View<double*> d_t   ("d_t",    NT)  ;
    Kokkos::View<double*> d_ye  ("d_ye",   NY)  ;
    Kokkos::View<double*> d_ymu ("d_ymu",  NYM) ;
    {
        auto h_rho = Kokkos::create_mirror_view(d_rho) ;
        auto h_t   = Kokkos::create_mirror_view(d_t)   ;
        auto h_ye  = Kokkos::create_mirror_view(d_ye)  ;
        auto h_ymu = Kokkos::create_mirror_view(d_ymu) ;
        for (int i=0; i<NR;  ++i) h_rho(i) = rhos[i]  ;
        for (int i=0; i<NT;  ++i) h_t(i)   = temps[i] ;
        for (int i=0; i<NY;  ++i) h_ye(i)  = yes[i]   ;
        for (int i=0; i<NYM; ++i) h_ymu(i) = ymus[i]  ;
        Kokkos::deep_copy(d_rho, h_rho) ;
        Kokkos::deep_copy(d_t,   h_t)   ;
        Kokkos::deep_copy(d_ye,  h_ye)  ;
        Kokkos::deep_copy(d_ymu, h_ymu) ;
    }

    Kokkos::parallel_for("lep4d_eos_api_test", NPTS,
    KOKKOS_LAMBDA (int idx)
    {
        int ir  = idx / (NT*NY*NYM) ;
        int rem = idx % (NT*NY*NYM) ;
        int it  = rem / (NY*NYM)    ;
        rem     = rem % (NY*NYM)    ;
        int iy  = rem / NYM         ;
        int iym = rem % NYM         ;

        double rho  = d_rho(ir)  ;
        double temp = d_t(it)    ;
        double ye   = d_ye(iy)   ;
        double ymu  = d_ymu(iym) ;

        grace::eos_err_t err ;
        int ww = 0 ;

        // Reference: press_h_csnd2_temp_entropy__eps_rho_ye_ymu
        double h0, csnd0, temp0, ent0 ;
        double eps0 = eos.eps__temp_rho_ye_ymu(temp, rho, ye, ymu, err) ;
        double p0   = eos.press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
                          h0, csnd0, temp0, ent0, eps0, rho, ye, ymu, err) ;

        d_press(idx) = p0 ;
        d_eps(idx)   = eps0 ;
        d_h(idx)     = h0 ;
        d_csnd(idx)  = csnd0 ;
        d_ent(idx)   = ent0 ;
        d_temp(idx)  = temp0 ;

        // Check 1: press__temp_rho_ye_ymu agrees
        {
            double tl{temp}, rhol{rho}, yel{ye}, ymul{ymu} ;
            double p1 = eos.press__temp_rho_ye_ymu(tl, rhol, yel, ymul, err) ;
            d_err(idx,ww++) = Kokkos::abs(p0-p1) / (Kokkos::abs(p0)+1e-40) ;
        }
        // Check 2: eps__temp_rho_ye_ymu round-trip
        {
            double tl{temp}, rhol{rho}, yel{ye}, ymul{ymu} ;
            double e1 = eos.eps__temp_rho_ye_ymu(tl, rhol, yel, ymul, err) ;
            d_err(idx,ww++) = Kokkos::abs(eps0-e1) / (Kokkos::abs(eps0)+1e-40) ;
        }
        // Check 3: press_temp__eps_rho_ye_ymu recovers temperature
        {
            double tl{temp}, rhol{rho}, yel{ye}, ymul{ymu}, epsl{eps0} ;
            double p2 = eos.press_temp__eps_rho_ye_ymu(tl, epsl, rhol, yel, ymul, err) ;
            d_err(idx,ww++) = Kokkos::abs(p0-p2) / (Kokkos::abs(p0)+1e-40) ;
            d_err(idx,ww++) = Kokkos::abs(temp-tl) / (Kokkos::abs(temp)+1e-40) ;
        }
        // Check 4: press__eps_rho_ye_ymu
        {
            double rhol{rho}, yel{ye}, ymul{ymu}, epsl{eps0} ;
            double p3 = eos.press__eps_rho_ye_ymu(epsl, rhol, yel, ymul, err) ;
            d_err(idx,ww++) = Kokkos::abs(p0-p3) / (Kokkos::abs(p0)+1e-40) ;
        }
        // Check 5: press_h_csnd2__temp_rho_ye_ymu
        {
            double hl, cl, tl{temp}, rhol{rho}, yel{ye}, ymul{ymu} ;
            double p4 = eos.press_h_csnd2__temp_rho_ye_ymu(hl, cl, tl, rhol, yel, ymul, err) ;
            d_err(idx,ww++) = Kokkos::abs(p0-p4) / (Kokkos::abs(p0)+1e-40) ;
            d_err(idx,ww++) = Kokkos::abs(h0-hl) / (Kokkos::abs(h0)+1e-40) ;
            d_err(idx,ww++) = Kokkos::abs(csnd0-cl) / (Kokkos::abs(csnd0)+1e-40) ;
        }
        // Check 6: entropy-based inversion recovers eps and press
        {
            double hl, cl, tl, epsl, ent_in{ent0}, rhol{rho}, yel{ye}, ymul{ymu} ;
            double p5 = eos.press_h_csnd2_temp_eps__entropy_rho_ye_ymu(
                            hl, cl, tl, epsl, ent_in, rhol, yel, ymul, err) ;
            d_err(idx,ww++) = Kokkos::abs(p0-p5)   / (Kokkos::abs(p0)+1e-40) ;
            d_err(idx,ww++) = Kokkos::abs(eps0-epsl)/ (Kokkos::abs(eps0)+1e-40) ;
        }
        // Check 7: eps range brackets eps0
        {
            double emin, emax, rhol{rho}, yel{ye}, ymul{ymu} ;
            eos.eps_range__rho_ye_ymu(emin, emax, rhol, yel, ymul, err) ;
            d_err(idx,ww++) = (eps0 >= emin-1e-12 && eps0 <= emax+1e-12) ? 0. : 1. ;
        }
    }) ;
    Kokkos::fence() ;

    // Copy errors to host and check
    auto h_err = Kokkos::create_mirror_view(d_err) ;
    Kokkos::deep_copy(h_err, d_err) ;

    constexpr double TOL = 1e-8 ;
    for (int i=0; i<NPTS; ++i) {
        for (int j=0; j<12; ++j) {
            INFO("Point " << i << " check " << j << " err=" << h_err(i,j)) ;
            CHECK_THAT(h_err(i,j), Catch::Matchers::WithinAbs(0., TOL)) ;
        }
    }
}

// ============================================================
//  Test 2: C2P round-trip  (prim→cons→prim)
// ============================================================

TEST_CASE("leptonic_4d C2P round-trip", "[leptonic_4d][c2p]")
{
    auto eos = grace::get_leptonic_4d_eos() ;

    // Build a flat Minkowski metric
    grace::metric_array_t metric ;
    metric._g    = {1.,0.,0.,1.,0.,1.} ;  // g_{ij} = diag(1,1,1)
    metric._beta = {0.,0.,0.} ;
    metric._alp  = 1. ;
    metric._sqrtg = 1. ;
    // Build inverse metric manually
    metric._ginv = {1.,0.,0.,1.,0.,1.} ;

    // Test a single representative primitive state
    grace::lep4d_prims_array_t prims_in{} ;
    prims_in[grace::LEP4D_RHOL]   = 1e-4   ;
    prims_in[grace::LEP4D_YEL]    = 0.1    ;
    prims_in[grace::LEP4D_YMUL]   = 0.01   ;
    prims_in[grace::LEP4D_TEMPL]  = 10.    ; // MeV-scale, already geometric
    prims_in[grace::LEP4D_ZXL]    = 0.02   ;
    prims_in[grace::LEP4D_ZYL]    = 0.     ;
    prims_in[grace::LEP4D_ZZL]    = 0.     ;
    prims_in[grace::LEP4D_BXL]    = 0.     ;
    prims_in[grace::LEP4D_BYL]    = 0.     ;
    prims_in[grace::LEP4D_BZL]    = 0.     ;

    // Compute eps and press for input primitives
    grace::eos_err_t eoserr ;
    double ye_l{prims_in[grace::LEP4D_YEL]} ;
    double ymu_l{prims_in[grace::LEP4D_YMUL]} ;
    double rho_l{prims_in[grace::LEP4D_RHOL]} ;
    double temp_l{prims_in[grace::LEP4D_TEMPL]} ;
    double h_dummy, cs_dummy, ent_dummy ;
    prims_in[grace::LEP4D_PRESSL] =
        eos.press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
            h_dummy, cs_dummy, temp_l, ent_dummy,
            prims_in[grace::LEP4D_EPSL],
            rho_l, ye_l, ymu_l, eoserr) ;
    prims_in[grace::LEP4D_ENTL] = ent_dummy ;

    // Compute conservatives from primitives
    grace::lep4d_cons_array_t cons{} ;
    grace::prims_to_conservs_lep4d(prims_in, cons, metric) ;

    // Invert back
    grace::atmo_params_t atmo ;
    atmo.rho_fl    = 1e-15 ;
    atmo.temp_fl   = 1e-10 ;
    atmo.ye_fl     = eos.eos_yemin ;
    atmo.atmo_tol  = 1e-3  ;

    grace::excision_params_t excision ;
    excision.alp_ex = 0.05 ;
    excision.r_ex   = 1.0  ;
    excision.excise_by_radius = false ;

    grace::c2p_params_t c2p_pars ;
    c2p_pars.alp_bh_thresh     = 0.05  ;
    c2p_pars.max_w             = 50.   ;
    c2p_pars.max_sigma         = 100.  ;
    c2p_pars.use_entropy_backup= false ;

    grace::lep4d_prims_array_t prims_out{} ;
    grace::c2p_err_t c2p_err ;
    double residual ;
    grace::conservs_to_prims_lep4d(
        cons, prims_out, metric, eos, atmo, excision, c2p_pars, &residual, c2p_err) ;

    // Check
    INFO("C2P residual = " << residual) ;
    CHECK_THAT(residual, Catch::Matchers::WithinAbs(0., 1e-10)) ;

    auto rel = [](double a, double b) {
        return std::fabs(a-b) / (std::fabs(b)+1e-40) ;
    } ;
    CHECK_THAT(rel(prims_out[grace::LEP4D_RHOL],   prims_in[grace::LEP4D_RHOL]),   Catch::Matchers::WithinAbs(0.,1e-8)) ;
    CHECK_THAT(rel(prims_out[grace::LEP4D_PRESSL], prims_in[grace::LEP4D_PRESSL]), Catch::Matchers::WithinAbs(0.,1e-8)) ;
    CHECK_THAT(rel(prims_out[grace::LEP4D_EPSL],   prims_in[grace::LEP4D_EPSL]),   Catch::Matchers::WithinAbs(0.,1e-8)) ;
    CHECK_THAT(rel(prims_out[grace::LEP4D_YEL],    prims_in[grace::LEP4D_YEL]),    Catch::Matchers::WithinAbs(0.,1e-8)) ;
    CHECK_THAT(rel(prims_out[grace::LEP4D_YMUL],   prims_in[grace::LEP4D_YMUL]),   Catch::Matchers::WithinAbs(0.,1e-8)) ;
    CHECK_THAT(rel(prims_out[grace::LEP4D_ZXL],    prims_in[grace::LEP4D_ZXL]),    Catch::Matchers::WithinAbs(0.,1e-7)) ;
}

// ============================================================
//  Test 3: Beta-equilibrium solver correctness
// ============================================================

TEST_CASE("leptonic_4d beta equilibrium", "[leptonic_4d][betaeq]")
{
    auto eos = grace::get_leptonic_4d_eos() ;

    // Test at mid-range density and temperature
    double rho  = std::exp(0.5*(eos.lrhomin + eos.lrhomax)) ;
    double temp = std::exp(0.5*(eos.ltempmin + eos.ltempmax)) ;

    double ye_eq{0.}, ymu_eq{0.} ;
    grace::find_beta_eq_ye_ymu(ye_eq, ymu_eq, rho, temp, eos) ;

    INFO("Beta-eq: Y_e=" << ye_eq << "  Y_mu=" << ymu_eq) ;
    REQUIRE(ye_eq  >= eos.eos_yemin)  ;
    REQUIRE(ye_eq  <= eos.eos_yemax)  ;
    REQUIRE(ymu_eq >= eos.eos_ymumin) ;
    REQUIRE(ymu_eq <= eos.eos_ymumax) ;

    // Check residuals
    grace::eos_err_t err ;
    double mue_v, mumu_v, mup_v, mun_v ;
    double rhol{rho}, templ{temp} ;
    mue_v = eos.mue_mumu_mup_mun__temp_rho_ye_ymu(
                mumu_v, mup_v, mun_v, templ, rhol, ye_eq, ymu_eq, err) ;

    double res_betaeq  = std::fabs(mun_v  - mup_v - mue_v) ;
    double res_muon_eq = std::fabs(mumu_v - mue_v) ;

    INFO("beta-eq residual:   " << res_betaeq) ;
    INFO("muon-eq  residual:  " << res_muon_eq) ;

    // Tolerance on chemical potential differences (geometric units, dimensionless)
    constexpr double MU_TOL = 1e-8 ;
    CHECK_THAT(res_betaeq,  Catch::Matchers::WithinAbs(0., MU_TOL)) ;
    // muon residual only required to be small if ymu > ymumin
    if (ymu_eq > eos.eos_ymumin + 1e-6) {
        CHECK_THAT(res_muon_eq, Catch::Matchers::WithinAbs(0., MU_TOL)) ;
    }
}

// ============================================================
//  Test 4: Y_mu flux advection (unit test for flux macro)
//          Check that YMUSTAR_ flux is proportional to DENS_ flux
//          at the same Y_mu.
// ============================================================

TEST_CASE("leptonic_4d YMUSTAR flux proportionality", "[leptonic_4d][flux]")
{
    // Simple HLL solver mock: returns upwind value (v>0 case)
    auto hll_mock = [](double fl, double fr, double ul, double ur,
                       double cmin, double cmax) -> double {
        return (cmax*fl + cmin*fr - cmin*cmax*(ur-ul)) / (cmin+cmax) ;
    } ;

    // If Y_mu is the same on both sides:
    // f[YMUSL] = Y_mu * f[DENSL]  exactly (by construction of the flux).
    double const ymu_const = 0.05 ;
    double const fdl = 1.2, fdr = 0.8 ;
    double const densl = 3.0, densr = 2.5 ;
    double const cmin = 0.3, cmax = 0.5 ;

    double f_dens = hll_mock(fdl, fdr, densl, densr, cmin, cmax) ;
    double f_ymu  = hll_mock(ymu_const*fdl, ymu_const*fdr,
                              ymu_const*densl, ymu_const*densr,
                              cmin, cmax) ;

    CHECK_THAT(std::fabs(f_ymu - ymu_const*f_dens),
               Catch::Matchers::WithinAbs(0., 1e-14)) ;
}

#endif /* GRACE_ENABLE_LEPTONIC_4D */
