/**
 * @file grmhd_lep4d_flux_patch.hh
 * @brief  Helper macros and inline functions for computing GRMHD fluxes
 *         in the 4D leptonic EOS mode (GRACE_ENABLE_LEPTONIC_4D).
 *
 *         This patch supplies the additional Y_mu flux term that must be
 *         added alongside the existing Y_e flux in grmhd.hh whenever
 *         the leptonic-4D EOS is active.
 *
 *         Usage pattern in the flux kernel (grmhd.hh):
 *
 *           f[YESL]   = sqrtg * solver(yel*fdl,  yer*fdr,  yel*densl,  yer*densr,  cmin,cmax) ;
 *       #ifdef GRACE_ENABLE_LEPTONIC_4D
 *           f[YMUSL]  = sqrtg * solver(ymul*fdl, ymur*fdr, ymul*densl, ymur*densr, cmin,cmax) ;
 *       #endif
 *
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifndef GRACE_GRMHD_LEP4D_FLUX_PATCH_HH
#define GRACE_GRMHD_LEP4D_FLUX_PATCH_HH

#ifdef GRACE_ENABLE_LEPTONIC_4D

#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/eos/c2p_lep4d.hh>
#include <grace/physics/grmhd_helpers.hh>

namespace grace {

// ---------------------------------------------------------------------------
//  compute_mhd_fluxes_lep4d
//
//  Computes the conserved-variable flux vector for the leptonic-4D sector.
//  All quantities are in the local (undensitized, cell-volume) frame.
//  The caller is responsible for multiplying by sqrtg.
//
//  Returns the HLL flux array f (lep4d_cons_array_t layout).
// ---------------------------------------------------------------------------

/**
 * @brief Compute MHD flux for the 4D leptonic EOS sector.
 *
 * @tparam idir   Flux direction (0/1/2 = x/y/z).
 * @tparam eos_t  EOS type (leptonic_eos_4d_t or compatible).
 * @param primL   Left  state (lep4d_prims_array_t).
 * @param primR   Right state (lep4d_prims_array_t).
 * @param metric  Face metric.
 * @param f       Output flux vector (lep4d_cons_array_t).
 * @param vbar    Signal speed array [4] = {vl, vr, cp, cm}.
 * @param eos     EOS object.
 * @param cmin    Minimum signal speed (input or computed internally).
 * @param cmax    Maximum signal speed.
 */
template< size_t idir, typename eos_t, typename riemann_solver_t >
GRACE_HOST_DEVICE GRACE_ALWAYS_INLINE
void compute_mhd_fluxes_lep4d(
    lep4d_prims_array_t&   primL,
    lep4d_prims_array_t&   primR,
    metric_array_t const&  metric,
    lep4d_cons_array_t&    f,
    std::array<double,4>&  vbar,
    eos_t const&           eos,
    riemann_solver_t const& solver,
    double const cmin = 1.,
    double const cmax = 1. )
{
    // ----------------------------------------------------------------
    //  Pointers into primitive arrays
    // ----------------------------------------------------------------
    double const * const zl    = &(primL[LEP4D_ZXL])  ;
    double const * const zr    = &(primR[LEP4D_ZXL])  ;
    double const * const Bl    = &(primL[LEP4D_BXL])  ;
    double const * const Br    = &(primR[LEP4D_BXL])  ;
    double& rhol  = primL[LEP4D_RHOL]  ;
    double& rhor  = primR[LEP4D_RHOL]  ;
    double& sl    = primL[LEP4D_ENTL]  ;
    double& sr    = primR[LEP4D_ENTL]  ;
    double& tl    = primL[LEP4D_TEMPL] ;
    double& tr    = primR[LEP4D_TEMPL] ;
    double& yel   = primL[LEP4D_YEL]   ;
    double& yer   = primR[LEP4D_YEL]   ;
    double& ymul  = primL[LEP4D_YMUL]  ;
    double& ymur  = primR[LEP4D_YMUL]  ;

    // ----------------------------------------------------------------
    //  Lorentz factors
    // ----------------------------------------------------------------
    const double * const gdd = metric._g.data() ;
    double wl, wr ;
    grmhd_get_W(gdd, zl, &wl) ;
    grmhd_get_W(gdd, zr, &wr) ;

    // ----------------------------------------------------------------
    //  Pressure / eps / cs2  (4D EOS call)
    // ----------------------------------------------------------------
    eos_err_t eoserr ;
    double epsl, epsr, pl, pr, cs2l, cs2r ;
    pl = eos.press_eps_csnd2__temp_rho_ye_ymu(epsl, cs2l, tl, rhol, yel, ymul, eoserr) ;
    pr = eos.press_eps_csnd2__temp_rho_ye_ymu(epsr, cs2r, tr, rhor, yer, ymur, eoserr) ;

    // ----------------------------------------------------------------
    //  b^mu, b^2 on both sides
    // ----------------------------------------------------------------
    const double * const betau = metric._beta.data() ;
    double alp = metric.alp() ;
    double smallbl[4], smallbr[4] ;
    double b2l, b2r ;
    grmhd_get_smallbu_smallb2(betau, gdd, Bl, zl, wl, alp, &smallbl, &b2l) ;
    grmhd_get_smallbu_smallb2(betau, gdd, Br, zr, wr, alp, &smallbr, &b2r) ;

    double vtildel[3], vtilder[3] ;
    grmhd_get_vtildeu(betau, wl, zl, alp, &vtildel) ;
    grmhd_get_vtildeu(betau, wr, zr, alp, &vtilder) ;

    // ----------------------------------------------------------------
    //  Conserved variables and fluxes on each side
    // ----------------------------------------------------------------
    double Dl, Dr, taul, taur ;
    double stl[3], str[3] ;
    double fdl, fdr, ftl, ftr ;
    double fstl[3], fstr[3] ;
    double entsl, entsr, fenl, fenr ;

    double hhl = 1. + epsl + pl/rhol ;
    double hhr = 1. + epsr + pr/rhor ;

    grmhd_get_conserved_and_flux<idir>(
        wl, rhol, smallbl, b2l, alp, epsl, pl, betau, zl, gdd, sl,
        &Dl, &taul, &stl, &entsl, &fdl, &ftl, &fstl, &fenl) ;
    grmhd_get_conserved_and_flux<idir>(
        wr, rhor, smallbr, b2r, alp, epsr, pr, betau, zr, gdd, sr,
        &Dr, &taur, &str, &entsr, &fdr, &ftr, &fstr, &fenr) ;

    double densl = Dl, densr = Dr ;

    // ----------------------------------------------------------------
    //  Signal speeds (using existing GRACE utility)
    // ----------------------------------------------------------------
    // (cmin/cmax already supplied by the caller from the standard
    //  compute_mhd_fluxes path or computed here)

    // ----------------------------------------------------------------
    //  HLL / HLL-like Riemann solve for each conserved variable
    // ----------------------------------------------------------------
    f[LEP4D_DENSL] = metric.sqrtg() * solver(fdl, fdr, densl, densr, cmin, cmax) ;
    f[LEP4D_ENTSL] = metric.sqrtg() * solver(fenl, fenr, entsl, entsr, cmin, cmax) ;
    f[LEP4D_TAUL]  = metric.sqrtg() * solver(ftl, ftr, taul, taur, cmin, cmax) ;
    f[LEP4D_STXL]  = metric.sqrtg() * solver(fstl[0],fstr[0],stl[0],str[0],cmin,cmax) ;
    f[LEP4D_STYL]  = metric.sqrtg() * solver(fstl[1],fstr[1],stl[1],str[1],cmin,cmax) ;
    f[LEP4D_STZL]  = metric.sqrtg() * solver(fstl[2],fstr[2],stl[2],str[2],cmin,cmax) ;
    // Advected fractions: flux = (D-flux) * Y
    f[LEP4D_YESL]  = metric.sqrtg() * solver(yel*fdl, yer*fdr, yel*densl, yer*densr, cmin, cmax) ;
    f[LEP4D_YMUSL] = metric.sqrtg() * solver(ymul*fdl, ymur*fdr, ymul*densl, ymur*densr, cmin, cmax) ;
}

// ---------------------------------------------------------------------------
//  Scatter 4D leptonic flux components back to the standard evolved
//  state array (maps lep4d_cons_array_t → state view indices).
// ---------------------------------------------------------------------------
#define SCATTER_LEP4D_FLUX_TO_STATE(fview, flep4d, idir, q, i, j, k)    \
fview(VEC(i,j,k), DENS_,        idir, q) = flep4d[LEP4D_DENSL]  ;       \
fview(VEC(i,j,k), TAU_,         idir, q) = flep4d[LEP4D_TAUL]   ;       \
fview(VEC(i,j,k), SX_,          idir, q) = flep4d[LEP4D_STXL]   ;       \
fview(VEC(i,j,k), SY_,          idir, q) = flep4d[LEP4D_STYL]   ;       \
fview(VEC(i,j,k), SZ_,          idir, q) = flep4d[LEP4D_STZL]   ;       \
fview(VEC(i,j,k), YESTAR_,      idir, q) = flep4d[LEP4D_YESL]   ;       \
fview(VEC(i,j,k), ENTROPYSTAR_, idir, q) = flep4d[LEP4D_ENTSL]  ;       \
fview(VEC(i,j,k), YMUSTAR_,     idir, q) = flep4d[LEP4D_YMUSL]

} /* namespace grace */

#endif /* GRACE_ENABLE_LEPTONIC_4D */
#endif /* GRACE_GRMHD_LEP4D_FLUX_PATCH_HH */
