/**
 * @file kastaun_c2p_lep4d.hh
 * @brief  Conservative-to-primitive inversion for GRMHD with a 4D leptonic EOS
 *         (rho, T, Y_e, Y_mu).  Extends the existing kastaun_c2p.hh to handle
 *         the muon fraction as an additional evolved conserved variable YMUSTAR_.
 *
 * Design:
 *  - The muon fraction Y_mu is advected like Y_e (YMUSTAR_ / DENS_).
 *  - The Kastaun root-finding residual is identical to the standard MHD one
 *    (it only depends on the total EOS h); the 4D EOS is called with both
 *    Y_e and Y_mu passed through explicitly.
 *  - A companion entropy-backup c2p is also provided.
 *
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifndef GRACE_C2P_KASTAUN_LEP4D_HH
#define GRACE_C2P_KASTAUN_LEP4D_HH

#include <grace_config.h>
#include <grace/utils/device.h>
#include <grace/utils/metric_utils.hh>
#include <grace/utils/rootfinding.hh>
#include <grace/physics/eos/eos_base.hh>
#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/grmhd_helpers.hh>
#include <grace/physics/eos/c2p.hh>

#include <Kokkos_Core.hpp>

namespace grace {

// ===========================================================================
//  Extended conservative / primitive index enums for the 4D leptonic sector.
//  These **extend** the existing GRMHD_CONS_LOC_INDICES / GRMHD_PRIMS_LOC_INDICES
//  without conflicting with them.
// ===========================================================================

/**
 * @brief Conservative array indices for the leptonic-4D sector.
 *        The first N entries mirror GRMHD_CONS_LOC_INDICES exactly so that
 *        the same grmhd_cons_array_t functions can be used for common ops.
 */
enum LEP4D_CONS_LOC_INDICES : int {
    LEP4D_DENSL  = 0,
    LEP4D_STXL,
    LEP4D_STYL,
    LEP4D_STZL,
    LEP4D_TAUL,
    LEP4D_YESL,       ///< Y_e * sqrt(g) * D
    LEP4D_ENTSL,      ///< entropy * sqrt(g) * D
    LEP4D_BSXL,
    LEP4D_BSYL,
    LEP4D_BSZL,
    LEP4D_YMUSL,      ///< Y_mu * sqrt(g) * D   (new)
    LEP4D_NUM_CONS_LOC
} ;

/**
 * @brief Primitive array indices for the leptonic-4D sector.
 */
enum LEP4D_PRIMS_LOC_INDICES : int {
    LEP4D_RHOL  = 0,
    LEP4D_PRESSL,
    LEP4D_ZXL,
    LEP4D_ZYL,
    LEP4D_ZZL,
    LEP4D_YEL,
    LEP4D_YMUL,       ///< muon fraction Y_mu (new)
    LEP4D_TEMPL,
    LEP4D_EPSL,
    LEP4D_ENTL,
    LEP4D_BXL,
    LEP4D_BYL,
    LEP4D_BZL,
    LEP4D_NUM_PRIMS_LOC
} ;

using lep4d_cons_array_t  = std::array<double, LEP4D_NUM_CONS_LOC>  ;
using lep4d_prims_array_t = std::array<double, LEP4D_NUM_PRIMS_LOC> ;

// ===========================================================================
//  fbrack_lep4d_t  –  identical to fbrack_t but self-contained.
//  (Bracket search for the initial mu0 when r^2 >= h_min^2.)
// ===========================================================================
struct fbrack_lep4d_t {
    KOKKOS_FUNCTION fbrack_lep4d_t(double _rsqr, double _bsqr,
                                   double _rbsqr, double _h0)
        : rsqr(_rsqr), bsqr(_bsqr), rbsqr(_rbsqr), h0sqr(_h0*_h0) {}

    KOKKOS_INLINE_FUNCTION double x__mu(double mu)  const { return 1./(1.+mu*bsqr) ; }
    KOKKOS_INLINE_FUNCTION double rfsqr__mu(double mu) const {
        double x = x__mu(mu) ;
        return x*(rsqr*x + mu*(x+1.)*rbsqr) ;
    }
    KOKKOS_INLINE_FUNCTION double operator()(double mu) const {
        return mu*std::sqrt(h0sqr + rfsqr__mu(mu)) - 1. ;
    }
    KOKKOS_INLINE_FUNCTION
    std::pair<double,double> operator_with_deriv(double mu) const {
        double rf2 = rfsqr__mu(mu) ;
        double sq  = std::sqrt(h0sqr + rf2) ;
        return { mu*sq - 1., sq } ; // simplified derivative
    }
    KOKKOS_INLINE_FUNCTION void bracket(double& mua, double& mub) const {
        mua = 0. ;
        mub = 1./std::sqrt(h0sqr) ;
    }
    double rsqr, bsqr, rbsqr, h0sqr ;
} ;

// ===========================================================================
//  froot_lep4d_t  –  Kastaun root functor for 4D leptonic EOS.
//  Mirrors froot_t in kastaun_c2p.hh but calls the 4D EOS interface.
// ===========================================================================
template< typename eos_t >
struct froot_lep4d_t {

    KOKKOS_FUNCTION
    froot_lep4d_t(eos_t _eos, double _d, double _qtot,
                  double _r2, double _rbsqr, double _bsqr, double _brosqr,
                  double _ye, double _ymu, double _h0)
        : eos(_eos), d(_d), qtot(_qtot), rsqr(_r2),
          rbsqr(_rbsqr), bsqr(_bsqr), brosqr(_brosqr),
          ye(_ye), ymu(_ymu), h0sqr(_h0*_h0), vsqrmax(1.-1./(50.*50.))
    { wmax = 1./std::sqrt(1.-vsqrmax) ; }

    KOKKOS_INLINE_FUNCTION double x__mu(double mu)     const { return 1./(1.+mu*bsqr) ; }
    KOKKOS_INLINE_FUNCTION double rfsqr__mu_x(double mu, double x) const {
        return x*(rsqr*x + mu*(x+1.)*rbsqr) ;
    }
    KOKKOS_INLINE_FUNCTION double qf__mu_x(double mu, double x) const {
        double mux = mu*x ;
        return qtot - 0.5*(bsqr + mux*mux*brosqr) ;
    }
    KOKKOS_INLINE_FUNCTION double
    eps_raw(double mu, double qf, double rfsqr_, double w) const {
        return w*(qf - mu*rfsqr_*(1. - mu*w/(1.+w))) ;
    }

    KOKKOS_INLINE_FUNCTION void
    get_eps_range(double& emin, double& emax, double rho) const {
        double yel{ye}, ymul{ymu} ;
        eos_err_t eoserr ;
        double rhol{rho} ;
        eos.eps_range__rho_ye_ymu(emin, emax, rhol, yel, ymul, eoserr) ;
    }

    KOKKOS_INLINE_FUNCTION double operator()(double mu) const {
        double x      = x__mu(mu) ;
        double rfsqr_ = rfsqr__mu_x(mu,x) ;
        double qf     = qf__mu_x(mu,x) ;
        double vsqr   = rfsqr_*mu*mu ;
        double w      = (vsqr > vsqrmax) ? wmax : 1./std::sqrt(1.-vsqr) ;

        double rhomax = eos.density_maximum() ;
        double rhomin = eos.density_minimum() ;
        double rho    = Kokkos::fmin(rhomax, Kokkos::fmax(rhomin, d/w)) ;

        double ep  = eps_raw(mu, qf, rfsqr_, w) ;
        double emin, emax ;
        get_eps_range(emin, emax, rho) ;
        ep = Kokkos::fmin(emax, Kokkos::fmax(emin, ep)) ;

        double yel{ye}, ymul{ymu} ;
        eos_err_t eoserr ;
        double press = eos.press__eps_rho_ye_ymu(ep, rho, yel, ymul, eoserr) ;

        double a   = press / (rho*(1.+ep)) ;
        double h   = (1.+ep)*(1.+a) ;
        double hbw = Kokkos::fmax((1.+a)*(1.+qf-mu*rfsqr_), h/w) ;
        return mu - 1./(hbw + rfsqr_*mu) ;
    }

    /**
     * @brief Compute all primitives given the root mu.
     *        Errors are accumulated in the c2p signal bitset.
     */
    KOKKOS_INLINE_FUNCTION double
    compute_primitives(double mu, c2p_sig_t& err,
                       double& x, double& w, double& rho,
                       double& press, double& eps,
                       double& temp, double& entropy) const
    {
        x           = x__mu(mu) ;
        double rfsqr_ = rfsqr__mu_x(mu,x) ;
        double qf     = qf__mu_x(mu,x) ;
        double vsqr   = rfsqr_*mu*mu ;

        if (vsqr > vsqrmax) {
            vsqr = vsqrmax ;
            w    = wmax    ;
            err.set(c2p_sig_enum_t::C2P_VEL_TOO_HIGH) ;
        } else {
            w = 1./std::sqrt(1.-vsqr) ;
        }

        double rhomax = eos.density_maximum() ;
        double rhomin = eos.density_minimum() ;
        rho = d/w ;
        if (rho >= rhomax) { err.set(c2p_sig_enum_t::C2P_RHO_TOO_HIGH) ; rho = rhomax ; }
        else if (rho <= rhomin) { err.set(c2p_sig_enum_t::C2P_RHO_TOO_LOW) ; rho = rhomin ; }

        eps = eps_raw(mu, qf, rfsqr_, w) ;
        double emin, emax ;
        get_eps_range(emin, emax, rho) ;
        emax = Kokkos::fmin(emax, eos.get_c2p_eps_max()) ;
        if (eps >= emax) { err.set(c2p_sig_enum_t::C2P_EPS_TOO_HIGH) ; eps = emax*0.999 ; }
        else if (eps < emin) { err.set(c2p_sig_enum_t::C2P_EPS_TOO_LOW) ; eps = emin ; }

        double hh, csnd2 ;
        eos_err_t eoserr ;
        double yel{ye}, ymul{ymu} ;
        press = eos.press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
            hh, csnd2, temp, entropy, eps, rho, yel, ymul, eoserr) ;

        double a      = press/(rho*(1.+eps)) ;
        double h      = (1.+eps)*(1.+a) ;
        double hbw    = Kokkos::fmax((1.+a)*(1.+qf-mu*rfsqr_), h/w) ;
        double newmu  = 1./(hbw + rfsqr_*mu) ;
        return mu - newmu ;
    }

    eos_t eos ;
    double d, qtot, ye, ymu ;
    double rsqr, bsqr, rbsqr, brosqr, h0sqr ;
    double vsqrmax, wmax ;
} ;

// ===========================================================================
//  kastaun_c2p_lep4d_t  –  main c2p struct for 4D leptonic EOS.
//  Drop-in replacement for kastaun_c2p_t but operating on lep4d arrays.
// ===========================================================================
template< typename eos_t >
struct kastaun_c2p_lep4d_t {

    GRACE_HOST_DEVICE
    kastaun_c2p_lep4d_t(
        eos_t const&            _eos,
        metric_array_t const&   _metric,
        lep4d_cons_array_t&     conservs
    ) : eos(_eos), metric(_metric), h0(_eos.enthalpy_minimum())
    {
        double const B2 = metric.square_vec({conservs[LEP4D_BSXL],
                                             conservs[LEP4D_BSYL],
                                             conservs[LEP4D_BSZL]}) ;
        conservs[LEP4D_DENSL] = Kokkos::fmax(0., conservs[LEP4D_DENSL]) ;
        D = conservs[LEP4D_DENSL] ;

        q  = conservs[LEP4D_TAUL] / D ;
        r  = { conservs[LEP4D_STXL]/D,
               conservs[LEP4D_STYL]/D,
               conservs[LEP4D_STZL]/D } ;

        Btilde = { conservs[LEP4D_BSXL]/std::sqrt(D),
                   conservs[LEP4D_BSYL]/std::sqrt(D),
                   conservs[LEP4D_BSZL]/std::sqrt(D) } ;
        B      = { conservs[LEP4D_BSXL],
                   conservs[LEP4D_BSYL],
                   conservs[LEP4D_BSZL] } ;

        r2    = metric.square_covec(r) ;
        Btilde2 = metric.square_vec(Btilde) ;
        r_dot_Btilde  = r[0]*Btilde[0]+r[1]*Btilde[1]+r[2]*Btilde[2] ;
        r_dot_Btilde2 = r_dot_Btilde*r_dot_Btilde ;

        r = metric.raise(r) ;

        ye  = conservs[LEP4D_YESL]  / D ;
        ymu = conservs[LEP4D_YMUSL] / D ;

        // Clamp fractions to table range
        ye  = Kokkos::fmin(eos.eos_yemax,  Kokkos::fmax(eos.eos_yemin,  ye))  ;
        ymu = Kokkos::fmin(eos.eos_ymumax, Kokkos::fmax(eos.eos_ymumin, ymu)) ;

        v02 = r2 / (h0*h0 + r2) ;
    }

    /**
     * @brief Invert conservatives → primitives.
     * @param prims  Output primitives (lep4d_prims_array_t).
     * @param c2p_errors  Error bitset (in/out).
     * @return Residual |f(mu_final)| of the root equation.
     */
    double GRACE_HOST_DEVICE
    invert(lep4d_prims_array_t& prims, c2p_sig_t& c2p_errors)
    {
        constexpr double tolerance = 1e-15 ;

        // --- bracket mu0 ---
        double mu0 = 1./h0 ;
        if (r2 >= h0*h0) {
            fbrack_lep4d_t g(r2, Btilde2, r_dot_Btilde2, h0) ;
            double mua, mub ;
            g.bracket(mua, mub) ;
            int rerr ;
            mu0 = utils::rootfind_newton_raphson(mua, mub, g, 30, 1e-10, rerr) ;
            if (rerr != 0) mu0 = 1./h0 ;
            else           mu0 *= 1.+1e-10 ;
        }

        // total energy density helper
        double qtot = q + 0.5*(Btilde2/D) ;  // q̃ as in Kastaun (2021)
        double brosqr = r_dot_Btilde2 ;
        froot_lep4d_t<eos_t> fmu(eos, D, qtot, r2, r_dot_Btilde2, Btilde2,
                                   brosqr, ye, ymu, h0) ;
        double mu = utils::brent(fmu, 0., mu0, tolerance) ;

        // recover primitives
        double x, w, eps, rho, press, temp, entropy ;
        double residual = fmu.compute_primitives(
            mu, c2p_errors, x, w, rho, press, eps, temp, entropy) ;

        prims[LEP4D_RHOL]   = rho     ;
        prims[LEP4D_PRESSL] = press   ;
        prims[LEP4D_EPSL]   = eps     ;
        prims[LEP4D_TEMPL]  = temp    ;
        prims[LEP4D_ENTL]   = entropy ;
        prims[LEP4D_YEL]    = ye      ;
        prims[LEP4D_YMUL]   = ymu     ;
        prims[LEP4D_BXL]    = B[0]    ;
        prims[LEP4D_BYL]    = B[1]    ;
        prims[LEP4D_BZL]    = B[2]    ;

        // z-vector (undensitized momentum per unit enthalpy density)
        for (int ii=0; ii<3; ++ii)
            prims[LEP4D_ZXL+ii] = w * mu * x * (r[ii] + mu*r_dot_Btilde*Btilde[ii]) ;

        return std::fabs(residual) ;
    }

  private:
    eos_t eos ;
    metric_array_t metric ;
    double r2, q, Btilde2, D, ye, ymu, r_dot_Btilde, r_dot_Btilde2, h0, v02 ;
    std::array<double,3> r, Btilde, B ;
} ;

// ===========================================================================
//  Entropy-backup c2p for 4D leptonic EOS
//  (mirrors ent_froot_t / entropy_fix_c2p_t in ent_based_c2p.hh)
// ===========================================================================
template< typename eos_t >
struct ent_froot_lep4d_t {

    KOKKOS_FUNCTION
    ent_froot_lep4d_t(eos_t _eos, double _d,
                      double _rsqr, double _rbsqr, double _bsqr,
                      double _s, double _ye, double _ymu, double _h0)
        : eos(_eos), d(_d), s(_s), ye(_ye), ymu(_ymu),
          rsqr(_rsqr), bsqr(_bsqr), rbsqr(_rbsqr)
    {
        double zsqrmax = rsqr/(_h0*_h0) ;
        wmax = std::sqrt(1. + zsqrmax) ;
        vsqrmax = zsqrmax / (1. + zsqrmax) ;
    }

    KOKKOS_INLINE_FUNCTION double x__mu(double mu)     const { return 1./(1.+mu*bsqr) ; }
    KOKKOS_INLINE_FUNCTION double rfsqr__mu_x(double mu, double x) const {
        return x*(rsqr*x + mu*(x+1.)*rbsqr) ;
    }

    KOKKOS_INLINE_FUNCTION double operator()(double mu) const {
        double x      = x__mu(mu) ;
        double rfsqr_ = rfsqr__mu_x(mu,x) ;
        double vsqr   = rfsqr_*mu*mu ;
        double w      = (vsqr > vsqrmax) ? wmax : 1./std::sqrt(1.-vsqr) ;

        double rhomax = eos.density_maximum() ;
        double rhomin = eos.density_minimum() ;
        double rho    = Kokkos::fmin(rhomax, Kokkos::fmax(rhomin, d/w)) ;

        double hh, csnd2, temp, eps ;
        eos_err_t eoserr ;
        double yel{ye}, ymul{ymu}, sl{s} ;
        eos.press_h_csnd2_temp_eps__entropy_rho_ye_ymu(
            hh, csnd2, temp, eps, sl, rho, yel, ymul, eoserr) ;

        double a   = eos.press__eps_rho_ye_ymu(eps, rho, yel, ymul, eoserr) / (rho*(1.+eps)) ;
        double h   = (1.+eps)*(1.+a) ;
        double hbw = h/w ;
        return mu - 1./(hbw + rfsqr_*mu) ;
    }

    KOKKOS_INLINE_FUNCTION double
    compute_primitives(double mu, c2p_sig_t& err,
                       double& x, double& w, double& rho,
                       double& press, double& eps, double& temp) const
    {
        x           = x__mu(mu) ;
        double rfsqr_ = rfsqr__mu_x(mu,x) ;
        double vsqr   = rfsqr_*mu*mu ;

        if (vsqr > vsqrmax) {
            vsqr = vsqrmax ;
            w    = wmax    ;
            err.set(c2p_sig_enum_t::C2P_VEL_TOO_HIGH) ;
        } else {
            w = 1./std::sqrt(1.-vsqr) ;
        }

        double rhomax = eos.density_maximum() ;
        double rhomin = eos.density_minimum() ;
        rho = d/w ;
        if (rho >= rhomax) { err.set(c2p_sig_enum_t::C2P_RHO_TOO_HIGH) ; rho = rhomax ; }
        else if (rho <= rhomin) { err.set(c2p_sig_enum_t::C2P_RHO_TOO_LOW) ; rho = rhomin ; }

        double hh, csnd2 ;
        eos_err_t eoserr ;
        double yel{ye}, ymul{ymu}, sl{s} ;
        press = eos.press_h_csnd2_temp_eps__entropy_rho_ye_ymu(
            hh, csnd2, temp, eps, sl, rho, yel, ymul, eoserr) ;

        if (eoserr.test(EOS_ERROR_T::EOS_ENTROPY_TOO_LOW))  err.set(C2P_ENT_TOO_LOW)  ;
        if (eoserr.test(EOS_ERROR_T::EOS_ENTROPY_TOO_HIGH)) err.set(C2P_ENT_TOO_HIGH) ;

        double a   = press/(rho*(1.+eps)) ;
        double h   = (1.+eps)*(1.+a) ;
        double hbw = h/w ;
        return mu - 1./(hbw + rfsqr_*mu) ;
    }

    eos_t eos ;
    double d, s, ye, ymu ;
    double rsqr, bsqr, rbsqr ;
    double vsqrmax, wmax ;
} ;

template< typename eos_t >
struct entropy_fix_c2p_lep4d_t {

    GRACE_HOST_DEVICE
    entropy_fix_c2p_lep4d_t(
        eos_t const&          _eos,
        metric_array_t const& _metric,
        lep4d_cons_array_t&   conservs
    ) : eos(_eos), metric(_metric), h0(_eos.enthalpy_minimum())
    {
        conservs[LEP4D_DENSL] = Kokkos::fmax(0., conservs[LEP4D_DENSL]) ;
        D  = conservs[LEP4D_DENSL] ;
        s  = conservs[LEP4D_ENTSL] / D ;
        r  = { conservs[LEP4D_STXL]/D,
               conservs[LEP4D_STYL]/D,
               conservs[LEP4D_STZL]/D } ;
        Btilde = { conservs[LEP4D_BSXL]/std::sqrt(D),
                   conservs[LEP4D_BSYL]/std::sqrt(D),
                   conservs[LEP4D_BSZL]/std::sqrt(D) } ;
        B      = { conservs[LEP4D_BSXL],
                   conservs[LEP4D_BSYL],
                   conservs[LEP4D_BSZL] } ;
        r2    = metric.square_covec(r) ;
        Btilde2 = metric.square_vec(Btilde) ;
        r_dot_Btilde  = r[0]*Btilde[0]+r[1]*Btilde[1]+r[2]*Btilde[2] ;
        r_dot_Btilde2 = r_dot_Btilde*r_dot_Btilde ;
        r = metric.raise(r) ;
        ye  = Kokkos::fmin(eos.eos_yemax,  Kokkos::fmax(eos.eos_yemin,  conservs[LEP4D_YESL]/D)) ;
        ymu = Kokkos::fmin(eos.eos_ymumax, Kokkos::fmax(eos.eos_ymumin, conservs[LEP4D_YMUSL]/D)) ;
        v02 = r2/(h0*h0 + r2) ;
    }

    double GRACE_HOST_DEVICE
    invert(lep4d_prims_array_t& prims, c2p_sig_t& c2p_errors)
    {
        constexpr double tolerance = 1e-15 ;
        double mu0 = 1./h0 ;
        if (r2 >= h0*h0) {
            fbrack_lep4d_t g(r2, Btilde2, r_dot_Btilde2, h0) ;
            double mua, mub ;
            g.bracket(mua, mub) ;
            int rerr ;
            mu0 = utils::rootfind_newton_raphson(mua, mub, g, 30, 1e-10, rerr) ;
            if (rerr != 0) mu0 = 1./h0 ; else mu0 *= 1.+1e-10 ;
        }

        ent_froot_lep4d_t<eos_t> fmu(eos, D, r2, r_dot_Btilde2, Btilde2,
                                      s, ye, ymu, h0) ;
        double mu = utils::brent(fmu, 0., mu0, tolerance) ;

        double x, w, eps, rho, press, temp ;
        double residual = fmu.compute_primitives(
            mu, c2p_errors, x, w, rho, press, eps, temp) ;

        prims[LEP4D_RHOL]   = rho   ;
        prims[LEP4D_PRESSL] = press ;
        prims[LEP4D_EPSL]   = eps   ;
        prims[LEP4D_TEMPL]  = temp  ;
        prims[LEP4D_ENTL]   = s     ;
        prims[LEP4D_YEL]    = ye    ;
        prims[LEP4D_YMUL]   = ymu   ;
        prims[LEP4D_BXL]    = B[0]  ;
        prims[LEP4D_BYL]    = B[1]  ;
        prims[LEP4D_BZL]    = B[2]  ;
        for (int ii=0; ii<3; ++ii)
            prims[LEP4D_ZXL+ii] = w*mu*x*(r[ii]+mu*r_dot_Btilde*Btilde[ii]) ;

        return std::fabs(residual) ;
    }

  private:
    eos_t eos ;
    metric_array_t metric ;
    double r2, s, Btilde2, D, ye, ymu, r_dot_Btilde, r_dot_Btilde2, h0, v02 ;
    std::array<double,3> r, Btilde, B ;
} ;

} /* namespace grace */
#endif /* GRACE_C2P_KASTAUN_LEP4D_HH */
