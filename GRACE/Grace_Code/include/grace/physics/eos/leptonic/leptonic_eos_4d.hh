/**
 * @file leptonic_eos_4d.hh
 * @author Based on Margherita leptonic_eos.hh by Harry Ho-Yin Ng and Elias Roland Most
 * @brief  4D tabulated EOS with muon fraction Y_mu as a fourth independent variable,
 *         consistent with the GRACE EOS framework (CRTP / eos_base_t) and Kokkos.
 * @date   2025
 *
 * @copyright This file is part of the General Relativistic Astrophysics
 * Code for Exascale (GRACE).
 * GRACE is an evolution framework that uses Finite Volume methods to
 * simulate relativistic spacetimes and plasmas.
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
 */

#ifndef GRACE_PHYSICS_EOS_LEPTONIC_4D_HH
#define GRACE_PHYSICS_EOS_LEPTONIC_4D_HH

#include <grace_config.h>
#include <grace/utils/grace_utils.hh>
#include <grace/utils/bitset.hh>
#include <grace/utils/rootfinding.hh>
#include <grace/physics/eos/eos_base.hh>

#include <Kokkos_Core.hpp>

#include <array>
#include <cmath>
#include <limits>

namespace grace {

// ============================================================
//  Forward declarations
// ============================================================
class leptonic_eos_4d_t ;

// ============================================================
//  4-D table interpolator
//  Axes: log(rho), log(T), Y_e, Y_mu  (uniform spacing assumed
//  for every axis, identical layout to tabeos_linterp_t)
// ============================================================
struct tabeos_linterp_4d_t {

    tabeos_linterp_4d_t() = default ;

    tabeos_linterp_4d_t(
        Kokkos::View<double *****> tabs,   ///< [irho,iT,iye,iymu,ivar]
        Kokkos::View<double *>     ar,     ///< log(rho) axis
        Kokkos::View<double *>     at,     ///< log(T)   axis
        Kokkos::View<double *>     ay,     ///< Y_e      axis
        Kokkos::View<double *>     aym     ///< Y_mu     axis
    ) : _tables(tabs), _logrho(ar), _logT(at), _ye(ay), _ymu(aym)
    {
        idr  = 1./(_logrho(1) - _logrho(0)) ;
        idt  = 1./(_logT(1)   - _logT(0))   ;
        idy  = 1./(_ye(1)     - _ye(0))      ;
        idym = 1./(_ymu(1)    - _ymu(0))     ;
    }

    KOKKOS_INLINE_FUNCTION double lrho (int i) const { return _logrho(i) ; }
    KOKKOS_INLINE_FUNCTION double ltemp(int i) const { return _logT(i)   ; }
    KOKKOS_INLINE_FUNCTION double ye   (int i) const { return _ye(i)     ; }
    KOKKOS_INLINE_FUNCTION double ymu  (int i) const { return _ymu(i)    ; }

    KOKKOS_INLINE_FUNCTION
    double interp(double lrho, double ltemp, double ye, double ymu, int idx) const
    {
        std::array<double,1> res ;
        std::array<int,1>    _i{idx} ;
        interp<1>(lrho,ltemp,ye,ymu,_i,res) ;
        return res[0] ;
    }

    template< int N >
    KOKKOS_INLINE_FUNCTION
    void interp(double lrho, double ltemp, double ye, double ymu,
                std::array<int,N> const& idx,
                std::array<double,N>&    res) const
    {
        for( int iv=0; iv<N; ++iv) res[iv] = 0.0 ;
        int ir, it, iy, iym ;
        getidx(lrho,ltemp,ye,ymu, ir,it,iy,iym) ;
        double wr[2], wt[2], wy[2], wym[2] ;
        getw(lrho ,ir ,_logrho,idr ,wr ) ;
        getw(ltemp,it ,_logT  ,idt ,wt ) ;
        getw(ye   ,iy ,_ye    ,idy ,wy ) ;
        getw(ymu  ,iym,_ymu   ,idym,wym) ;
        for( int ii=0; ii<2; ++ii)
        for( int jj=0; jj<2; ++jj)
        for( int kk=0; kk<2; ++kk)
        for( int ll=0; ll<2; ++ll) {
            double w = wr[ii]*wt[jj]*wy[kk]*wym[ll] ;
            for( int iv=0; iv<N; ++iv)
                res[iv] += w * _tables(ir+ii, it+jj, iy+kk, iym+ll, idx[iv]) ;
        }
    }

    KOKKOS_INLINE_FUNCTION
    void getidx(double lr, double lt, double ye_, double ym_,
                int& ir, int& it, int& iy, int& iym) const
    {
        ir  = Kokkos::max(0UL,
              Kokkos::min( static_cast<size_t>((lr  - _logrho(0))*idr ),
                           _logrho.extent(0)-2 )) ;
        it  = Kokkos::max(0UL,
              Kokkos::min( static_cast<size_t>((lt  - _logT(0)  )*idt ),
                           _logT.extent(0)-2   )) ;
        iy  = Kokkos::max(0UL,
              Kokkos::min( static_cast<size_t>((ye_ - _ye(0)    )*idy ),
                           _ye.extent(0)-2     )) ;
        iym = Kokkos::max(0UL,
              Kokkos::min( static_cast<size_t>((ym_ - _ymu(0)   )*idym),
                           _ymu.extent(0)-2    )) ;
    }

    KOKKOS_INLINE_FUNCTION
    void getw(double x, int i,
              Kokkos::View<const double*> ax, double ih, double* w) const
    {
        double lam = (x - ax(i)) * ih ;
        w[0] = 1. - lam ;
        w[1] = lam      ;
    }

    Kokkos::View<const double*> _logrho, _logT, _ye, _ymu ;
    Kokkos::View<double *****>  _tables ;
    double idr, idt, idy, idym ;
} ;

// ============================================================
//  leptonic_eos_4d_t
//  Concrete EOS type: 4D tabulated (rho, T, Y_e, Y_mu).
//  Inherits the standard GRACE CRTP interface from eos_base_t.
// ============================================================

/**
 * @brief 4D leptonic tabulated EOS including muon fraction.
 * \ingroup eos
 *
 * Independent thermodynamic axes:  rho, T, Y_e, Y_mu.
 * The stored baryon table contains the same variables as
 * tabulated_eos_t::TEOS_VIDX.  Two additional per-species
 * lepton tables (ELE / MUON) carry partial contributions
 * mirroring EOS_Leptonic from Margherita.
 *
 * All methods follow the GRACE eos_base_t CRTP convention:
 *   - arguments are passed by (non-const) reference so that
 *     the implementation may clamp them to table bounds,
 *   - error codes are OR-ed into an eos_err_t bitset,
 *   - all hot-path methods are KOKKOS_INLINE_FUNCTION / GRACE_HOST_DEVICE.
 */
class leptonic_eos_4d_t
    : public eos_base_t<leptonic_eos_4d_t>
{
    using err_t  = eos_err_t ;
    using base_t = eos_base_t<leptonic_eos_4d_t> ;

  public:

    // ----------------------------------------------------------
    //  Additional EOS error bits beyond the base set
    // ----------------------------------------------------------
    enum LEP4D_ERROR_T : int {
        EOS_YMU_TOO_LOW  = EOS_NUM_ERRORS,   ///< Y_mu below table minimum
        EOS_YMU_TOO_HIGH,                     ///< Y_mu above table maximum
        EOS_4D_NUM_ERRORS
    } ;

    // ----------------------------------------------------------
    //  Table variable indices – baryon table (same as tabulated_eos_t)
    // ----------------------------------------------------------
    enum TEOS_VIDX : int {
        TABPRESS = 0,
        TABEPS,
        TABCSND2,
        TABENTROPY,
        TABMUE,
        TABMUP,
        TABMUN,
        TABXA,
        TABXH,
        TABXN,
        TABXP,
        TABABAR,
        TABZBAR,
        N_TAB_VARS_BARYON
    } ;

    // ----------------------------------------------------------
    //  Table variable indices – electron lepton table
    //  Mirrors EOS_Leptonic::EELE
    // ----------------------------------------------------------
    enum ELE_VIDX : int {
        TABMUELE = 0,      ///< electron chemical potential mu_e [geometric]
        TABYLE_MINUS,      ///< Y_{e^-}
        TABYLE_PLUS,       ///< Y_{e^+}
        TABPRESS_E_MINUS,  ///< partial pressure e^-
        TABPRESS_E_PLUS,   ///< partial pressure e^+
        TABEPS_E_MINUS,    ///< partial eps e^-
        TABEPS_E_PLUS,     ///< partial eps e^+
        TABS_E_MINUS,      ///< partial entropy e^-
        TABS_E_PLUS,       ///< partial entropy e^+
        N_TAB_VARS_ELE
    } ;

    // ----------------------------------------------------------
    //  Table variable indices – muon lepton table
    //  Mirrors EOS_Leptonic::EMUON
    // ----------------------------------------------------------
    enum MUON_VIDX : int {
        TABMUMU = 0,       ///< muon chemical potential mu_mu [geometric]
        TABYMU_MINUS,      ///< Y_{mu^-}
        TABYMU_PLUS,       ///< Y_{mu^+}
        TABPRESS_MU_MINUS, ///< partial pressure mu^-
        TABPRESS_MU_PLUS,  ///< partial pressure mu^+
        TABEPS_MU_MINUS,   ///< partial eps mu^-
        TABEPS_MU_PLUS,    ///< partial eps mu^+
        TABS_MU_MINUS,     ///< partial entropy mu^-
        TABS_MU_PLUS,      ///< partial entropy mu^+
        N_TAB_VARS_MUON
    } ;

    // ----------------------------------------------------------
    //  Cold-slice table indices (same layout as tabulated_eos_t)
    // ----------------------------------------------------------
    enum COLD_VIDX : int {
        CTABTEMP = 0,
        CTABYE,
        CTABYMU,
        CTABPRESS,
        CTABEPS,
        CTABCSND2,
        CTABENTROPY,
        N_CTAB_VARS
    } ;

    // ----------------------------------------------------------
    //  Default / constructors
    // ----------------------------------------------------------
    leptonic_eos_4d_t() = default ;

    /**
     * @brief Full constructor – called from read_leptonic_4d_table().
     */
    leptonic_eos_4d_t(
        // 4-D baryon table  [irho,iT,iye,iymu,ivar]
        Kokkos::View<double *****, grace::default_space> tab_baryon,
        Kokkos::View<double *,     grace::default_space> logrho,
        Kokkos::View<double *,     grace::default_space> logT,
        Kokkos::View<double *,     grace::default_space> ye_ax,
        Kokkos::View<double *,     grace::default_space> ymu_ax,
        // 4-D lepton tables [irho,iT,iye,iymu,ivar]
        Kokkos::View<double *****, grace::default_space> tab_ele,
        Kokkos::View<double *****, grace::default_space> tab_muon,
        // Cold slice  [irho, ivar]
        Kokkos::View<double **,    grace::default_space> cold_tab,
        Kokkos::View<double *,     grace::default_space> cold_logrho,
        // Thermodynamic range parameters
        double rhomax,   double rhomin,
        double tempmax,  double tempmin,
        double yemax,    double yemin,
        double ymumax,   double ymumin,
        double baryon_mass,
        double energy_shift,
        double c2p_epsmin, double c2p_epsmax,
        double c2p_hmin,   double c2p_hmax,
        double c2p_temp_atm,
        double c2p_ye_atm,
        double c2p_ymu_atm,
        bool   atmo_is_beta_eq
    )
    : tables(tab_baryon, logrho, logT, ye_ax, ymu_ax)
    , tables_ele (tab_ele,  logrho, logT, ye_ax, ymu_ax)
    , tables_muon(tab_muon, logrho, logT, ye_ax, ymu_ax)
    , cold_table(cold_tab, cold_logrho)
    , nrho(logrho.size()), nT(logT.size())
    , nye(ye_ax.size()),   nymu(ymu_ax.size())
    , energy_shift(energy_shift)
    , eos_ymumax(ymumax), eos_ymumin(ymumin)
    , _c2p_ymu_atm(c2p_ymu_atm)
    , base_t( rhomax, rhomin,
              tempmax, tempmin,
              yemax, yemin,
              baryon_mass,
              c2p_epsmin, c2p_epsmax,
              c2p_hmin,   c2p_hmax,
              c2p_temp_atm,
              c2p_ye_atm,
              atmo_is_beta_eq,
              false )
    {
        lrhomin  = logrho[0] ; lrhomax  = logrho[logrho.size()-1] ;
        ltempmin = logT[0]   ; ltempmax = logT[logT.size()-1]     ;
    }

    // ----------------------------------------------------------
    //  Public accessors – Y_mu limits
    // ----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION double get_c2p_ymu_min() const { return eos_ymumin      ; }
    KOKKOS_INLINE_FUNCTION double get_c2p_ymu_max() const { return eos_ymumax      ; }
    KOKKOS_INLINE_FUNCTION double get_c2p_ymu_atm() const { return _c2p_ymu_atm   ; }

    // ===========================================================
    //  CRTP implementation methods (called by eos_base_t)
    //  Signature convention identical to tabulated_eos_t.
    //  Where a 3D EOS takes (eps/rho/ye), here we add ymu.
    // ===========================================================

    // -----------------------------------------------------------
    //  press given eps, rho, ye, ymu
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    press__eps_rho_ye_impl(double& eps, double& rho, double& ye, err_t& err) const
    {
        // Fallback: call 4-D version with ymu = atm (needed for base API compat)
        double ymu = _c2p_ymu_atm ;
        return press__eps_rho_ye_ymu(eps, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION double
    press__eps_rho_ye_ymu(double& eps, double& rho, double& ye, double& ymu,
                          err_t& err) const
    {
        limit_rho(rho, err) ;
        limit_ye (ye,  err) ;
        limit_ymu(ymu, err) ;
        double lrho  = Kokkos::log(rho) ;
        double ltemp = ltemp__eps_lrho_ye_ymu(eps, lrho, ye, ymu, err) ;
        return press__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu) ;
    }

    // -----------------------------------------------------------
    //  press + temp given eps, rho, ye, ymu
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    press_temp__eps_rho_ye_impl(double& temp, double& eps, double& rho,
                                double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        return press_temp__eps_rho_ye_ymu(temp, eps, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION double
    press_temp__eps_rho_ye_ymu(double& temp, double& eps, double& rho,
                               double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err) ;
        limit_ye (ye,  err) ;
        limit_ymu(ymu, err) ;
        double lrho  = Kokkos::log(rho) ;
        double ltemp = ltemp__eps_lrho_ye_ymu(eps, lrho, ye, ymu, err) ;
        temp = Kokkos::exp(ltemp) ;
        return press__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu) ;
    }

    // -----------------------------------------------------------
    //  press given temp, rho, ye, ymu
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    press__temp_rho_ye_impl(double& temp, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        return press__temp_rho_ye_ymu(temp, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION double
    press__temp_rho_ye_ymu(double& temp, double& rho, double& ye, double& ymu,
                           err_t& err) const
    {
        limit_rho (rho,  err) ;
        limit_ye  (ye,   err) ;
        limit_ymu (ymu,  err) ;
        limit_temp(temp, err) ;
        double lrho  = Kokkos::log(rho)  ;
        double ltemp = Kokkos::log(temp) ;
        return press__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu) ;
    }

    // -----------------------------------------------------------
    //  eps given temp, rho, ye, ymu
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    eps__temp_rho_ye_impl(double& temp, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        return eps__temp_rho_ye_ymu(temp, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION double
    eps__temp_rho_ye_ymu(double& temp, double& rho, double& ye, double& ymu,
                         err_t& err) const
    {
        limit_rho (rho,  err) ;
        limit_ye  (ye,   err) ;
        limit_ymu (ymu,  err) ;
        limit_temp(temp, err) ;
        double lrho  = Kokkos::log(rho)  ;
        double ltemp = Kokkos::log(temp) ;
        return Kokkos::exp( tables.interp(lrho,ltemp,ye,ymu,TABEPS) ) - energy_shift ;
    }

    // -----------------------------------------------------------
    //  eps range  (min/max over temperature at fixed rho,ye,ymu)
    // -----------------------------------------------------------
    void KOKKOS_FUNCTION
    eps_range__rho_ye(double& epsmin, double& epsmax,
                      double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        eps_range__rho_ye_ymu(epsmin, epsmax, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION void
    eps_range__rho_ye_ymu(double& epsmin, double& epsmax,
                          double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err) ;
        limit_ye (ye,  err) ;
        limit_ymu(ymu, err) ;
        double lrho = Kokkos::log(rho) ;
        epsmin = Kokkos::exp( tables.interp(lrho,ltempmin,ye,ymu,TABEPS) ) - energy_shift ;
        epsmax = Kokkos::exp( tables.interp(lrho,ltempmax,ye,ymu,TABEPS) ) - energy_shift ;
    }

    // -----------------------------------------------------------
    //  entropy range
    // -----------------------------------------------------------
    void KOKKOS_FUNCTION
    entropy_range__rho_ye(double& smin, double& smax,
                          double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        entropy_range__rho_ye_ymu(smin, smax, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION void
    entropy_range__rho_ye_ymu(double& smin, double& smax,
                              double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err) ;
        limit_ye (ye,  err) ;
        limit_ymu(ymu, err) ;
        double lrho = Kokkos::log(rho) ;
        smin = tables.interp(lrho, ltempmin, ye, ymu, TABENTROPY) ;
        smax = tables.interp(lrho, ltempmax, ye, ymu, TABENTROPY) ;
    }

    // -----------------------------------------------------------
    //  Full hot-path: press + h + cs2 + temp + entropy  (eps-based)
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    press_h_csnd2_temp_entropy__eps_rho_ye_impl(
        double& h, double& csnd2, double& temp, double& entropy,
        double& eps, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        return press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
            h, csnd2, temp, entropy, eps, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION double
    press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
        double& h, double& csnd2, double& temp, double& entropy,
        double& eps, double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err) ;
        limit_ye (ye,  err) ;
        limit_ymu(ymu, err) ;
        double lrho  = Kokkos::log(rho) ;
        double ltemp = ltemp__eps_lrho_ye_ymu(eps, lrho, ye, ymu, err) ;
        temp = Kokkos::exp(ltemp) ;

        constexpr std::array<int,4> vidx{ TABPRESS, TABCSND2, TABENTROPY } ;
        std::array<double,3> res ;
        std::array<int,3> _v{ TABPRESS, TABCSND2, TABENTROPY } ;
        tables.interp<3>(lrho, ltemp, ye, ymu, _v, res) ;

        double press = Kokkos::exp(res[0]) ;
        csnd2   = res[1] ;
        entropy = res[2] ;
        h = 1. + eps + press / rho ;
        return press ;
    }

    // -----------------------------------------------------------
    //  Full hot-path: press + h + cs2 + temp + eps  (entropy-based)
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    press_h_csnd2_temp_eps__entropy_rho_ye_impl(
        double& h, double& csnd2, double& temp, double& eps,
        double& entropy, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        return press_h_csnd2_temp_eps__entropy_rho_ye_ymu(
            h, csnd2, temp, eps, entropy, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION double
    press_h_csnd2_temp_eps__entropy_rho_ye_ymu(
        double& h, double& csnd2, double& temp, double& eps,
        double& entropy, double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err) ;
        limit_ye (ye,  err) ;
        limit_ymu(ymu, err) ;
        double lrho  = Kokkos::log(rho) ;
        double ltemp = ltemp__entropy_lrho_ye_ymu(entropy, lrho, ye, ymu) ;
        temp = Kokkos::exp(ltemp) ;

        std::array<int,3> _v{ TABPRESS, TABCSND2, TABEPS } ;
        std::array<double,3> res ;
        tables.interp<3>(lrho, ltemp, ye, ymu, _v, res) ;

        double press = Kokkos::exp(res[0]) ;
        csnd2 = res[1] ;
        eps   = Kokkos::exp(res[2]) - energy_shift ;
        h     = 1. + eps + press / rho ;
        return press ;
    }

    // -----------------------------------------------------------
    //  Convenience:  press + eps + cs2 given temp, rho, ye, ymu
    // -----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION double
    press_eps_csnd2__temp_rho_ye_ymu(double& eps, double& csnd2,
                                      double& temp, double& rho,
                                      double& ye, double& ymu, err_t& err) const
    {
        limit_rho (rho,  err) ;
        limit_ye  (ye,   err) ;
        limit_ymu (ymu,  err) ;
        limit_temp(temp, err) ;
        double lrho  = Kokkos::log(rho)  ;
        double ltemp = Kokkos::log(temp) ;
        std::array<int,3> _v{ TABPRESS, TABEPS, TABCSND2 } ;
        std::array<double,3> res ;
        tables.interp<3>(lrho, ltemp, ye, ymu, _v, res) ;
        double press = Kokkos::exp(res[0]) ;
        eps   = Kokkos::exp(res[1]) - energy_shift ;
        csnd2 = res[2] ;
        return press ;
    }

    // -----------------------------------------------------------
    //  Convenience: press + eps + cs2 + entropy given temp
    //  (used by flux routines) – uses base API default ymu
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    press_eps_csnd2_entropy__temp_rho_ye_impl(
        double& eps, double& csnd2, double& entropy,
        double& temp, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        return press_eps_csnd2_entropy__temp_rho_ye_ymu(
            eps, csnd2, entropy, temp, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION double
    press_eps_csnd2_entropy__temp_rho_ye_ymu(
        double& eps, double& csnd2, double& entropy,
        double& temp, double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho (rho,  err) ;
        limit_ye  (ye,   err) ;
        limit_ymu (ymu,  err) ;
        limit_temp(temp, err) ;
        double lrho  = Kokkos::log(rho)  ;
        double ltemp = Kokkos::log(temp) ;
        std::array<int,4> _v{ TABPRESS, TABEPS, TABCSND2, TABENTROPY } ;
        std::array<double,4> res ;
        tables.interp<4>(lrho, ltemp, ye, ymu, _v, res) ;
        double press = Kokkos::exp(res[0]) ;
        eps     = Kokkos::exp(res[1]) - energy_shift ;
        csnd2   = res[2] ;
        entropy = res[3] ;
        return press ;
    }

    // -----------------------------------------------------------
    //  press + h + cs2  (economy interface – eps-based)
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    press_h_csnd2__eps_rho_ye_impl(double& h, double& csnd2,
                                    double& eps, double& rho,
                                    double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        double temp_dummy ;
        double entropy_dummy ;
        return press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
            h, csnd2, temp_dummy, entropy_dummy, eps, rho, ye, ymu, err) ;
    }

    // -----------------------------------------------------------
    //  press + h + cs2  (economy interface – temp-based)
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    press_h_csnd2__temp_rho_ye_impl(double& h, double& csnd2,
                                     double& temp, double& rho,
                                     double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm ;
        return press_h_csnd2__temp_rho_ye_ymu(h, csnd2, temp, rho, ye, ymu, err) ;
    }

    KOKKOS_INLINE_FUNCTION double
    press_h_csnd2__temp_rho_ye_ymu(double& h, double& csnd2,
                                    double& temp, double& rho,
                                    double& ye, double& ymu, err_t& err) const
    {
        limit_rho (rho,  err) ;
        limit_ye  (ye,   err) ;
        limit_ymu (ymu,  err) ;
        limit_temp(temp, err) ;
        double lrho  = Kokkos::log(rho)  ;
        double ltemp = Kokkos::log(temp) ;
        std::array<int,3> _v{ TABPRESS, TABCSND2, TABEPS } ;
        std::array<double,3> res ;
        tables.interp<3>(lrho, ltemp, ye, ymu, _v, res) ;
        double press = Kokkos::exp(res[0]) ;
        double eps   = Kokkos::exp(res[2]) - energy_shift ;
        csnd2 = res[1] ;
        h     = 1. + eps + press / rho ;
        return press ;
    }

    // -----------------------------------------------------------
    //  Cold-slice accessors
    // -----------------------------------------------------------
    double GRACE_HOST_DEVICE
    temp_cold__rho_impl(double& rho, err_t& err) const
    {
        limit_rho(rho, err) ;
        double lrho = Kokkos::log(rho) ;
        return Kokkos::exp(cold_table.interp(lrho, CTABTEMP)) ;
    }

    double GRACE_HOST_DEVICE
    entropy_cold__rho_impl(double& rho, err_t& err) const
    {
        limit_rho(rho, err) ;
        double lrho = Kokkos::log(rho) ;
        return cold_table.interp(lrho, CTABENTROPY) ;
    }

    // -----------------------------------------------------------
    //  Muon chemical potential  mu_mu(rho, T, ye, ymu)
    // -----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION double
    mumu__temp_rho_ye_ymu(double& temp, double& rho, double& ye, double& ymu,
                          err_t& err) const
    {
        limit_rho (rho,  err) ;
        limit_ye  (ye,   err) ;
        limit_ymu (ymu,  err) ;
        limit_temp(temp, err) ;
        double lrho  = Kokkos::log(rho)  ;
        double ltemp = Kokkos::log(temp) ;
        return tables_muon.interp(lrho, ltemp, ye, ymu, TABMUMU) ;
    }

    // -----------------------------------------------------------
    //  Electron chemical potential  mu_e(rho, T, ye, ymu)
    // -----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION double
    mue__temp_rho_ye_ymu(double& temp, double& rho, double& ye, double& ymu,
                         err_t& err) const
    {
        limit_rho (rho,  err) ;
        limit_ye  (ye,   err) ;
        limit_ymu (ymu,  err) ;
        limit_temp(temp, err) ;
        double lrho  = Kokkos::log(rho)  ;
        double ltemp = Kokkos::log(temp) ;
        return tables.interp(lrho, ltemp, ye, ymu, TABMUE) ;
    }

    // -----------------------------------------------------------
    //  Full chemical potentials + composition
    // -----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION double
    mue_mumu_mup_mun__temp_rho_ye_ymu(
        double& mumu, double& mup, double& mun,
        double& temp, double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho (rho,  err) ;
        limit_ye  (ye,   err) ;
        limit_ymu (ymu,  err) ;
        limit_temp(temp, err) ;
        double lrho  = Kokkos::log(rho)  ;
        double ltemp = Kokkos::log(temp) ;
        std::array<int,4> _v{ TABMUE, TABMUP, TABMUN, -1 } ;
        // mumu from muon table
        mumu = tables_muon.interp(lrho, ltemp, ye, ymu, TABMUMU) ;
        // baryon table: mue, mup, mun
        std::array<int,3> _vb{ TABMUE, TABMUP, TABMUN } ;
        std::array<double,3> res ;
        tables.interp<3>(lrho, ltemp, ye, ymu, _vb, res) ;
        mup = res[1] ;
        mun = res[2] ;
        return res[0] ; // returns mue
    }

    // -----------------------------------------------------------
    //  Beta-equilibrium solvers
    //  (needed for atmosphere initialisation and ID converter)
    // -----------------------------------------------------------

    /**
     * @brief Find Y_e, Y_mu in full (neutrino-less) beta equilibrium
     *        at given rho, T.  Solves:
     *          mu_n - mu_p - mu_e = 0   (charge neutrality + e-capture)
     *          mu_mu - mu_e = 0          (muon threshold)
     *          Y_e = Y_e(Y_mu, rho, T)  from table
     * Host-only due to Brent iteration.
     */
    void
    press_eps_ye_ymu__beta_eq__rho_temp(
        double& press, double& eps,
        double& ye, double& ymu,
        double& rho, double& temp, err_t& err) const ;

    // ===========================================================
    //  Public data members (needed by GPU kernels via capture)
    // ===========================================================
    tabeos_linterp_4d_t tables      ; ///< 4D baryon table interpolator
    tabeos_linterp_4d_t tables_ele  ; ///< 4D electron lepton table interpolator
    tabeos_linterp_4d_t tables_muon ; ///< 4D muon lepton table interpolator

    // Cold-slice 1-D table  (same struct as tabulated_eos_t)
    struct cold_eos_linterp_1d_t {
        cold_eos_linterp_1d_t() = default ;
        cold_eos_linterp_1d_t(
            Kokkos::View<double **> tabs,
            Kokkos::View<double *>  ar
        ) : _tables(tabs), _logrho(ar)
        { idr = 1./(_logrho(1)-_logrho(0)) ; }

        KOKKOS_INLINE_FUNCTION double interp(double lrho, int idx) const {
            int i = Kokkos::max(0UL,
                    Kokkos::min( static_cast<size_t>((lrho-_logrho(0))*idr),
                                 _logrho.extent(0)-2 )) ;
            double lam = (lrho - _logrho(i)) * idr ;
            return (1.-lam)*_tables(i,idx) + lam*_tables(i+1,idx) ;
        }
        Kokkos::View<const double*> _logrho ;
        Kokkos::View<double **>     _tables  ;
        double idr ;
    } cold_table ;

    int nrho, nT, nye, nymu ;
    double energy_shift ;
    double lrhomin, lrhomax ;
    double ltempmin, ltempmax ;
    double eos_ymumin, eos_ymumax ;

  private:

    double _c2p_ymu_atm ;

    // ----------------------------------------------------------
    //  Clamping helpers
    // ----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION void limit_rho(double& rho, err_t& err) const {
        if ( rho < this->eos_rhomin ) { rho = this->eos_rhomin ; err.set(EOS_RHO_TOO_LOW)  ; }
        if ( rho > this->eos_rhomax ) { rho = this->eos_rhomax ; err.set(EOS_RHO_TOO_HIGH) ; }
    }
    KOKKOS_INLINE_FUNCTION void limit_ye(double& ye, err_t& err) const {
        if ( ye < this->eos_yemin ) { ye = this->eos_yemin ; err.set(EOS_YE_TOO_LOW)  ; }
        if ( ye > this->eos_yemax ) { ye = this->eos_yemax ; err.set(EOS_YE_TOO_HIGH) ; }
    }
    KOKKOS_INLINE_FUNCTION void limit_ymu(double& ymu, err_t& err) const {
        if ( ymu < eos_ymumin ) { ymu = eos_ymumin ; err.set(EOS_YE_TOO_LOW)  ; } // reuse ye slot
        if ( ymu > eos_ymumax ) { ymu = eos_ymumax ; err.set(EOS_YE_TOO_HIGH) ; }
    }
    KOKKOS_INLINE_FUNCTION void limit_temp(double& temp, err_t& err) const {
        double tmin = Kokkos::exp(ltempmin) ;
        double tmax = Kokkos::exp(ltempmax) ;
        if ( temp < tmin ) { temp = tmin ; err.set(EOS_TEMPERATURE_TOO_LOW)  ; }
        if ( temp > tmax ) { temp = tmax ; err.set(EOS_TEMPERATURE_TOO_HIGH) ; }
    }

    // ----------------------------------------------------------
    //  Root-find for log(T) given eps – 4D path
    // ----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION double
    ltemp__eps_lrho_ye_ymu(double& eps, double lrho, double ye, double ymu,
                            err_t& err) const
    {
        double leps = Kokkos::log(eps + energy_shift) ;
        double lepsmin = tables.interp(lrho, ltempmin, ye, ymu, TABEPS) ;
        double lepsmax = tables.interp(lrho, ltempmax, ye, ymu, TABEPS) ;
        if ( leps <= lepsmin ) {
            eps = Kokkos::exp(lepsmin) - energy_shift ;
            err.set(EOS_EPS_TOO_LOW) ;
            return ltempmin ;
        }
        if ( leps >= lepsmax ) {
            eps = Kokkos::exp(lepsmax) - energy_shift ;
            err.set(EOS_EPS_TOO_HIGH) ;
            return ltempmax ;
        }
        auto rootfun = [this,lrho,ye,ymu,leps](double lt){
            return tables.interp(lrho,lt,ye,ymu,TABEPS) - leps ;
        } ;
        return utils::brent(rootfun, ltempmin, ltempmax, 1e-14) ;
    }

    // ----------------------------------------------------------
    //  Root-find for log(T) given entropy – 4D path
    // ----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION double
    ltemp__entropy_lrho_ye_ymu(double entropy, double lrho,
                                double ye, double ymu) const
    {
        auto rootfun = [this,lrho,ye,ymu,entropy](double lt){
            return tables.interp(lrho,lt,ye,ymu,TABENTROPY) - entropy ;
        } ;
        return utils::brent(rootfun, ltempmin, ltempmax, 1e-14) ;
    }

    // ----------------------------------------------------------
    //  Direct pressure lookup  (no unit conversion needed –
    //  tables already store log(press) in geometric units)
    // ----------------------------------------------------------
    KOKKOS_INLINE_FUNCTION double
    press__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                              double ye,   double ymu) const
    {
        return Kokkos::exp(tables.interp(lrho, ltemp, ye, ymu, TABPRESS)) ;
    }

} ; // class leptonic_eos_4d_t

} /* namespace grace */

#endif /* GRACE_PHYSICS_EOS_LEPTONIC_4D_HH */
