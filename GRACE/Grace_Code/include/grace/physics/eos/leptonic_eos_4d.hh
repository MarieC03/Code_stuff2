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
#include <grace/physics/eos/eos_base.hh>
#include <grace/utils/bitset.hh>
#include <grace/utils/grace_utils.hh>
#include <grace/utils/rootfinding.hh>

#include <Kokkos_Core.hpp>

#include <array>
#include <cmath>
#include <limits>

namespace grace {

class leptonic_eos_4d_t;

struct tabeos_linterp_4d_t {

    tabeos_linterp_4d_t() = default;

    tabeos_linterp_4d_t(
        Kokkos::View<double *****> tabs,
        Kokkos::View<double *> ar,
        Kokkos::View<double *> at,
        Kokkos::View<double *> ay,
        Kokkos::View<double *> aym)
        : _tables(tabs), _logrho(ar), _logT(at), _ye(ay), _ymu(aym)
    {
        idr = 1. / (_logrho(1) - _logrho(0));
        idt = 1. / (_logT(1) - _logT(0));
        idy = 1. / (_ye(1) - _ye(0));
        idym = 1. / (_ymu(1) - _ymu(0));
    }

    KOKKOS_INLINE_FUNCTION double lrho(int i) const { return _logrho(i); }
    KOKKOS_INLINE_FUNCTION double ltemp(int i) const { return _logT(i); }
    KOKKOS_INLINE_FUNCTION double ye(int i) const { return _ye(i); }
    KOKKOS_INLINE_FUNCTION double ymu(int i) const { return _ymu(i); }

    KOKKOS_INLINE_FUNCTION
    double interp(double lrho, double ltemp, double ye, double ymu, int idx) const
    {
        std::array<double, 1> res{};
        std::array<int, 1> ids{idx};
        interp<1>(lrho, ltemp, ye, ymu, ids, res);
        return res[0];
    }

    template<int N>
    KOKKOS_INLINE_FUNCTION void
    interp(double lrho, double ltemp, double ye, double ymu,
           std::array<int, N> const& idx,
           std::array<double, N>& res) const
    {
        for (int iv = 0; iv < N; ++iv) {
            res[iv] = 0.0;
        }
        int ir{}, it{}, iy{}, iym{};
        getidx(lrho, ltemp, ye, ymu, ir, it, iy, iym);

        double wr[2], wt[2], wy[2], wym[2];
        getw(lrho, ir, _logrho, idr, wr);
        getw(ltemp, it, _logT, idt, wt);
        getw(ye, iy, _ye, idy, wy);
        getw(ymu, iym, _ymu, idym, wym);

        for (int ii = 0; ii < 2; ++ii) {
            for (int jj = 0; jj < 2; ++jj) {
                for (int kk = 0; kk < 2; ++kk) {
                    for (int ll = 0; ll < 2; ++ll) {
                        double const w = wr[ii] * wt[jj] * wy[kk] * wym[ll];
                        for (int iv = 0; iv < N; ++iv) {
                            res[iv] += w * _tables(ir + ii, it + jj, iy + kk, iym + ll, idx[iv]);
                        }
                    }
                }
            }
        }
    }

    KOKKOS_INLINE_FUNCTION void
    getidx(double lr, double lt, double ye_in, double ym_in,
           int& ir, int& it, int& iy, int& iym) const
    {
        ir = Kokkos::max(
            0UL,
            Kokkos::min(static_cast<size_t>((lr - _logrho(0)) * idr), _logrho.extent(0) - 2));
        it = Kokkos::max(
            0UL,
            Kokkos::min(static_cast<size_t>((lt - _logT(0)) * idt), _logT.extent(0) - 2));
        iy = Kokkos::max(
            0UL,
            Kokkos::min(static_cast<size_t>((ye_in - _ye(0)) * idy), _ye.extent(0) - 2));
        iym = Kokkos::max(
            0UL,
            Kokkos::min(static_cast<size_t>((ym_in - _ymu(0)) * idym), _ymu.extent(0) - 2));
    }

    KOKKOS_INLINE_FUNCTION void
    getw(double x, int i, Kokkos::View<const double*> ax, double ih, double* w) const
    {
        double const lam = (x - ax(i)) * ih;
        w[0] = 1.0 - lam;
        w[1] = lam;
    }

    Kokkos::View<const double*> _logrho, _logT, _ye, _ymu;
    Kokkos::View<double *****> _tables;
    double idr{}, idt{}, idy{}, idym{};
};

class leptonic_eos_4d_t : public eos_base_t<leptonic_eos_4d_t>
{
    using err_t = eos_err_t;
    using base_t = eos_base_t<leptonic_eos_4d_t>;

  public:
    using error_type = eos_err_t;
    static constexpr bool has_ye = true;

    enum LEP4D_ERROR_T : int {
        EOS_YMU_TOO_LOW = EOS_NUM_ERRORS,
        EOS_YMU_TOO_HIGH,
        EOS_4D_NUM_ERRORS
    };

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
    };

    enum ELE_VIDX : int {
        TABMUELE = 0,
        TABYLE_MINUS,
        TABYLE_PLUS,
        TABPRESS_E_MINUS,
        TABPRESS_E_PLUS,
        TABEPS_E_MINUS,
        TABEPS_E_PLUS,
        TABS_E_MINUS,
        TABS_E_PLUS,
        N_TAB_VARS_ELE
    };

    enum MUON_VIDX : int {
        TABMUMU = 0,
        TABYMU_MINUS,
        TABYMU_PLUS,
        TABPRESS_MU_MINUS,
        TABPRESS_MU_PLUS,
        TABEPS_MU_MINUS,
        TABEPS_MU_PLUS,
        TABS_MU_MINUS,
        TABS_MU_PLUS,
        N_TAB_VARS_MUON
    };

    enum COLD_VIDX : int {
        CTABTEMP = 0,
        CTABYE,
        CTABYMU,
        CTABPRESS,
        CTABEPS,
        CTABCSND2,
        CTABENTROPY,
        N_CTAB_VARS
    };

    leptonic_eos_4d_t() = default;

    leptonic_eos_4d_t(
        Kokkos::View<double *****, grace::default_space> tab_baryon,
        Kokkos::View<double *, grace::default_space> logrho,
        Kokkos::View<double *, grace::default_space> logT,
        Kokkos::View<double *, grace::default_space> ye_ax,
        Kokkos::View<double *, grace::default_space> ymu_ax,
        Kokkos::View<double *****, grace::default_space> tab_ele,
        Kokkos::View<double *****, grace::default_space> tab_muon,
        Kokkos::View<double **, grace::default_space> cold_tab,
        Kokkos::View<double *, grace::default_space> cold_logrho,
        double rhomax,
        double rhomin,
        double tempmax,
        double tempmin,
        double yemax,
        double yemin,
        double ymumax,
        double ymumin,
        double baryon_mass,
        double energy_shift_in,
        double c2p_epsmin,
        double c2p_epsmax,
        double c2p_hmin,
        double c2p_hmax,
        double c2p_temp_atm,
        double c2p_ye_atm,
        double c2p_ymu_atm,
        bool atmo_is_beta_eq)
        : tables(tab_baryon, logrho, logT, ye_ax, ymu_ax)
        , tables_ele(tab_ele, logrho, logT, ye_ax, ymu_ax)
        , tables_muon(tab_muon, logrho, logT, ye_ax, ymu_ax)
        , cold_table(cold_tab, cold_logrho)
        , nrho(logrho.size())
        , nT(logT.size())
        , nye(ye_ax.size())
        , nymu(ymu_ax.size())
        , energy_shift(energy_shift_in)
        , lrhomin(logrho[0])
        , lrhomax(logrho[logrho.size() - 1])
        , ltempmin(logT[0])
        , ltempmax(logT[logT.size() - 1])
        , eos_ymumin(ymumin)
        , eos_ymumax(ymumax)
        , _c2p_ymu_atm(c2p_ymu_atm)
        , base_t(rhomax, rhomin,
                 tempmax, tempmin,
                 yemax, yemin,
                 baryon_mass,
                 c2p_epsmin, c2p_epsmax,
                 c2p_hmin, c2p_hmax,
                 c2p_temp_atm,
                 c2p_ye_atm,
                 atmo_is_beta_eq,
                 false)
    {}

    KOKKOS_INLINE_FUNCTION double get_c2p_ymu_min() const { return eos_ymumin; }
    KOKKOS_INLINE_FUNCTION double get_c2p_ymu_max() const { return eos_ymumax; }
    KOKKOS_INLINE_FUNCTION double get_c2p_ymu_atm() const { return _c2p_ymu_atm; }

    void set_atmosphere_fractions(double ye_atm, double ymu_atm)
    {
        this->c2p_ye_atm = ye_atm;
        _c2p_ymu_atm = ymu_atm;
    }

    double GRACE_HOST_DEVICE
    press__eps_rho_ye_impl(double& eps, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return press__eps_rho_ye_ymu(eps, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    press__eps_rho_ye_ymu(double& eps, double& rho, double& ye, double& ymu,
                          err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = ltemp__eps_lrho_ye_ymu(eps, lrho, ye, ymu, err);
        return press__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu);
    }

    double GRACE_HOST_DEVICE
    press_temp__eps_rho_ye_impl(double& temp, double& eps, double& rho,
                                double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return press_temp__eps_rho_ye_ymu(temp, eps, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    press_temp__eps_rho_ye_ymu(double& temp, double& eps, double& rho,
                               double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = ltemp__eps_lrho_ye_ymu(eps, lrho, ye, ymu, err);
        temp = Kokkos::exp(ltemp);
        return press__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu);
    }

    double GRACE_HOST_DEVICE
    press__temp_rho_ye_impl(double& temp, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return press__temp_rho_ye_ymu(temp, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    press__temp_rho_ye_ymu(double& temp, double& rho, double& ye, double& ymu,
                           err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);
        return press__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu);
    }

    double GRACE_HOST_DEVICE
    eps__temp_rho_ye_impl(double& temp, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return eps__temp_rho_ye_ymu(temp, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    eps__temp_rho_ye_ymu(double& temp, double& rho, double& ye, double& ymu,
                         err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);
        return eps__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu);
    }

    void GRACE_HOST_DEVICE
    eps_range__rho_ye_impl(double& epsmin, double& epsmax,
                           double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        eps_range__rho_ye_ymu(epsmin, epsmax, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION void
    eps_range__rho_ye_ymu(double& epsmin, double& epsmax,
                          double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        double const lrho = Kokkos::log(rho);
        epsmin = eps__lrho_ltemp_ye_ymu(lrho, ltempmin, ye, ymu);
        epsmax = eps__lrho_ltemp_ye_ymu(lrho, ltempmax, ye, ymu);
    }

    void GRACE_HOST_DEVICE
    entropy_range__rho_ye_impl(double& smin, double& smax,
                               double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        entropy_range__rho_ye_ymu(smin, smax, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION void
    entropy_range__rho_ye_ymu(double& smin, double& smax,
                              double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        double const lrho = Kokkos::log(rho);
        smin = entropy__lrho_ltemp_ye_ymu(lrho, ltempmin, ye, ymu);
        smax = entropy__lrho_ltemp_ye_ymu(lrho, ltempmax, ye, ymu);
    }

    double GRACE_HOST_DEVICE
    press_h_csnd2_temp_entropy__eps_rho_ye_impl(
        double& h, double& csnd2, double& temp, double& entropy,
        double& eps, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
            h, csnd2, temp, entropy, eps, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
        double& h, double& csnd2, double& temp, double& entropy,
        double& eps, double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = ltemp__eps_lrho_ye_ymu(eps, lrho, ye, ymu, err);
        temp = Kokkos::exp(ltemp);
        double press{};
        thermo_press_eps_csnd2_entropy__lrho_ltemp_ye_ymu(
            press, eps, csnd2, entropy, lrho, ltemp, ye, ymu);
        h = 1.0 + eps + press / rho;
        return press;
    }

    double GRACE_HOST_DEVICE
    press_h_csnd2_temp_eps__entropy_rho_ye_impl(
        double& h, double& csnd2, double& temp, double& eps,
        double& entropy, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return press_h_csnd2_temp_eps__entropy_rho_ye_ymu(
            h, csnd2, temp, eps, entropy, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    press_h_csnd2_temp_eps__entropy_rho_ye_ymu(
        double& h, double& csnd2, double& temp, double& eps,
        double& entropy, double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_entropy_rho_ye_ymu(entropy, rho, ye, ymu, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = ltemp__entropy_lrho_ye_ymu(entropy, lrho, ye, ymu, err);
        temp = Kokkos::exp(ltemp);
        double press{};
        thermo_press_eps_csnd2_entropy__lrho_ltemp_ye_ymu(
            press, eps, csnd2, entropy, lrho, ltemp, ye, ymu);
        h = 1.0 + eps + press / rho;
        return press;
    }

    double GRACE_HOST_DEVICE
    press_eps_csnd2__temp_rho_ye_impl(double& eps, double& csnd2, double& temp,
                                      double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return press_eps_csnd2__temp_rho_ye_ymu(eps, csnd2, temp, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    press_eps_csnd2__temp_rho_ye_ymu(double& eps, double& csnd2,
                                     double& temp, double& rho,
                                     double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);
        double dummy_entropy{};
        double press{};
        thermo_press_eps_csnd2_entropy__lrho_ltemp_ye_ymu(
            press, eps, csnd2, dummy_entropy, lrho, ltemp, ye, ymu);
        return press;
    }

    double GRACE_HOST_DEVICE
    eps_csnd2_entropy__temp_rho_ye_impl(double& csnd2, double& entropy,
                                        double& temp, double& rho,
                                        double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return eps_csnd2_entropy__temp_rho_ye_ymu(csnd2, entropy, temp, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    eps_csnd2_entropy__temp_rho_ye_ymu(double& csnd2, double& entropy,
                                       double& temp, double& rho,
                                       double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);
        double press{};
        double eps{};
        thermo_press_eps_csnd2_entropy__lrho_ltemp_ye_ymu(
            press, eps, csnd2, entropy, lrho, ltemp, ye, ymu);
        return eps;
    }

    double GRACE_HOST_DEVICE
    press_eps_csnd2_entropy__temp_rho_ye_impl(
        double& eps, double& csnd2, double& entropy,
        double& temp, double& rho, double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return press_eps_csnd2_entropy__temp_rho_ye_ymu(
            eps, csnd2, entropy, temp, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    press_eps_csnd2_entropy__temp_rho_ye_ymu(
        double& eps, double& csnd2, double& entropy,
        double& temp, double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);
        double press{};
        thermo_press_eps_csnd2_entropy__lrho_ltemp_ye_ymu(
            press, eps, csnd2, entropy, lrho, ltemp, ye, ymu);
        return press;
    }

    double GRACE_HOST_DEVICE
    press_h_csnd2__eps_rho_ye_impl(double& h, double& csnd2,
                                   double& eps, double& rho,
                                   double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        double temp{};
        double entropy{};
        return press_h_csnd2_temp_entropy__eps_rho_ye_ymu(
            h, csnd2, temp, entropy, eps, rho, ye, ymu, err);
    }

    double GRACE_HOST_DEVICE
    press_h_csnd2__temp_rho_ye_impl(double& h, double& csnd2,
                                    double& temp, double& rho,
                                    double& ye, err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        return press_h_csnd2__temp_rho_ye_ymu(h, csnd2, temp, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    press_h_csnd2__temp_rho_ye_ymu(double& h, double& csnd2,
                                   double& temp, double& rho,
                                   double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);
        double press{};
        double eps{};
        double entropy{};
        thermo_press_eps_csnd2_entropy__lrho_ltemp_ye_ymu(
            press, eps, csnd2, entropy, lrho, ltemp, ye, ymu);
        h = 1.0 + eps + press / rho;
        return press;
    }

    double GRACE_ALWAYS_INLINE GRACE_HOST_DEVICE
    eps_h_csnd2_temp_entropy__press_rho_ye_impl(
        double& h, double& csnd2, double& temp, double& entropy,
        double& press, double& rho, double& ye, err_t& err) const
    {
        Kokkos::abort(
            "The function eps_h_csnd2_temp_entropy__press_rho_ye_impl is not well defined for a leptonic tabulated EOS.");
    }

    double GRACE_HOST_DEVICE
    press_cold__rho_impl(double& rho, err_t& err) const
    {
        limit_rho(rho, err);
        double const lrho = Kokkos::log(rho);
        return Kokkos::exp(cold_table.interp(lrho, CTABPRESS));
    }

    double GRACE_HOST_DEVICE
    rho__press_cold_impl(double& press_cold, err_t& err) const
    {
        double const lp = cold_lpress__press_limited(press_cold, err);
        auto rootfun = [this, lp](double lrho) {
            return cold_table.interp(lrho, CTABPRESS) - lp;
        };
        double const lrmin = cold_table._logrho(0);
        double const lrmax = cold_table._logrho(cold_table._logrho.size() - 1);
        return Kokkos::exp(utils::brent(rootfun, lrmin, lrmax, 1e-14));
    }

    double GRACE_HOST_DEVICE
    rho__energy_cold_impl(double& e_cold, err_t& err) const
    {
        int const n = cold_table._logrho.size();
        double const eps_min = Kokkos::exp(cold_table._tables(0, CTABEPS)) - energy_shift;
        double const eps_max = Kokkos::exp(cold_table._tables(n - 1, CTABEPS)) - energy_shift;
        double const e_min = (1.0 + eps_min) * Kokkos::exp(cold_table._logrho(0));
        double const e_max = (1.0 + eps_max) * Kokkos::exp(cold_table._logrho(n - 1));
        if (e_cold < e_min) {
            e_cold = e_min;
            err.set(EOS_EPS_TOO_LOW);
            return Kokkos::exp(cold_table._logrho(0));
        }
        if (e_cold > e_max) {
            e_cold = e_max;
            err.set(EOS_EPS_TOO_HIGH);
            return Kokkos::exp(cold_table._logrho(n - 1));
        }
        auto rootfun = [this, e_cold](double lrho) {
            double const eps = Kokkos::exp(cold_table.interp(lrho, CTABEPS)) - energy_shift;
            return (1.0 + eps) * Kokkos::exp(lrho) - e_cold;
        };
        return Kokkos::exp(utils::brent(
            rootfun, cold_table._logrho(0), cold_table._logrho(n - 1), 1e-14));
    }

    double GRACE_HOST_DEVICE
    energy_cold__press_cold_impl(double& press_cold, err_t& err) const
    {
        double const lp = cold_lpress__press_limited(press_cold, err);
        auto rootfun = [this, lp](double lrho) {
            return cold_table.interp(lrho, CTABPRESS) - lp;
        };
        double const lrmin = cold_table._logrho(0);
        double const lrmax = cold_table._logrho(cold_table._logrho.size() - 1);
        double const lrho = utils::brent(rootfun, lrmin, lrmax, 1e-14);
        double const eps = Kokkos::exp(cold_table.interp(lrho, CTABEPS)) - energy_shift;
        return Kokkos::exp(lrho) * (1.0 + eps);
    }

    double GRACE_HOST_DEVICE
    eps_cold__rho_impl(double& rho, err_t& err) const
    {
        limit_rho(rho, err);
        double const lrho = Kokkos::log(rho);
        return Kokkos::exp(cold_table.interp(lrho, CTABEPS)) - energy_shift;
    }

    double GRACE_HOST_DEVICE
    ye_cold__rho_impl(double& rho, err_t& err) const
    {
        limit_rho(rho, err);
        double const lrho = Kokkos::log(rho);
        return cold_table.interp(lrho, CTABYE);
    }

    double GRACE_HOST_DEVICE
    ye_cold__press_impl(double& press, err_t& err) const
    {
        double const rho = rho__press_cold_impl(press, err);
        double const lrho = Kokkos::log(rho);
        return cold_table.interp(lrho, CTABYE);
    }

    double GRACE_HOST_DEVICE
    temp_cold__rho_impl(double& rho, err_t& err) const
    {
        limit_rho(rho, err);
        double const lrho = Kokkos::log(rho);
        return Kokkos::exp(cold_table.interp(lrho, CTABTEMP));
    }

    double GRACE_HOST_DEVICE
    entropy_cold__rho_impl(double& rho, err_t& err) const
    {
        limit_rho(rho, err);
        double const lrho = Kokkos::log(rho);
        return cold_table.interp(lrho, CTABENTROPY);
    }

    double GRACE_HOST_DEVICE
    mue_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_ye_impl(
        double& mup, double& mun, double& Xa, double& Xh, double& Xn, double& Xp,
        double& Abar, double& Zbar, double& temp, double& rho, double& ye,
        err_t& err) const
    {
        double ymu = _c2p_ymu_atm;
        double mumu{};
        return mue_mumu_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_ye_ymu(
            mumu, mup, mun, Xa, Xh, Xn, Xp, Abar, Zbar, temp, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    mumu__temp_rho_ye_ymu(double& temp, double& rho, double& ye, double& ymu,
                          err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);
        return tables_muon.interp(lrho, ltemp, ye, ymu, TABMUMU);
    }

    KOKKOS_INLINE_FUNCTION double
    mue__temp_rho_ye_ymu(double& temp, double& rho, double& ye, double& ymu,
                         err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);
        return tables_ele.interp(lrho, ltemp, ye, ymu, TABMUELE);
    }

    KOKKOS_INLINE_FUNCTION double
    mue_mumu_mup_mun__temp_rho_ye_ymu(
        double& mumu, double& mup, double& mun,
        double& temp, double& rho, double& ye, double& ymu, err_t& err) const
    {
        double Xa{}, Xh{}, Xn{}, Xp{}, Abar{}, Zbar{};
        return mue_mumu_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_ye_ymu(
            mumu, mup, mun, Xa, Xh, Xn, Xp, Abar, Zbar,
            temp, rho, ye, ymu, err);
    }

    KOKKOS_INLINE_FUNCTION double
    mue_mumu_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_ye_ymu(
        double& mumu, double& mup, double& mun,
        double& Xa, double& Xh, double& Xn, double& Xp,
        double& Abar, double& Zbar,
        double& temp, double& rho, double& ye, double& ymu, err_t& err) const
    {
        limit_rho(rho, err);
        limit_ye(ye, err);
        limit_ymu(ymu, err);
        limit_temp(temp, err);
        double const lrho = Kokkos::log(rho);
        double const ltemp = Kokkos::log(temp);

        std::array<int, 8> baryon_ids{
            TABMUP, TABMUN, TABXA, TABXH, TABXN, TABXP, TABABAR, TABZBAR};
        std::array<double, 8> baryon_vals{};
        tables.interp<8>(lrho, ltemp, ye, ymu, baryon_ids, baryon_vals);

        mup = baryon_vals[0];
        mun = baryon_vals[1];
        Xa = baryon_vals[2];
        Xh = baryon_vals[3];
        Xn = baryon_vals[4];
        Xp = baryon_vals[5];
        Abar = baryon_vals[6];
        Zbar = baryon_vals[7];

        mumu = tables_muon.interp(lrho, ltemp, ye, ymu, TABMUMU);
        return tables_ele.interp(lrho, ltemp, ye, ymu, TABMUELE);
    }

    void press_eps_ye_ymu__beta_eq__rho_temp(
        double& press, double& eps,
        double& ye, double& ymu,
        double& rho, double& temp, err_t& err) const
    {
        limit_rho(rho, err);
        limit_temp(temp, err);

        double const yemin = this->eos_yemin;
        double const yemax = this->eos_yemax;
        double const ymumin = eos_ymumin;
        double const ymumax = eos_ymumax;

        auto solve_ye_for_ymu = [&](double ymu_trial) {
            auto dmu_beta = [&](double ye_trial) {
                err_t local_err;
                double rhol = rho;
                double templ = temp;
                double yel = ye_trial;
                double ymul = ymu_trial;
                double mumu_dummy{}, mup{}, mun{};
                double const mue = mue_mumu_mup_mun__temp_rho_ye_ymu(
                    mumu_dummy, mup, mun, templ, rhol, yel, ymul, local_err);
                return mun - mup - mue;
            };
            return utils::brent(dmu_beta, yemin, yemax, 1e-14);
        };

        auto mu_equilibrium = [&](double ymu_trial) {
            double const ye_trial = solve_ye_for_ymu(ymu_trial);
            err_t local_err;
            double rhol = rho;
            double templ = temp;
            double yel = ye_trial;
            double ymul = ymu_trial;
            double mumu{}, mup{}, mun{};
            double const mue = mue_mumu_mup_mun__temp_rho_ye_ymu(
                mumu, mup, mun, templ, rhol, yel, ymul, local_err);
            return mumu - mue;
        };

        double const fmin = mu_equilibrium(ymumin);
        double const fmax = mu_equilibrium(ymumax);

        if (fmin * fmax > 0.0) {
            ymu = (Kokkos::fabs(fmin) < Kokkos::fabs(fmax)) ? ymumin : ymumax;
        } else {
            ymu = utils::brent(mu_equilibrium, ymumin, ymumax, 1e-14);
        }
        ye = solve_ye_for_ymu(ymu);

        err_t local_err;
        press = press__temp_rho_ye_ymu(temp, rho, ye, ymu, local_err);
        eps = eps__temp_rho_ye_ymu(temp, rho, ye, ymu, local_err);
        for (size_t iw = 0; iw < err_t::kWords; ++iw) {
            err.words[iw] |= local_err.words[iw];
        }
    }

    tabeos_linterp_4d_t tables;
    tabeos_linterp_4d_t tables_ele;
    tabeos_linterp_4d_t tables_muon;

    struct cold_eos_linterp_1d_t {
        cold_eos_linterp_1d_t() = default;

        cold_eos_linterp_1d_t(
            Kokkos::View<double **> tabs,
            Kokkos::View<double *> ar)
            : _tables(tabs), _logrho(ar)
        {
            idr = 1. / (_logrho(1) - _logrho(0));
        }

        KOKKOS_INLINE_FUNCTION double interp(double lrho, int idx) const
        {
            int const i = Kokkos::max(
                0UL,
                Kokkos::min(
                    static_cast<size_t>((lrho - _logrho(0)) * idr),
                    _logrho.extent(0) - 2));
            double const lam = (lrho - _logrho(i)) * idr;
            return (1.0 - lam) * _tables(i, idx) + lam * _tables(i + 1, idx);
        }

        Kokkos::View<const double*> _logrho;
        Kokkos::View<double **> _tables;
        double idr{};
    } cold_table;

    int nrho{}, nT{}, nye{}, nymu{};
    double energy_shift{};
    double lrhomin{}, lrhomax{};
    double ltempmin{}, ltempmax{};
    double eos_ymumin{}, eos_ymumax{};

  private:
    double _c2p_ymu_atm{};

    KOKKOS_INLINE_FUNCTION void limit_rho(double& rho, err_t& err) const
    {
        if (rho < this->eos_rhomin) {
            rho = this->eos_rhomin;
            err.set(EOS_RHO_TOO_LOW);
        }
        if (rho > this->eos_rhomax) {
            rho = this->eos_rhomax;
            err.set(EOS_RHO_TOO_HIGH);
        }
    }

    KOKKOS_INLINE_FUNCTION void limit_ye(double& ye, err_t& err) const
    {
        if (ye < this->eos_yemin) {
            ye = this->eos_yemin;
            err.set(EOS_YE_TOO_LOW);
        }
        if (ye > this->eos_yemax) {
            ye = this->eos_yemax;
            err.set(EOS_YE_TOO_HIGH);
        }
    }

    KOKKOS_INLINE_FUNCTION void limit_ymu(double& ymu, err_t& err) const
    {
        if (ymu < eos_ymumin) {
            ymu = eos_ymumin;
            err.set(EOS_YE_TOO_LOW);
        }
        if (ymu > eos_ymumax) {
            ymu = eos_ymumax;
            err.set(EOS_YE_TOO_HIGH);
        }
    }

    KOKKOS_INLINE_FUNCTION void limit_temp(double& temp, err_t& err) const
    {
        double const tmin = Kokkos::exp(ltempmin);
        double const tmax = Kokkos::exp(ltempmax);
        if (temp < tmin) {
            temp = tmin;
            err.set(EOS_TEMPERATURE_TOO_LOW);
        }
        if (temp > tmax) {
            temp = tmax;
            err.set(EOS_TEMPERATURE_TOO_HIGH);
        }
    }

    KOKKOS_INLINE_FUNCTION void
    limit_entropy_rho_ye_ymu(double& entropy, double& rho, double& ye,
                             double& ymu, err_t& err) const
    {
        double smin{}, smax{};
        entropy_range__rho_ye_ymu(smin, smax, rho, ye, ymu, err);
        if (entropy < smin) {
            entropy = smin;
            err.set(EOS_ENTROPY_TOO_LOW);
        }
        if (entropy > smax) {
            entropy = smax;
            err.set(EOS_ENTROPY_TOO_HIGH);
        }
    }

    KOKKOS_INLINE_FUNCTION double
    ltemp__eps_lrho_ye_ymu(double& eps, double lrho, double ye, double ymu,
                           err_t& err) const
    {
        double const leps = Kokkos::log(eps + energy_shift);
        double const lepsmin = Kokkos::log(
            eps__lrho_ltemp_ye_ymu(lrho, ltempmin, ye, ymu) + energy_shift);
        double const lepsmax = Kokkos::log(
            eps__lrho_ltemp_ye_ymu(lrho, ltempmax, ye, ymu) + energy_shift);
        if (leps <= lepsmin) {
            eps = Kokkos::exp(lepsmin) - energy_shift;
            err.set(EOS_EPS_TOO_LOW);
            return ltempmin;
        }
        if (leps >= lepsmax) {
            eps = Kokkos::exp(lepsmax) - energy_shift;
            err.set(EOS_EPS_TOO_HIGH);
            return ltempmax;
        }
        auto rootfun = [this, lrho, ye, ymu, leps](double lt) {
            return Kokkos::log(eps__lrho_ltemp_ye_ymu(lrho, lt, ye, ymu) + energy_shift) - leps;
        };
        return utils::brent(rootfun, ltempmin, ltempmax, 1e-14);
    }

    KOKKOS_INLINE_FUNCTION double
    ltemp__entropy_lrho_ye_ymu(double& entropy, double lrho,
                               double ye, double ymu, err_t& err) const
    {
        double const smin = entropy__lrho_ltemp_ye_ymu(lrho, ltempmin, ye, ymu);
        double const smax = entropy__lrho_ltemp_ye_ymu(lrho, ltempmax, ye, ymu);
        if (entropy <= smin) {
            entropy = smin;
            err.set(EOS_ENTROPY_TOO_LOW);
            return ltempmin;
        }
        if (entropy >= smax) {
            entropy = smax;
            err.set(EOS_ENTROPY_TOO_HIGH);
            return ltempmax;
        }
        auto rootfun = [this, lrho, ye, ymu, entropy](double lt) {
            return entropy__lrho_ltemp_ye_ymu(lrho, lt, ye, ymu) - entropy;
        };
        return utils::brent(rootfun, ltempmin, ltempmax, 1e-14);
    }

    KOKKOS_INLINE_FUNCTION double
    press__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                             double ye, double ymu) const
    {
        return Kokkos::exp(tables.interp(lrho, ltemp, ye, ymu, TABPRESS))
             + electron_pressure__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu)
             + muon_pressure__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu);
    }

    KOKKOS_INLINE_FUNCTION double
    eps__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                           double ye, double ymu) const
    {
        return Kokkos::exp(tables.interp(lrho, ltemp, ye, ymu, TABEPS)) - energy_shift
             + electron_eps__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu)
             + muon_eps__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu);
    }

    KOKKOS_INLINE_FUNCTION double
    entropy__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                               double ye, double ymu) const
    {
        return tables.interp(lrho, ltemp, ye, ymu, TABENTROPY)
             + electron_entropy__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu)
             + muon_entropy__lrho_ltemp_ye_ymu(lrho, ltemp, ye, ymu);
    }

    KOKKOS_INLINE_FUNCTION double
    electron_pressure__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                                         double ye, double ymu) const
    {
        std::array<int, 2> ids{TABPRESS_E_MINUS, TABPRESS_E_PLUS};
        std::array<double, 2> vals{};
        tables_ele.interp<2>(lrho, ltemp, ye, ymu, ids, vals);
        return vals[0] + vals[1];
    }

    KOKKOS_INLINE_FUNCTION double
    electron_eps__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                                    double ye, double ymu) const
    {
        std::array<int, 2> ids{TABEPS_E_MINUS, TABEPS_E_PLUS};
        std::array<double, 2> vals{};
        tables_ele.interp<2>(lrho, ltemp, ye, ymu, ids, vals);
        return vals[0] + vals[1];
    }

    KOKKOS_INLINE_FUNCTION double
    electron_entropy__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                                        double ye, double ymu) const
    {
        std::array<int, 2> ids{TABS_E_MINUS, TABS_E_PLUS};
        std::array<double, 2> vals{};
        tables_ele.interp<2>(lrho, ltemp, ye, ymu, ids, vals);
        return vals[0] + vals[1];
    }

    KOKKOS_INLINE_FUNCTION double
    muon_pressure__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                                     double ye, double ymu) const
    {
        std::array<int, 2> ids{TABPRESS_MU_MINUS, TABPRESS_MU_PLUS};
        std::array<double, 2> vals{};
        tables_muon.interp<2>(lrho, ltemp, ye, ymu, ids, vals);
        return vals[0] + vals[1];
    }

    KOKKOS_INLINE_FUNCTION double
    muon_eps__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                                double ye, double ymu) const
    {
        std::array<int, 2> ids{TABEPS_MU_MINUS, TABEPS_MU_PLUS};
        std::array<double, 2> vals{};
        tables_muon.interp<2>(lrho, ltemp, ye, ymu, ids, vals);
        return vals[0] + vals[1];
    }

    KOKKOS_INLINE_FUNCTION double
    muon_entropy__lrho_ltemp_ye_ymu(double lrho, double ltemp,
                                    double ye, double ymu) const
    {
        std::array<int, 2> ids{TABS_MU_MINUS, TABS_MU_PLUS};
        std::array<double, 2> vals{};
        tables_muon.interp<2>(lrho, ltemp, ye, ymu, ids, vals);
        return vals[0] + vals[1];
    }

    KOKKOS_INLINE_FUNCTION void
    thermo_press_eps_csnd2_entropy__lrho_ltemp_ye_ymu(
        double& press, double& eps, double& csnd2, double& entropy,
        double lrho, double ltemp, double ye, double ymu) const
    {
        std::array<int, 4> baryon_ids{TABPRESS, TABEPS, TABCSND2, TABENTROPY};
        std::array<double, 4> baryon_vals{};
        tables.interp<4>(lrho, ltemp, ye, ymu, baryon_ids, baryon_vals);

        std::array<int, 6> ele_ids{
            TABPRESS_E_MINUS, TABPRESS_E_PLUS,
            TABEPS_E_MINUS, TABEPS_E_PLUS,
            TABS_E_MINUS, TABS_E_PLUS};
        std::array<double, 6> ele_vals{};
        tables_ele.interp<6>(lrho, ltemp, ye, ymu, ele_ids, ele_vals);

        std::array<int, 6> mu_ids{
            TABPRESS_MU_MINUS, TABPRESS_MU_PLUS,
            TABEPS_MU_MINUS, TABEPS_MU_PLUS,
            TABS_MU_MINUS, TABS_MU_PLUS};
        std::array<double, 6> mu_vals{};
        tables_muon.interp<6>(lrho, ltemp, ye, ymu, mu_ids, mu_vals);

        press = Kokkos::exp(baryon_vals[0])
              + ele_vals[0] + ele_vals[1]
              + mu_vals[0] + mu_vals[1];
        eps = Kokkos::exp(baryon_vals[1]) - energy_shift
            + ele_vals[2] + ele_vals[3]
            + mu_vals[2] + mu_vals[3];
        csnd2 = baryon_vals[2];
        entropy = baryon_vals[3]
                + ele_vals[4] + ele_vals[5]
                + mu_vals[4] + mu_vals[5];
    }

    KOKKOS_INLINE_FUNCTION double
    cold_lpress__press_limited(double& press_cold, err_t& err) const
    {
        double const pmin = Kokkos::exp(cold_table._tables(0, CTABPRESS));
        int const n = cold_table._logrho.size();
        double const pmax = Kokkos::exp(cold_table._tables(n - 1, CTABPRESS));
        if (press_cold < pmin) {
            err.set(EOS_RHO_TOO_LOW);
            press_cold = pmin * (1.0 + 1e-10);
        }
        if (press_cold > pmax) {
            err.set(EOS_RHO_TOO_HIGH);
            press_cold = pmax * (1.0 - 1e-10);
        }
        return Kokkos::log(press_cold);
    }
};

} /* namespace grace */

#endif /* GRACE_PHYSICS_EOS_LEPTONIC_4D_HH */
