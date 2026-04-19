/**
 * @file grmhd_lep4d_id.hh
 * @brief  Initial-data helpers for the leptonic 4D EOS.
 * @date   2025
 * @copyright This file is part of GRACE. GPL-3 or later.
 */

#ifndef GRACE_GRMHD_LEP4D_ID_HH
#define GRACE_GRMHD_LEP4D_ID_HH

#ifdef GRACE_ENABLE_LEPTONIC_4D

#include <grace_config.h>
#include <grace/amr/amr_functions.hh>
#include <grace/data_structures/variable_indices.hh>
#include <grace/data_structures/variables.hh>
#include <grace/physics/eos/eos_storage_lep4d.hh>
#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/utils/grace_utils.hh>
#include <grace/utils/rootfinding.hh>

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace grace {

GRACE_HOST_DEVICE inline void
find_beta_eq_ye_ymu(double& ye_eq, double& ymu_eq,
                    double rho, double temp,
                    leptonic_eos_4d_t const& eos)
{
    rho = Kokkos::fmax(eos.density_minimum(), Kokkos::fmin(eos.density_maximum(), rho));
    temp = Kokkos::fmax(Kokkos::exp(eos.ltempmin),
                        Kokkos::fmin(Kokkos::exp(eos.ltempmax), temp));

    double const yemin = eos.get_c2p_ye_min();
    double const yemax = eos.get_c2p_ye_max();
    double const ymumin = eos.get_c2p_ymu_min();
    double const ymumax = eos.get_c2p_ymu_max();

    auto solve_ye_for_ymu = [&](double ymu_trial) {
        auto residual_inner = [&](double ye_trial) {
            eos_err_t lerr;
            double rhol{rho}, templ{temp}, yel{ye_trial}, ymul{ymu_trial};
            double mumu_dummy{}, mup{}, mun{};
            double const mue = eos.mue_mumu_mup_mun__temp_rho_ye_ymu(
                mumu_dummy, mup, mun, templ, rhol, yel, ymul, lerr);
            return mun - mup - mue;
        };
        double const fin = residual_inner(yemin);
        double const fax = residual_inner(yemax);
        if (!::isfinite(fin)) return ::isfinite(fax) ? yemax : yemin;
        if (!::isfinite(fax)) return yemin;
        if (fin * fax > 0.0) {
            return (Kokkos::fabs(fin) < Kokkos::fabs(fax)) ? yemin : yemax;
        }
        return utils::brent(residual_inner, yemin, yemax, 1e-14);
    };

    auto residual_outer = [&](double ymu_trial) {
        double const ye_trial = solve_ye_for_ymu(ymu_trial);
        eos_err_t lerr;
        double rhol{rho}, templ{temp}, yel{ye_trial}, ymul{ymu_trial};
        double mumu{}, mup_dummy{}, mun_dummy{};
        double const mue = eos.mue_mumu_mup_mun__temp_rho_ye_ymu(
            mumu, mup_dummy, mun_dummy, templ, rhol, yel, ymul, lerr);
        return mumu - mue;
    };

    double const fmin = residual_outer(ymumin);
    double const fmax = residual_outer(ymumax);
    if (!::isfinite(fmin)) {
        ymu_eq = ::isfinite(fmax) ? ymumax : ymumin;
    } else if (!::isfinite(fmax)) {
        ymu_eq = ymumin;
    } else if (fmin * fmax > 0.0) {
        ymu_eq = (Kokkos::fabs(fmin) < Kokkos::fabs(fmax)) ? ymumin : ymumax;
    } else {
        ymu_eq = utils::brent(residual_outer, ymumin, ymumax, 1e-14);
    }
    ye_eq = solve_ye_for_ymu(ymu_eq);
}

inline void
set_ymustar_from_beta_eq(leptonic_eos_4d_t const& eos,
                         variable_list_impl_t& vars)
{
    using namespace Kokkos;

    DECLARE_GRID_EXTENTS;
    auto& state = vars.getstate();
    auto& aux = vars.getaux();

    double const rho_fl = grace::get_param<double>("grmhd", "atmosphere", "rho_fl");
    double const atmo_tol = grace::get_param<double>("grmhd", "atmosphere", "atmo_tol");
    double const ye_atm = eos.ye_atmosphere();
    double const ymu_atm = eos.get_c2p_ymu_atm();

    parallel_for(
        GRACE_EXECUTION_TAG("ID", "lep4d_set_ymustar"),
        MDRangePolicy<Rank<GRACE_NSPACEDIM + 1>, default_execution_space>(
            {VEC(0, 0, 0), 0},
            {VEC(nx + 2 * ngz, ny + 2 * ngz, nz + 2 * ngz), nq}),
        KOKKOS_LAMBDA(VEC(int const& i, int const& j, int const& k), int const& q) {
            double const dens = state(VEC(i, j, k), DENS_, q);
            if (dens <= 0.0) {
                aux(VEC(i, j, k), YE_, q) = ye_atm;
                state(VEC(i, j, k), YESTAR_, q) = 0.0;
                aux(VEC(i, j, k), YMU_, q) = ymu_atm;
                state(VEC(i, j, k), YMUSTAR_, q) = 0.0;
                return;
            }

            double ye_loc = aux(VEC(i, j, k), YE_, q);
            double temp_loc = aux(VEC(i, j, k), TEMP_, q);
            double rho_loc = aux(VEC(i, j, k), RHO_, q);

            ye_loc = Kokkos::fmin(eos.get_c2p_ye_max(), Kokkos::fmax(eos.get_c2p_ye_min(), ye_loc));
            temp_loc = Kokkos::fmin(Kokkos::exp(eos.ltempmax),
                                    Kokkos::fmax(Kokkos::exp(eos.ltempmin), temp_loc));
            rho_loc = Kokkos::fmin(eos.density_maximum(),
                                   Kokkos::fmax(eos.density_minimum(), rho_loc));

            double ymu_loc = ymu_atm;
            if (rho_loc > rho_fl * (1.0 + atmo_tol)) {
                auto residual_mu = [&](double ymu_trial) {
                    eos_err_t lerr;
                    double rhol{rho_loc}, templ{temp_loc}, yel{ye_loc}, ymul{ymu_trial};
                    double const mumu = eos.mumu__temp_rho_ye_ymu(templ, rhol, yel, ymul, lerr);
                    double const mue = eos.mue__temp_rho_ye_ymu(templ, rhol, yel, ymul, lerr);
                    if (lerr.any() || !::isfinite(mumu) || !::isfinite(mue)) {
                        return std::numeric_limits<double>::quiet_NaN();
                    }
                    return mumu - mue;
                };

                double const fmin = residual_mu(eos.get_c2p_ymu_min());
                double const fmax = residual_mu(eos.get_c2p_ymu_max());
                if (!::isfinite(fmin)) {
                    ymu_loc = ::isfinite(fmax) ? eos.get_c2p_ymu_max() : ymu_atm;
                } else if (!::isfinite(fmax)) {
                    ymu_loc = eos.get_c2p_ymu_min();
                } else if (fmin == 0.0) {
                    ymu_loc = eos.get_c2p_ymu_min();
                } else if (fmax == 0.0) {
                    ymu_loc = eos.get_c2p_ymu_max();
                } else if (fmin * fmax > 0.0) {
                    ymu_loc = (Kokkos::fabs(fmin) < Kokkos::fabs(fmax))
                                  ? eos.get_c2p_ymu_min()
                                  : eos.get_c2p_ymu_max();
                } else {
                    ymu_loc = utils::brent(
                        residual_mu,
                        eos.get_c2p_ymu_min(),
                        eos.get_c2p_ymu_max(),
                        1e-14);
                }
            } else {
                ye_loc = ye_atm;
                ymu_loc = ymu_atm;
                aux(VEC(i, j, k), YE_, q) = ye_loc;
                state(VEC(i, j, k), YESTAR_, q) = dens * ye_loc;
            }

            aux(VEC(i, j, k), YMU_, q) = ymu_loc;
            state(VEC(i, j, k), YMUSTAR_, q) = dens * ymu_loc;
        });

    fence();
}

void set_grmhd_initial_data_lep4d();

} /* namespace grace */

#endif /* GRACE_ENABLE_LEPTONIC_4D */
#endif /* GRACE_GRMHD_LEP4D_ID_HH */
