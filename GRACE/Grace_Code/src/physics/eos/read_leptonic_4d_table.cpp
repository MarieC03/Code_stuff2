/**
 * @file read_leptonic_4d_table.cpp
 * @brief  Implementation of the 4D leptonic EOS table reader for GRACE.
 *         Reads the native GRACE leptonic 4D HDF5 format and converts the
 *         FIL/Margherita-derived baryon, electron, and muon sectors into the
 *         standard GRACE EOS objects.
 * @date   2025
 *
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/eos/read_leptonic_4d_table.hh>
#include <grace/physics/eos/tabulated_eos.hh>
#include <grace/system/grace_system.hh>
#include <grace/utils/grace_utils.hh>

#define H5_USE_16_API 1
#include <hdf5.h>

#include <Kokkos_Core.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#define HDF5_CALL(ret, fn_call) \
    do { (ret) = (fn_call); ASSERT((ret) >= 0, "HDF5 call failed: " #fn_call); } while (0)

#define READ_4D_EOS_HDF5(NAME, VAR, TYPE, MEM) \
    do { \
        hid_t _ds; \
        HDF5_CALL(_ds, H5Dopen(file, NAME, H5P_DEFAULT)); \
        HDF5_CALL(h5err, H5Dread(_ds, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR)); \
        HDF5_CALL(h5err, H5Dclose(_ds)); \
    } while (0)

#define READ_4D_EOSTABLE_HDF5(NAME, OFF) \
    do { \
        hsize_t offset5[2] = {OFF, 0}; \
        H5Sselect_hyperslab(mem5, H5S_SELECT_SET, offset5, NULL, var5, NULL); \
        READ_4D_EOS_HDF5(NAME, alltables_temp, H5T_NATIVE_DOUBLE, mem5); \
    } while (0)

namespace grace {

static void read_leptonic_cold_table(
    const std::string& filename,
    Kokkos::View<double**, grace::default_execution_space>& d_data,
    Kokkos::View<double*, grace::default_execution_space>& d_rho)
{
    std::ifstream f(filename);
    ASSERT(f.is_open(), "Cannot open leptonic cold table: " << filename);

    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    std::istringstream iss_n(line);
    size_t nrows{};
    iss_n >> nrows;
    ASSERT(iss_n, "Failed to read nrows in leptonic cold table");

    constexpr int NCOLS = leptonic_eos_4d_t::COLD_VIDX::N_CTAB_VARS;

    Kokkos::realloc(d_data, nrows, NCOLS);
    Kokkos::realloc(d_rho, nrows);
    auto h_data = Kokkos::create_mirror_view(d_data);
    auto h_rho = Kokkos::create_mirror_view(d_rho);

    for (size_t i = 0; i < nrows; ++i) {
        ASSERT(std::getline(f, line), "Unexpected EOF in leptonic cold table at row " << i);
        std::istringstream iss(line);
        double logrho_val{};
        iss >> logrho_val;
        h_rho(i) = logrho_val;
        for (int j = 0; j < NCOLS; ++j) {
            double val{};
            iss >> val;
            h_data(i, j) = val;
        }
    }

    Kokkos::deep_copy(d_data, h_data);
    Kokkos::deep_copy(d_rho, h_rho);
}

grace::leptonic_eos_4d_t read_leptonic_4d_table()
{
    using namespace grace::physical_constants;

    auto const fname = grace::get_param<std::string>("eos", "leptonic_4d", "table_filename");
    auto const cold_fname = grace::get_param<std::string>("eos", "leptonic_4d", "cold_table_filename");
    auto const tab_format = grace::get_param<std::string>("eos", "leptonic_4d", "table_format");

    bool const do_energy_shift =
        grace::get_param<bool>("eos", "leptonic_4d", "do_energy_shift");
    bool const use_muonic_eos =
        grace::get_param<bool>("eos", "leptonic_4d", "use_muonic_eos");
    bool const atm_beta_eq =
        grace::get_param<bool>("eos", "leptonic_4d", "atmosphere_beta_eq");

    GRACE_INFO("Reading 4D leptonic EOS table: {}", fname);
    GRACE_INFO("Format: {}", tab_format);
    ASSERT(tab_format == "leptonic_4d_native",
           "Unsupported leptonic_4d table_format: " << tab_format);

    herr_t h5err;
    hid_t file;
    HDF5_CALL(file, H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));

    int nrho{}, ntemp{}, nye{}, nymu{};
    READ_4D_EOS_HDF5("pointsrho", &nrho, H5T_NATIVE_INT, H5S_ALL);
    READ_4D_EOS_HDF5("pointstemp", &ntemp, H5T_NATIVE_INT, H5S_ALL);
    READ_4D_EOS_HDF5("pointsye", &nye, H5T_NATIVE_INT, H5S_ALL);
    READ_4D_EOS_HDF5("pointsymu", &nymu, H5T_NATIVE_INT, H5S_ALL);

    GRACE_INFO("4D EOS grid: nrho={} nT={} nye={} nymu={}", nrho, ntemp, nye, nymu);

    size_t const NTAB_B = leptonic_eos_4d_t::TEOS_VIDX::N_TAB_VARS_BARYON;
    size_t const NTAB_E = leptonic_eos_4d_t::ELE_VIDX::N_TAB_VARS_ELE;
    size_t const NTAB_M = leptonic_eos_4d_t::MUON_VIDX::N_TAB_VARS_MUON;
    size_t const NP4 = static_cast<size_t>(nrho) * ntemp * nye * nymu;

    auto* alltables_temp = new double[std::max({NTAB_B, NTAB_E, NTAB_M}) * NP4];
    hsize_t table_dims[2] = {std::max({NTAB_B, NTAB_E, NTAB_M}), NP4};
    hsize_t var5[2] = {1, NP4};
    hid_t mem5 = H5Screate_simple(2, table_dims, NULL);

    std::vector<double> logrho(nrho), logtemp(ntemp), yes(nye), ymus(nymu);
    double energy_shift_raw{0.0};
    double baryon_mass = mn_MeV * MeV_to_g;

    READ_4D_EOS_HDF5("logrho", logrho.data(), H5T_NATIVE_DOUBLE, H5S_ALL);
    READ_4D_EOS_HDF5("logtemp", logtemp.data(), H5T_NATIVE_DOUBLE, H5S_ALL);
    READ_4D_EOS_HDF5("ye", yes.data(), H5T_NATIVE_DOUBLE, H5S_ALL);
    READ_4D_EOS_HDF5("ymu", ymus.data(), H5T_NATIVE_DOUBLE, H5S_ALL);
    READ_4D_EOS_HDF5("energy_shift", &energy_shift_raw, H5T_NATIVE_DOUBLE, H5S_ALL);

    if (H5Lexists(file, "/mass_factor", H5P_DEFAULT)) {
        hid_t mb_ds;
        HDF5_CALL(mb_ds, H5Dopen(file, "mass_factor", H5P_DEFAULT));
        H5Dread(mb_ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &baryon_mass);
        GRACE_INFO("Baryon mass from file: {}", baryon_mass);
        HDF5_CALL(h5err, H5Dclose(mb_ds));
    }

    constexpr double G_si = physical_constants::G_si;
    constexpr double c_si = physical_constants::c_si;
    constexpr double Msun = physical_constants::Msun_si;
    constexpr double LENGTHGF = G_si * Msun / (c_si * c_si);
    constexpr double TIMEGF = LENGTHGF / c_si;
    constexpr double RHOGF = Msun / (LENGTHGF * LENGTHGF * LENGTHGF);
    constexpr double PRESSGF = 1.0 / (RHOGF * c_si * c_si);
    constexpr double EPSGF = 1.0 / (c_si * c_si);
    constexpr double MeV2geom = physical_constants::MeV_to_g;

    double const energy_shift = do_energy_shift ? energy_shift_raw * EPSGF : 0.0;

    for (int i = 0; i < nrho; ++i) {
        logrho[i] = logrho[i] * std::log(10.0) + std::log(RHOGF);
    }
    for (int i = 0; i < ntemp; ++i) {
        logtemp[i] = logtemp[i] * std::log(10.0);
    }

    using view5d = Kokkos::View<double*****, grace::default_execution_space>;
    view5d v_baryon("eos4d_baryon", nrho, ntemp, nye, nymu, NTAB_B);
    view5d v_ele("eos4d_ele", nrho, ntemp, nye, nymu, NTAB_E);
    view5d v_muon("eos4d_muon", nrho, ntemp, nye, nymu, NTAB_M);

    auto h_baryon = Kokkos::create_mirror_view(v_baryon);
    auto h_ele = Kokkos::create_mirror_view(v_ele);
    auto h_muon = Kokkos::create_mirror_view(v_muon);
    Kokkos::deep_copy(h_muon, 0.0);

    auto fill_view = [&](auto& hv, size_t ntab) {
        for (size_t iv = 0; iv < ntab; ++iv) {
            for (int lm = 0; lm < nymu; ++lm) {
                for (int k = 0; k < nye; ++k) {
                    for (int j = 0; j < ntemp; ++j) {
                        for (int i = 0; i < nrho; ++i) {
                            size_t const indold =
                                i + nrho * (j + ntemp * (k + nye * (lm + nymu * iv)));
                            hv(i, j, k, lm, iv) = alltables_temp[indold];
                        }
                    }
                }
            }
        }
    };

    table_dims[0] = NTAB_B;
    H5Sset_extent_simple(mem5, 2, table_dims, NULL);
    READ_4D_EOSTABLE_HDF5("logpress", leptonic_eos_4d_t::TEOS_VIDX::TABPRESS);
    READ_4D_EOSTABLE_HDF5("logenergy", leptonic_eos_4d_t::TEOS_VIDX::TABEPS);
    READ_4D_EOSTABLE_HDF5("entropy", leptonic_eos_4d_t::TEOS_VIDX::TABENTROPY);
    READ_4D_EOSTABLE_HDF5("cs2", leptonic_eos_4d_t::TEOS_VIDX::TABCSND2);
    READ_4D_EOSTABLE_HDF5("mu_e", leptonic_eos_4d_t::TEOS_VIDX::TABMUE);
    READ_4D_EOSTABLE_HDF5("mu_p", leptonic_eos_4d_t::TEOS_VIDX::TABMUP);
    READ_4D_EOSTABLE_HDF5("mu_n", leptonic_eos_4d_t::TEOS_VIDX::TABMUN);
    READ_4D_EOSTABLE_HDF5("Xa", leptonic_eos_4d_t::TEOS_VIDX::TABXA);
    READ_4D_EOSTABLE_HDF5("Xh", leptonic_eos_4d_t::TEOS_VIDX::TABXH);
    READ_4D_EOSTABLE_HDF5("Xn", leptonic_eos_4d_t::TEOS_VIDX::TABXN);
    READ_4D_EOSTABLE_HDF5("Xp", leptonic_eos_4d_t::TEOS_VIDX::TABXP);
    READ_4D_EOSTABLE_HDF5("Abar", leptonic_eos_4d_t::TEOS_VIDX::TABABAR);
    READ_4D_EOSTABLE_HDF5("Zbar", leptonic_eos_4d_t::TEOS_VIDX::TABZBAR);
    fill_view(h_baryon, NTAB_B);

    table_dims[0] = NTAB_E;
    H5Sset_extent_simple(mem5, 2, table_dims, NULL);
    READ_4D_EOSTABLE_HDF5("ele/mu_ele", leptonic_eos_4d_t::ELE_VIDX::TABMUELE);
    READ_4D_EOSTABLE_HDF5("ele/yle_minus", leptonic_eos_4d_t::ELE_VIDX::TABYLE_MINUS);
    READ_4D_EOSTABLE_HDF5("ele/yle_plus", leptonic_eos_4d_t::ELE_VIDX::TABYLE_PLUS);
    READ_4D_EOSTABLE_HDF5("ele/press_e_minus", leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_MINUS);
    READ_4D_EOSTABLE_HDF5("ele/press_e_plus", leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_PLUS);
    READ_4D_EOSTABLE_HDF5("ele/eps_e_minus", leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_MINUS);
    READ_4D_EOSTABLE_HDF5("ele/eps_e_plus", leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_PLUS);
    READ_4D_EOSTABLE_HDF5("ele/s_e_minus", leptonic_eos_4d_t::ELE_VIDX::TABS_E_MINUS);
    READ_4D_EOSTABLE_HDF5("ele/s_e_plus", leptonic_eos_4d_t::ELE_VIDX::TABS_E_PLUS);
    fill_view(h_ele, NTAB_E);

    if (use_muonic_eos) {
        table_dims[0] = NTAB_M;
        H5Sset_extent_simple(mem5, 2, table_dims, NULL);
        READ_4D_EOSTABLE_HDF5("muon/mu_mu", leptonic_eos_4d_t::MUON_VIDX::TABMUMU);
        READ_4D_EOSTABLE_HDF5("muon/ymu_minus", leptonic_eos_4d_t::MUON_VIDX::TABYMU_MINUS);
        READ_4D_EOSTABLE_HDF5("muon/ymu_plus", leptonic_eos_4d_t::MUON_VIDX::TABYMU_PLUS);
        READ_4D_EOSTABLE_HDF5("muon/press_mu_minus", leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_MINUS);
        READ_4D_EOSTABLE_HDF5("muon/press_mu_plus", leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_PLUS);
        READ_4D_EOSTABLE_HDF5("muon/eps_mu_minus", leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_MINUS);
        READ_4D_EOSTABLE_HDF5("muon/eps_mu_plus", leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_PLUS);
        READ_4D_EOSTABLE_HDF5("muon/s_mu_minus", leptonic_eos_4d_t::MUON_VIDX::TABS_MU_MINUS);
        READ_4D_EOSTABLE_HDF5("muon/s_mu_plus", leptonic_eos_4d_t::MUON_VIDX::TABS_MU_PLUS);
        fill_view(h_muon, NTAB_M);
    }

    H5Sclose(mem5);
    HDF5_CALL(h5err, H5Fclose(file));
    delete[] alltables_temp;

    double hmin{std::numeric_limits<double>::max()};
    double hmax{std::numeric_limits<double>::lowest()};
    double epsmin{std::numeric_limits<double>::max()};
    double epsmax_val{std::numeric_limits<double>::lowest()};

    for (int lm = 0; lm < nymu; ++lm) {
        for (int k = 0; k < nye; ++k) {
            for (int j = 0; j < ntemp; ++j) {
                for (int i = 0; i < nrho; ++i) {
                    h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABPRESS) =
                        h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABPRESS) * std::log(10.0)
                      + std::log(PRESSGF);

                    double const leps_raw =
                        h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABEPS);
                    double const epsT = std::pow(10.0, leps_raw) * EPSGF;
                    h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABEPS) = std::log(epsT);

                    double& cs2 = h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABCSND2);
                    cs2 *= LENGTHGF * LENGTHGF / (TIMEGF * TIMEGF);
                    if (cs2 < 0.0) {
                        cs2 = 0.0;
                    }

                    double const mb = baryon_mass;
                    h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABMUE) *= MeV2geom / mb;
                    h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABMUP) *= MeV2geom / mb;
                    h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABMUN) *= MeV2geom / mb;
                }
            }
        }
    }

    auto convert_lepton_view = [&](auto& hv, bool is_muon_table) {
        for (int lm = 0; lm < nymu; ++lm) {
            for (int k = 0; k < nye; ++k) {
                for (int j = 0; j < ntemp; ++j) {
                    for (int i = 0; i < nrho; ++i) {
                        double const mb = baryon_mass;
                        if (is_muon_table) {
                            hv(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABMUMU) *= MeV2geom / mb;
                            hv(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_MINUS) *= PRESSGF;
                            hv(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_PLUS) *= PRESSGF;
                            hv(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_MINUS) *= EPSGF;
                            hv(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_PLUS) *= EPSGF;
                        } else {
                            hv(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABMUELE) *= MeV2geom / mb;
                            hv(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_MINUS) *= PRESSGF;
                            hv(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_PLUS) *= PRESSGF;
                            hv(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_MINUS) *= EPSGF;
                            hv(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_PLUS) *= EPSGF;
                        }
                    }
                }
            }
        }
    };
    convert_lepton_view(h_ele, false);
    convert_lepton_view(h_muon, true);

    for (int lm = 0; lm < nymu; ++lm) {
        for (int k = 0; k < nye; ++k) {
            for (int j = 0; j < ntemp; ++j) {
                for (int i = 0; i < nrho; ++i) {
                    double const rhoL = std::exp(logrho[i]);
                    double const press_b =
                        std::exp(h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABPRESS));
                    double const eps_b =
                        std::exp(h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABEPS)) - energy_shift;
                    double const press_e =
                        h_ele(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_MINUS) +
                        h_ele(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_PLUS);
                    double const eps_e =
                        h_ele(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_MINUS) +
                        h_ele(i, j, k, lm, leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_PLUS);
                    double const press_mu =
                        h_muon(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_MINUS) +
                        h_muon(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_PLUS);
                    double const eps_mu =
                        h_muon(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_MINUS) +
                        h_muon(i, j, k, lm, leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_PLUS);

                    double const pressL = press_b + press_e + press_mu;
                    double const epsL = eps_b + eps_e + eps_mu;
                    double const hL = 1.0 + epsL + pressL / rhoL;

                    epsmin = std::min(epsmin, epsL);
                    epsmax_val = std::max(epsmax_val, epsL);
                    hmin = std::min(hmin, hL);
                    hmax = std::max(hmax, hL);

                    double& cs2 = h_baryon(i, j, k, lm, leptonic_eos_4d_t::TEOS_VIDX::TABCSND2);
                    cs2 /= hL;
                    if (cs2 > 0.9999999) {
                        cs2 = 0.9999999;
                    }
                }
            }
        }
    }

    Kokkos::deep_copy(v_baryon, h_baryon);
    Kokkos::deep_copy(v_ele, h_ele);
    Kokkos::deep_copy(v_muon, h_muon);

    using view1d = Kokkos::View<double*, grace::default_execution_space>;
    view1d v_logrho("eos4d_logrho", nrho);
    view1d v_logT("eos4d_logT", ntemp);
    view1d v_ye("eos4d_ye", nye);
    view1d v_ymu("eos4d_ymu", nymu);
    auto h_lr = Kokkos::create_mirror_view(v_logrho);
    auto h_lt = Kokkos::create_mirror_view(v_logT);
    auto h_ye = Kokkos::create_mirror_view(v_ye);
    auto h_ym = Kokkos::create_mirror_view(v_ymu);
    for (int i = 0; i < nrho; ++i) {
        h_lr(i) = logrho[i];
    }
    for (int i = 0; i < ntemp; ++i) {
        h_lt(i) = logtemp[i];
    }
    for (int i = 0; i < nye; ++i) {
        h_ye(i) = yes[i];
    }
    for (int i = 0; i < nymu; ++i) {
        h_ym(i) = ymus[i];
    }
    Kokkos::deep_copy(v_logrho, h_lr);
    Kokkos::deep_copy(v_logT, h_lt);
    Kokkos::deep_copy(v_ye, h_ye);
    Kokkos::deep_copy(v_ymu, h_ym);

    Kokkos::View<double**, grace::default_execution_space> cold_tabs;
    Kokkos::View<double*, grace::default_execution_space> cold_lrho;
    read_leptonic_cold_table(cold_fname, cold_tabs, cold_lrho);

    double const rhomax = std::exp(logrho[nrho - 1]);
    double const rhomin = std::exp(logrho[0]);
    double const tempmax = std::exp(logtemp[ntemp - 1]);
    double const tempmin = std::exp(logtemp[0]);
    double const yemax = yes[nye - 1];
    double const yemin = yes[0];
    double const ymumax = ymus[nymu - 1];
    double const ymumin = ymus[0];

    auto usr_epsmax = grace::get_param<double>("eos", "eps_maximum");
    if (usr_epsmax < epsmax_val) {
        epsmax_val = usr_epsmax;
    }

    double temp_floor = grace::get_param<double>("grmhd", "atmosphere", "temp_fl");
    double rho_floor = grace::get_param<double>("grmhd", "atmosphere", "rho_fl");
    if (temp_floor < tempmin) {
        GRACE_WARN("Requested leptonic atmosphere temperature is below table bound {}.", tempmin);
        temp_floor = tempmin * (1.0 + 1e-5);
    }
    if (rho_floor < rhomin) {
        ERROR("Requested leptonic atmosphere density is below table bound.");
    }

    double ye_atm = yemin;
    double ymu_atm = ymumin;

    GRACE_INFO("4D leptonic EOS rho [{:.4e}, {:.4e}]  T [{:.4e}, {:.4e}]  "
               "Ye [{:.3f}, {:.3f}]  Ymu [{:.3f}, {:.3f}]",
               rhomin, rhomax, tempmin, tempmax, yemin, yemax, ymumin, ymumax);

    auto eos = leptonic_eos_4d_t(
        v_baryon, v_logrho, v_logT, v_ye, v_ymu,
        v_ele, v_muon,
        cold_tabs, cold_lrho,
        rhomax, rhomin,
        tempmax, tempmin,
        yemax, yemin,
        ymumax, ymumin,
        baryon_mass, energy_shift,
        epsmin, epsmax_val,
        hmin, hmax,
        temp_floor,
        ye_atm, ymu_atm,
        atm_beta_eq);

    if (atm_beta_eq) {
        GRACE_INFO("Atmosphere: finding beta-equilibrium Y_e, Y_mu at rho_fl, T_fl");

        tabeos_linterp_4d_t host_baryon(h_baryon, h_lr, h_lt, h_ye, h_ym);
        tabeos_linterp_4d_t host_ele(h_ele, h_lr, h_lt, h_ye, h_ym);
        tabeos_linterp_4d_t host_muon(h_muon, h_lr, h_lt, h_ye, h_ym);

        auto solve_ye_for_ymu = [&](double ymu_trial) {
            auto dmu_beta = [&](double ye_trial) {
                double const lrhoL = std::log(rho_floor);
                double const ltempL = std::log(temp_floor);
                double const mue = host_ele.interp(
                    lrhoL, ltempL, ye_trial, ymu_trial,
                    leptonic_eos_4d_t::ELE_VIDX::TABMUELE);
                double const mup = host_baryon.interp(
                    lrhoL, ltempL, ye_trial, ymu_trial,
                    leptonic_eos_4d_t::TEOS_VIDX::TABMUP);
                double const mun = host_baryon.interp(
                    lrhoL, ltempL, ye_trial, ymu_trial,
                    leptonic_eos_4d_t::TEOS_VIDX::TABMUN);
                return mun - mup - mue;
            };
            return utils::brent(dmu_beta, yemin, yemax, 1e-14);
        };

        auto mu_equilibrium = [&](double ymu_trial) {
            double const ye_trial = solve_ye_for_ymu(ymu_trial);
            double const lrhoL = std::log(rho_floor);
            double const ltempL = std::log(temp_floor);
            double const mue = host_ele.interp(
                lrhoL, ltempL, ye_trial, ymu_trial,
                leptonic_eos_4d_t::ELE_VIDX::TABMUELE);
            double const mumu = host_muon.interp(
                lrhoL, ltempL, ye_trial, ymu_trial,
                leptonic_eos_4d_t::MUON_VIDX::TABMUMU);
            return mumu - mue;
        };

        double const fmin = mu_equilibrium(ymumin);
        double const fmax = mu_equilibrium(ymumax);
        if (fmin * fmax > 0.0) {
            ymu_atm = (std::fabs(fmin) < std::fabs(fmax)) ? ymumin : ymumax;
        } else {
            ymu_atm = utils::brent(mu_equilibrium, ymumin, ymumax, 1e-14);
        }
        ye_atm = solve_ye_for_ymu(ymu_atm);
        eos.set_atmosphere_fractions(ye_atm, ymu_atm);
    }

    return eos;
}

} /* namespace grace */
