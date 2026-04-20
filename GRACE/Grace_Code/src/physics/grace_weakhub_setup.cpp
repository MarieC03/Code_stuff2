/**
 * @file grace_weakhub_setp.cpp
 * @author Marie Cassing (mcassing@itp.uni-frankfurt.de)
 * @brief Weakhub library setup for rates
 * @date 2026-03-02
 * 
 * @copyright This file is part of of the General Relativistic Astrophysics
 * Code for Exascale.
 * GRACE is an evolution framework that uses Finite Volume
 * methods to simulate relativistic spacetimes and plasmas
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
 * 
 */

#include "grace_weakhub_table.hh"

#include <stdexcept>
#include <string>
#include <vector>

#ifdef GRACE_USE_HDF5
// h5 api?
#define H5_USE_16_API 1
#include <hdf5.h>
#endif

namespace grace { namespace weakhub {
namespace {
  device_handle g_handle;
  bool g_initialized = false;

  void validate_weakhub_runtime_configuration() {
    const auto species_mode = grace::get_param<std::string>("m1", "eas", "species_mode");
    const auto eos_type = grace::get_param<std::string>("eos", "eos_type");

    if (eos_type != "tabulated" && eos_type != "leptonic_4d") {
      throw std::runtime_error(
          "Weakhub neutrino rates require a Ye-dependent tabulated EOS or the leptonic_4d EOS.");
    }

    if (eos_type == "leptonic_4d") {
      const bool use_muonic_eos =
          grace::get_param<bool>("eos", "leptonic_4d", "use_muonic_eos");
#ifdef M1_NU_FIVESPECIES
      if (!use_muonic_eos) {
        throw std::runtime_error(
            "The 5-species Weakhub path requires eos.leptonic_4d.use_muonic_eos = true.");
      }
#else
      if (use_muonic_eos) {
        throw std::runtime_error(
            "Muonic EOS contributions require the 5-species M1 build. "
            "Disable eos.leptonic_4d.use_muonic_eos or use a 5-species build.");
      }
#endif
    }

#ifdef M1_NU_THREESPECIES
    if (species_mode != "build_default" && species_mode != "three") {
      throw std::runtime_error(
          "This GRACE executable was built with three M1 neutrino species. "
          "Set m1.eas.species_mode to \"three\" or \"build_default\".");
    }
#elif defined(M1_NU_FIVESPECIES)
    if (species_mode != "build_default" && species_mode != "five") {
      throw std::runtime_error(
          "This GRACE executable was built with five M1 neutrino species. "
          "Set m1.eas.species_mode to \"five\" or \"build_default\".");
    }
    if (eos_type != "leptonic_4d") {
      throw std::runtime_error(
          "The 5-species Weakhub path requires eos.eos_type = \"leptonic_4d\".");
    }
#else
    if (species_mode != "build_default" && species_mode != "one") {
      throw std::runtime_error(
          "This GRACE executable was built with one M1 neutrino species. "
          "Set m1.eas.species_mode to \"one\" or \"build_default\".");
    }
#endif
  }

#ifdef GRACE_USE_HDF5
  template <typename T>
  void read_dataset(hid_t file, const char* name, hid_t type, T* data) {
    hid_t dset = H5Dopen(file, name);
    if (dset < 0) throw std::runtime_error(std::string("Failed to open HDF5 dataset: ") + name);
    const herr_t ierr = H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset);
    if (ierr < 0) throw std::runtime_error(std::string("Failed to read HDF5 dataset: ") + name);
  }
#endif
}

bool weakhub_enabled_from_params() {
  bool use_weakhub = false;
  try { use_weakhub = grace::get_param<bool>("m1","eas","use_weakhub"); } catch(...) {}
  try {
    const auto kind = grace::get_param<std::string>("m1","eas","kind");
    if (kind == "neutrino_weakhub") use_weakhub = true;
    if (kind == "neutrino_analytic") use_weakhub = false;
  } catch(...) {}
  bool use_analytic = false;
  try { use_analytic = grace::get_param<bool>("m1","eas","use_analytic"); } catch(...) {}
  if (use_analytic) use_weakhub = false;
  return use_weakhub;
}

void initialize_weakhub_from_params() {
  if (g_initialized || !weakhub_enabled_from_params()) return;
#ifndef GRACE_USE_HDF5
  throw std::runtime_error("Weakhub requested but GRACE_USE_HDF5 is not enabled.");
#else
  validate_weakhub_runtime_configuration();
  const auto path = grace::get_param<std::string>("m1","eas","weakhub_table");
  hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) throw std::runtime_error("Failed to open Weakhub HDF5 table: " + path);

  int n_species = 0;
  int IVrho = 0, IVtemp = 0, IVye = 0, IVymu = 0;
  read_dataset(file, "n_species", H5T_NATIVE_INT, &n_species);
  read_dataset(file, "IVrho", H5T_NATIVE_INT, &IVrho);
  read_dataset(file, "IVtemp", H5T_NATIVE_INT, &IVtemp);
  read_dataset(file, "IVye", H5T_NATIVE_INT, &IVye);
  read_dataset(file, "IVymu", H5T_NATIVE_INT, &IVymu);

#ifdef M1_NU_THREESPECIES
  if (n_species != 3) {
    H5Fclose(file);
    throw std::runtime_error(
        "The current GRACE build uses 3 M1 neutrino species, but the selected weakhub table "
        "does not provide the required 3-spec layout.");
  }
#elif defined(M1_NU_FIVESPECIES)
  if (n_species != 6) {
    H5Fclose(file);
    throw std::runtime_error(
        "The current GRACE build uses 5 M1 neutrino species, but the selected weakhub table "
        "does not provide the required 6-spec muonic layout.");
  }
  if (IVymu <= 0) {
    H5Fclose(file);
    throw std::runtime_error(
        "The current GRACE build uses 5 M1 neutrino species, but the selected weakhub table "
        "has no Y_mu dimension.");
  }
#endif

  std::vector<double> logrho(IVrho), logtemp(IVtemp), ye(IVye), logymu(std::max(IVymu,1));
  read_dataset(file, "logrho_IVtable", H5T_NATIVE_DOUBLE, logrho.data());
  read_dataset(file, "logtemp_IVtable", H5T_NATIVE_DOUBLE, logtemp.data());
  read_dataset(file, "ye_IVtable", H5T_NATIVE_DOUBLE, ye.data());
  if (IVymu > 0) read_dataset(file, "logymu_IVtable", H5T_NATIVE_DOUBLE, logymu.data());

  const size_t table_size = static_cast<size_t>(IVrho) * IVtemp * IVye * std::max(IVymu,1) * n_species;
  std::vector<double> raw_ae(table_size), raw_an(table_size), raw_s(table_size);
  read_dataset(file, "kappa_a_en_grey_table", H5T_NATIVE_DOUBLE, raw_ae.data());
  read_dataset(file, "kappa_a_num_grey_table", H5T_NATIVE_DOUBLE, raw_an.data());
  read_dataset(file, "kappa_s_grey_table", H5T_NATIVE_DOUBLE, raw_s.data());
  H5Fclose(file);

  std::vector<double> ord_ae(table_size), ord_an(table_size), ord_s(table_size);
  for (int iv = 0; iv < n_species; ++iv)
    for (int l = 0; l < std::max(IVymu,1); ++l)
      for (int k = 0; k < IVye; ++k)
        for (int j = 0; j < IVtemp; ++j)
          for (int i = 0; i < IVrho; ++i) {
            const size_t indold = i + IVrho * (j + IVtemp * (k + IVye * (l + std::max(IVymu,1) * iv)));
            const size_t indnew = iv + n_species * (i + IVrho * (j + IVtemp * (k + IVye * l)));
            ord_ae[indnew] = raw_ae[indold];
            ord_an[indnew] = raw_an[indold];
            ord_s[indnew]  = raw_s[indold];
          }

  g_handle.n_species_table = n_species;
  g_handle.nrho = IVrho; g_handle.ntemp = IVtemp; g_handle.nye = IVye; g_handle.nymu = std::max(IVymu,1);
  g_handle.logrho_min = logrho.front(); g_handle.logrho_max = logrho.back();
  g_handle.logtemp_min = logtemp.front(); g_handle.logtemp_max = logtemp.back();
  g_handle.ye_min = ye.front(); g_handle.ye_max = ye.back();
  g_handle.logymu_min = logymu.front(); g_handle.logymu_max = logymu.back();

  g_handle.logrho_axis = Kokkos::View<double*>("weakhub_logrho", IVrho);
  g_handle.logtemp_axis = Kokkos::View<double*>("weakhub_logtemp", IVtemp);
  g_handle.ye_axis = Kokkos::View<double*>("weakhub_ye", IVye);
  g_handle.logymu_axis = Kokkos::View<double*>("weakhub_logymu", std::max(IVymu,1));
  g_handle.kappa_a_en_table = Kokkos::View<double*>("weakhub_ae", table_size);
  g_handle.kappa_a_num_table = Kokkos::View<double*>("weakhub_an", table_size);
  g_handle.kappa_s_table = Kokkos::View<double*>("weakhub_s", table_size);

  auto h_rho = Kokkos::create_mirror_view(g_handle.logrho_axis);
  auto h_temp = Kokkos::create_mirror_view(g_handle.logtemp_axis);
  auto h_ye = Kokkos::create_mirror_view(g_handle.ye_axis);
  auto h_ymu = Kokkos::create_mirror_view(g_handle.logymu_axis);
  auto h_ae = Kokkos::create_mirror_view(g_handle.kappa_a_en_table);
  auto h_an = Kokkos::create_mirror_view(g_handle.kappa_a_num_table);
  auto h_s = Kokkos::create_mirror_view(g_handle.kappa_s_table);
  for (int i = 0; i < IVrho; ++i) h_rho(i) = logrho[i];
  for (int i = 0; i < IVtemp; ++i) h_temp(i) = logtemp[i];
  for (int i = 0; i < IVye; ++i) h_ye(i) = ye[i];
  for (int i = 0; i < std::max(IVymu,1); ++i) h_ymu(i) = logymu[i];
  for (size_t i = 0; i < table_size; ++i) { h_ae(i) = ord_ae[i]; h_an(i) = ord_an[i]; h_s(i) = ord_s[i]; }
  Kokkos::deep_copy(g_handle.logrho_axis, h_rho);
  Kokkos::deep_copy(g_handle.logtemp_axis, h_temp);
  Kokkos::deep_copy(g_handle.ye_axis, h_ye);
  Kokkos::deep_copy(g_handle.logymu_axis, h_ymu);
  Kokkos::deep_copy(g_handle.kappa_a_en_table, h_ae);
  Kokkos::deep_copy(g_handle.kappa_a_num_table, h_an);
  Kokkos::deep_copy(g_handle.kappa_s_table, h_s);
  g_handle.valid = true;
  g_initialized = true;
#endif
}

const device_handle& get_device_handle() { return g_handle; }
bool is_initialized() { return g_initialized && g_handle.valid; }

}} // namespace grace::weakhub
