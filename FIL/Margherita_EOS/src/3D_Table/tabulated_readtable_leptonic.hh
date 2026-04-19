//
//  Copyright (C) 2021, Harry Ho-Yin Ng
//  Based on routines by Elias Roland Most
//                       Ludwig Jens Papenfort
//

// All lepton_tables are already in code unit for rho, prs, ent, eps
// and MeV for Temp, mu
// mu_l is in MeV defined as mu_l = mu_0 + m_l which includes rest mass of lepton 
// Yl_plus/minus is in dimenionless
// T, rho, Ymu is in log spacing, Yle is in unit spacing
  // Indices of lepton tables:
  // 1. mu_l, 2. yl_minus, 3. yl_plus
  // 4. prs_l_minus, 5. prs_l_plus, 6. eps_l_minus,
  // 7. eps_l_plus, 8. s_l_minus, 9. s_l_plus
// logrho, logT are using natural log


#include "tabulated.hh"
#include <fstream>
#include <iostream>

#define H5_USE_16_API 1
#include <hdf5.h>

#ifndef EOS_TABULATED_READTABLE_LEPTONIC_HH
#define EOS_TABULATED_READTABLE_LEPTONIC_HH

// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_EOS_HDF5(NAME, VAR, TYPE, MEM)                             \
  do {                                                                  \
    hid_t dataset;                                                      \
    HDF5_ERROR(dataset = H5Dopen(file, NAME));                          \
    HDF5_ERROR(H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR)); \
    HDF5_ERROR(H5Dclose(dataset));                                      \
  } while (0)
// The second reads a given variable into a hyperslab of the alltables_temp
// array
#define READ_EOSTABLE_HDF5(NAME, OFF)                                    \
  do {                                                                   \
    hsize_t offset[2] = {OFF, 0};                                        \
    H5Sselect_hyperslab(mem3, H5S_SELECT_SET, offset, NULL, var3, NULL); \
    READ_EOS_HDF5(NAME, alltables_temp, H5T_NATIVE_DOUBLE, mem3);        \
  } while (0)

void EOS_Tabulated::readtable_leptonic(const char *leptonic_table_name) 
{
  using namespace Margherita_constants;

  constexpr size_t NTABLES_ELE = EOS_Tabulated::EELE::NUM_VARS_ELE;

#ifndef STANDALONE
  CCTK_VInfo(CCTK_THORNSTRING, "*******************************");
  CCTK_VInfo(CCTK_THORNSTRING, "Reading Electronic EoS table file:");
  CCTK_VInfo(CCTK_THORNSTRING, "%s", leptonic_table_name);
  CCTK_VInfo(CCTK_THORNSTRING, "*******************************");
#endif

  hid_t file;
  if (!file_is_readable(leptonic_table_name)) {
#ifndef STANDALONE
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read electronic table %s \n", leptonic_table_name);
#else
    std::cout << "Cannot open table" << std::endl;
#endif
  }

  HDF5_ERROR(file = H5Fopen(leptonic_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

  // nyle is the total number of actual electron fraction in the table, since the original one is Yp (charge fraction)
  int nrho, ntemp, nyle;
  
  // Easier reader
  // Read size of tables
  READ_EOS_HDF5("nrho", &nrho, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("ntemp", &ntemp, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("nyle", &nyle, H5T_NATIVE_INT, H5S_ALL);

  // Allocate memory for tables
  double *logrho = new double[nrho];
  double *logtemp = new double[ntemp];
  double *yle = new double[nyle];


  // Read additional tables and variables
  READ_EOS_HDF5("logrho_table", logrho, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("logtemp_table", logtemp, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("yle_table", yle, H5T_NATIVE_DOUBLE, H5S_ALL);

  READ_EOS_HDF5("eos_ylemin", &eos_ylemin, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("eos_ylemax", &eos_ylemax, H5T_NATIVE_DOUBLE, H5S_ALL);

  std::cout << "\n nrho in the leptons table =\n"<< nrho;
  std::cout << "\n ntemp in the leptons table =\n"<< ntemp;
  std::cout << "\n nyle in the leptons table =\n"<< nyle;

  std::cout << "\n eos_ylemax\n" <<EOS_Tabulated::eos_ylemax;
  std::cout << "\n eos_ylemin\n" <<EOS_Tabulated::eos_ylemin;

 // Read the electronic table

  // Allocate memory of both
  double *electronic_eos_tables = new double[nrho * ntemp * nyle * NTABLES_ELE];

  // Electronic eos tables are 4-D arrays tho. but put all variables to 1D table

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims_ele[3] = {NTABLES_ELE, (hsize_t)nrho * ntemp * nyle * EOS_Tabulated::EELE::NUM_VARS_ELE};
  hsize_t var3_ele[3] = {1, (hsize_t)nrho * ntemp * nyle * EOS_Tabulated::EELE::NUM_VARS_ELE};
  hid_t mem3_ele = H5Screate_simple(3, table_dims_ele, NULL);
  READ_EOS_HDF5("electronic_eos_tables", electronic_eos_tables, H5T_NATIVE_DOUBLE, H5S_ALL);

  HDF5_ERROR(H5Fclose(file));
  HDF5_ERROR(H5Sclose(mem3_ele));
  // End reader

  auto alltables_ele =
      std::unique_ptr<double[]>(new double[nrho * ntemp * nyle * NTABLES_ELE]);

  for (int iv = EELE::MUELE; iv < EELE::NUM_VARS_ELE; iv++)
    for (int k = 0; k < nyle; k++)
      for (int j = 0; j < ntemp; j++)
        for (int i = 0; i < nrho; i++) {
          int indold = i + nrho * (j + ntemp * (k + nyle * iv));
          int indnew = iv + NTABLES_ELE * (i + nrho * (j + ntemp * k));
          alltables_ele[indnew] = electronic_eos_tables[indold];
        }

  delete[] electronic_eos_tables;

  double prs_min_local;
  double prs_max_local;
  double eps_min_local;
  double eps_max_local;

  eos_muelemin    = 1.0e90;
  eos_prs_ele_min = 1.0e90;
  eos_eps_ele_min = 1.0e90;
  // reader checkers
  for (int i = 0; i < nrho * ntemp * nyle; i++) {
    { // check ele
      int idx = EOS_Tabulated::EELE::MUELE + NTABLES_ELE * i;
      eos_muelemax = std::max(alltables_ele[idx], eos_muelemax);
      eos_muelemin = std::min(alltables_ele[idx], eos_muelemin);

      idx = EOS_Tabulated::EELE::PRESS_E_MINUS + NTABLES_ELE * i;
      int jdx = EOS_Tabulated::EELE::PRESS_E_PLUS + NTABLES_ELE * i;
      prs_max_local = alltables_ele[idx] + alltables_ele[jdx];
      prs_min_local = alltables_ele[idx] + alltables_ele[jdx];
      eos_prs_ele_max = std::max(prs_max_local, eos_prs_ele_max);
      eos_prs_ele_min = std::min(prs_min_local, eos_prs_ele_min);

      idx = EOS_Tabulated::EELE::EPS_E_MINUS + NTABLES_ELE * i;
      jdx = EOS_Tabulated::EELE::EPS_E_PLUS + NTABLES_ELE * i;
      eps_max_local = alltables_ele[idx] + alltables_ele[jdx];
      eps_min_local = alltables_ele[idx] + alltables_ele[jdx];
      eos_eps_ele_max = std::max(eps_max_local, eos_eps_ele_max);
      eos_eps_ele_min = std::min(eps_min_local, eos_eps_ele_min);
    }

  }
  std::cout << "\n eos_muelemax\n" << std::scientific <<eos_muelemax;
  std::cout << "\n eos_muelemin\n" << std::scientific <<eos_muelemin;
  std::cout << "\n eos_prs_ele_max\n" << std::scientific <<eos_prs_ele_max;
  std::cout << "\n eos_prs_ele_min\n" << std::scientific <<eos_prs_ele_min;
  std::cout << "\n eos_eps_ele_max\n" << std::scientific <<eos_eps_ele_max;
  std::cout << "\n eos_eps_ele_min\n" << std::scientific <<eos_eps_ele_min;

  auto logrho_ptr1 = std::unique_ptr<double[]>(new double[nrho]);
  auto logtemp_ptr1 = std::unique_ptr<double[]>(new double[ntemp]);
  auto yle_ptr1 = std::unique_ptr<double[]>(new double[nyle]);

  for (int i = 0; i < nrho; ++i) logrho_ptr1[i] = logrho[i];
  for (int i = 0; i < ntemp; ++i) logtemp_ptr1[i] = logtemp[i];
  for (int i = 0; i < nyle; ++i) yle_ptr1[i] = yle[i];

  auto num_points_ele =
      std::array<size_t, 3>{size_t(nrho), size_t(ntemp), size_t(nyle)};


  // Storing alltables to global alltables
  EOS_Tabulated::alltables_ele = linear_interp_uniform_ND_t<double, 3, EELE::NUM_VARS_ELE>(
      std::move(alltables_ele), std::move(num_points_ele), std::move(logrho_ptr1),
      std::move(logtemp_ptr1), std::move(yle_ptr1));

  free(logrho);
  free(logtemp);
  free(yle);
};

#endif
