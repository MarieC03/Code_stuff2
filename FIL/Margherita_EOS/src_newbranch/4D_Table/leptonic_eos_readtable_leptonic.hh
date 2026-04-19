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
// logrho, logT, logymu are using natural log


#include "leptonic_eos.hh"
#include <fstream>
#include <iostream>

#define H5_USE_16_API 1
#include <hdf5.h>

#ifndef EOS_LEPTONIC_READTABLE_HH
#define EOS_LEPTONIC_READTABLE_HH

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

void EOS_Leptonic::readtable_leptonic(const char *leptonic_eos_table_name) 
{
  using namespace Margherita_constants;

  constexpr size_t NTABLES_ELE = EOS_Leptonic::EELE::NUM_VARS_ELE;
  constexpr size_t NTABLES_MUON = EOS_Leptonic::EMUON::NUM_VARS_MUON;

#ifndef STANDALONE
  CCTK_VInfo(CCTK_THORNSTRING, "*******************************");
  CCTK_VInfo(CCTK_THORNSTRING, "Reading Leptonic EoS table file:");
  CCTK_VInfo(CCTK_THORNSTRING, "%s", leptonic_eos_table_name);
  CCTK_VInfo(CCTK_THORNSTRING, "*******************************");
#endif

  hid_t file;
  if (!file_is_readable(leptonic_eos_table_name)) {
#ifndef STANDALONE
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read leptonic_eos_table_name %s \n", leptonic_eos_table_name);
#else
    std::cout << "Cannot open table" << std::endl;
#endif
  }

  HDF5_ERROR(file = H5Fopen(leptonic_eos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

  // nyle is the total number of actual electron fraction in the table, since the original one is Yp (charge fraction)
  int nrho, ntemp, nymu, nyle;
  
  // Easier reader
  // Read size of tables
  READ_EOS_HDF5("nrho", &nrho, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("ntemp", &ntemp, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("nyle", &nyle, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("nymu", &nymu, H5T_NATIVE_INT, H5S_ALL);

  // Allocate memory for tables
  double *logrho = new double[nrho];
  double *logtemp = new double[ntemp];
  double *yle = new double[nyle];
  double *logymu = new double[nymu]; // Ymu is in log scale just like rho, temp


  // Read additional tables and variables
  READ_EOS_HDF5("logrho_table", logrho, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("logtemp_table", logtemp, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("yle_table", yle, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("ymu_table", logymu, H5T_NATIVE_DOUBLE, H5S_ALL);

  READ_EOS_HDF5("eos_ymumin", &eos_ymumin, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("eos_ymumax", &eos_ymumax, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("eos_ylemin", &eos_ylemin, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("eos_ylemax", &eos_ylemax, H5T_NATIVE_DOUBLE, H5S_ALL);

  if (eos_ylemax != eos_yemax) {
   std::cout << "\n EOS_Leptonic::eos_ylemax != EOS_Leptonic::eos_yemax ="<< eos_ylemax <<" "<<EOS_Leptonic::eos_yemax<<"\n";
   MPI_Abort(MPI_COMM_WORLD, 911);
  }
  if (eos_ylemin != eos_yemin) {
   std::cout << "\n EOS_Leptonic::eos_ylemin != EOS_Leptonic::eos_yemin ="<< eos_ylemin <<" "<<EOS_Leptonic::eos_yemin<<"\n";
   MPI_Abort(MPI_COMM_WORLD, 911);
  }
  
  std::cout << "\n nrho in the leptons table =\n"<< nrho;
  std::cout << "\n ntemp in the leptons table =\n"<< ntemp;
  std::cout << "\n nyle in the leptons table =\n"<< nyle;
  std::cout << "\n nymu in the leptons table =\n"<< nymu;


  std::cout << "\n eos_ylemax\n" <<EOS_Leptonic::eos_ylemax;
  std::cout << "\n eos_ylemin\n" <<EOS_Leptonic::eos_ylemin;
  std::cout << "\n eos_ymumax\n" <<EOS_Leptonic::eos_ymumax;
  std::cout << "\n eos_ymumin\n" <<EOS_Leptonic::eos_ymumin;

 // Read the electronic and muonic tables

  // Allocate memory of both
  double *electronic_eos_tables = new double[nrho * ntemp * nyle * NTABLES_ELE];
  double *muonic_eos_tables = new double[nrho * ntemp * nymu * NTABLES_MUON];

  // muonic/electronic eos tables are 4-D arrays tho. but put all variables to 1D table

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims_ele[3] = {NTABLES_ELE, (hsize_t)nrho * ntemp * nyle * EOS_Leptonic::EELE::NUM_VARS_ELE};
  hsize_t var3_ele[3] = {1, (hsize_t)nrho * ntemp * nyle * EOS_Leptonic::EELE::NUM_VARS_ELE};
  hid_t mem3_ele = H5Screate_simple(3, table_dims_ele, NULL);
  READ_EOS_HDF5("electronic_eos_tables", electronic_eos_tables, H5T_NATIVE_DOUBLE, H5S_ALL);

  hsize_t table_dims_muon[3] = {NTABLES_MUON, (hsize_t)nrho * ntemp * nymu * EOS_Leptonic::EMUON::NUM_VARS_MUON};
  hsize_t var3_muon[3] = {1, (hsize_t)nrho * ntemp * nymu * EOS_Leptonic::EMUON::NUM_VARS_MUON};
  hid_t mem3_muon = H5Screate_simple(3, table_dims_muon, NULL);
  READ_EOS_HDF5("muonic_eos_tables", muonic_eos_tables, H5T_NATIVE_DOUBLE, H5S_ALL);

  HDF5_ERROR(H5Fclose(file));
  HDF5_ERROR(H5Sclose(mem3_ele));
  HDF5_ERROR(H5Sclose(mem3_muon));
  // End reader

  constexpr size_t MUL_C = 1;
  constexpr size_t YL_MINUS_C = 2;
  constexpr size_t YL_PLUS_C = 3;
  constexpr size_t PRESS_L_MINUS_C = 4;
  constexpr size_t PRESS_L_PLUS_C = 5;
  constexpr size_t EPS_L_MINUS = 6;
  constexpr size_t EPS_L_PLUS = 7;
  constexpr size_t S_L_MINUS = 8;
  constexpr size_t S_L_PLUS = 9;

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

  auto alltables_muon =
      std::unique_ptr<double[]>(new double[nrho * ntemp * nymu * NTABLES_MUON]);

  for (int iv = EMUON::MUMU; iv < EMUON::NUM_VARS_MUON; iv++)
    for (int k = 0; k < nymu; k++)
      for (int j = 0; j < ntemp; j++)
        for (int i = 0; i < nrho; i++) {
          int indold = i + nrho * (j + ntemp * (k + nymu * iv));
          int indnew = iv + NTABLES_MUON * (i + nrho * (j + ntemp * k));
          alltables_muon[indnew] = muonic_eos_tables[indold];
        }

  delete[] electronic_eos_tables;
  delete[] muonic_eos_tables;


  double prs_min_local;
  double prs_max_local;
  double eps_min_local;
  double eps_max_local;

  eos_mumumin     = 1.0e90;
  eos_muelemin    = 1.0e90;
  eos_prs_mu_min  = 1.0e90;
  eos_prs_ele_min = 1.0e90;
  eos_eps_mu_min  = 1.0e90;
  eos_eps_ele_min = 1.0e90;
  // reader checkers
  for (int i = 0; i < nrho * ntemp * nymu; i++) {
    {  // check muon
      int idx = EOS_Leptonic::EMUON::MUMU + NTABLES_MUON * i;
      //std::cout << alltables_muon[idx] << "\n";
      eos_mumumax = std::max(alltables_muon[idx], eos_mumumax);
      eos_mumumin = std::min(alltables_muon[idx], eos_mumumin);
      idx = EOS_Leptonic::EMUON::PRESS_MU_MINUS + NTABLES_MUON * i;
      int jdx = EOS_Leptonic::EMUON::PRESS_MU_PLUS + NTABLES_MUON * i;
      prs_max_local = alltables_muon[idx] + alltables_muon[jdx];
      prs_min_local = alltables_muon[idx] + alltables_muon[jdx];
      eos_prs_mu_max = std::max(prs_max_local, eos_prs_mu_max);
      eos_prs_mu_min = std::min(prs_min_local, eos_prs_mu_min);
      idx = EOS_Leptonic::EMUON::EPS_MU_MINUS + NTABLES_MUON * i;
      jdx = EOS_Leptonic::EMUON::EPS_MU_PLUS + NTABLES_MUON * i;
      eps_max_local = alltables_muon[idx] + alltables_muon[jdx];
      eps_min_local = alltables_muon[idx] + alltables_muon[jdx];
      eos_eps_mu_max = std::max(eps_max_local, eos_eps_mu_max);
      eos_eps_mu_min = std::min(eps_min_local, eos_eps_mu_min);
    }
  }
  for (int i = 0; i < nrho * ntemp * nymu; i++) {
    { // check ele
      int idx = EOS_Leptonic::EELE::MUELE + NTABLES_ELE * i;
      eos_muelemax = std::max(alltables_ele[idx], eos_muelemax);
      eos_muelemin = std::min(alltables_ele[idx], eos_muelemin);

      idx = EOS_Leptonic::EELE::PRESS_E_MINUS + NTABLES_ELE * i;
      int jdx = EOS_Leptonic::EELE::PRESS_E_PLUS + NTABLES_ELE * i;
      prs_max_local = alltables_ele[idx] + alltables_ele[jdx];
      prs_min_local = alltables_ele[idx] + alltables_ele[jdx];
      eos_prs_ele_max = std::max(prs_max_local, eos_prs_ele_max);
      eos_prs_ele_min = std::min(prs_min_local, eos_prs_ele_min);

      idx = EOS_Leptonic::EELE::EPS_E_MINUS + NTABLES_ELE * i;
      jdx = EOS_Leptonic::EELE::EPS_E_PLUS + NTABLES_ELE * i;
      eps_max_local = alltables_ele[idx] + alltables_ele[jdx];
      eps_min_local = alltables_ele[idx] + alltables_ele[jdx];
      eos_eps_ele_max = std::max(eps_max_local, eos_eps_ele_max);
      eos_eps_ele_min = std::min(eps_min_local, eos_eps_ele_min);
    }

  }
  std::cout << "\n eos_mumumax\n" << std::scientific <<eos_mumumax;
  std::cout << "\n eos_mumumin\n" << std::scientific <<eos_mumumin;
  std::cout << "\n eos_muelemax\n" << std::scientific <<eos_muelemax;
  std::cout << "\n eos_muelemin\n" << std::scientific <<eos_muelemin;
  std::cout << "\n eos_prs_mu_max\n" << std::scientific <<eos_prs_mu_max;
  std::cout << "\n eos_prs_mu_min\n" << std::scientific <<eos_prs_mu_min;
  std::cout << "\n eos_eps_mu_max\n" << std::scientific <<eos_eps_mu_max;
  std::cout << "\n eos_eps_mu_min\n" << std::scientific <<eos_eps_mu_min;
  std::cout << "\n eos_prs_ele_max\n" << std::scientific <<eos_prs_ele_max;
  std::cout << "\n eos_prs_ele_min\n" << std::scientific <<eos_prs_ele_min;
  std::cout << "\n eos_eps_ele_max\n" << std::scientific <<eos_eps_ele_max;
  std::cout << "\n eos_eps_ele_min\n" << std::scientific <<eos_eps_ele_min;
  //std::cout << "Passed checker of mu mu\n";
  //MPI_Abort(MPI_COMM_WORLD, 911);
  
  // eos_hmin for Kc2p is  h_min = 1 + (eps+eps_l) + (P+P_l)/rho
  //c2p_h_mu_min    = 1.0e90;
  //c2p_h_ele_min   = 1.0e90;
  //double h_min_local;
  //for (int i = 0; i < nrho; i++) {
  //   h_min_local = 1.0e0 + eos_eps_mu_min + eos_prs_mu_min/ exp(logrho[i]);
  //   c2p_h_mu_min = std::min(c2p_h_mu_min, h_min_local);
  //   h_min_local = 1.0e0 + eos_eps_ele_min + eos_prs_ele_min/ exp(logrho[i]);
  //   c2p_h_ele_min = std::min(c2p_h_ele_min, h_min_local);
  //}
  //std::cout << "\n c2p_h_mu_min\n" << std::scientific  <<c2p_h_mu_min;
  //std::cout << "\n c2p_h_ele_min\n" << std::scientific <<c2p_h_ele_min;

  auto logrho_ptr1 = std::unique_ptr<double[]>(new double[nrho]);
  auto logtemp_ptr1 = std::unique_ptr<double[]>(new double[ntemp]);
  auto yle_ptr1 = std::unique_ptr<double[]>(new double[nyle]);

  for (int i = 0; i < nrho; ++i) logrho_ptr1[i] = logrho[i];
  for (int i = 0; i < ntemp; ++i) logtemp_ptr1[i] = logtemp[i];
  for (int i = 0; i < nyle; ++i) yle_ptr1[i] = yle[i];

  auto num_points_ele =
      std::array<size_t, 3>{size_t(nrho), size_t(ntemp), size_t(nyle)};

  auto logrho_ptr2 = std::unique_ptr<double[]>(new double[nrho]);
  auto logtemp_ptr2 = std::unique_ptr<double[]>(new double[ntemp]);
  auto logymu_ptr2 = std::unique_ptr<double[]>(new double[nymu]);

  for (int i = 0; i < nrho; ++i) logrho_ptr2[i] = logrho[i];
  for (int i = 0; i < ntemp; ++i) logtemp_ptr2[i] = logtemp[i];
  for (int i = 0; i < nymu; ++i) logymu_ptr2[i] = logymu[i];

  auto num_points_muon =
      std::array<size_t, 3>{size_t(nrho), size_t(ntemp), size_t(nymu)};


  // Storing alltables to global alltables
  EOS_Leptonic::alltables_ele = linear_interp_uniform_ND_t<double, 3, EELE::NUM_VARS_ELE>(
      std::move(alltables_ele), std::move(num_points_ele), std::move(logrho_ptr1),
      std::move(logtemp_ptr1), std::move(yle_ptr1));

  EOS_Leptonic::alltables_muon = linear_interp_uniform_ND_t<double, 3, EMUON::NUM_VARS_MUON>(
      std::move(alltables_muon), std::move(num_points_muon), std::move(logrho_ptr2),
      std::move(logtemp_ptr2), std::move(logymu_ptr2));

  // updating eos_h_min and eos_h_max:
  // since h_min(baryonic) < h_min(baryonic+leptonic), it is also okay for Kc2p

  

  // start of code test //
  ////code test Passed!  when T < T_thr, the muonic beta eqm condition has bad convergence for rho > 1e-3
  //// T_thr = 0.7 MeV for SFHo
  //// T_thr = 1.5 MeV for DD2, but DD2 tends to have lower density, those unphysical regions not relevant
  //typename EOS_Leptonic::error_type error1;
  //typename EOS_Leptonic::error_type error2;
  //// passed very high density region, but low density region, there is no muon -->
  //double temp_pt = 1;
  //double rho_pt  = 1e-14;
  ////double rho_pt  = 8.9850285983761747e-4;
  //double yle_pt, ymu_pt, eps_pt, yp_pt;
  //double press_pt;
  //double delta, temp_new_pt;


  //std::ofstream myfile ("press_yle_ymu_muonic_1mev_SFHo_merged_branch.log");
  //std::ofstream myfile2 ("press_yp_1mev_SFHo_merged_branch.log");

  //  for (int i = 0; i < nrho; ++i) {
  //    rho_pt = exp(logrho[i]);
  //    //std::cout << "\n rho yle ymuafter beta eqm muonic=\n"<< rho_pt<<" "<<yle_pt<<" "<<ymu_pt<<" \n";
  //    //std::cout << "\n rho after beta eqm muonic=\n"<< rho;
  //    //std::cout << "\n yle after beta eqm muonic=\n"<< yle;
  //    //std::cout << "\n ymu after beta eqm muonic=\n"<< ymu;

  //    press_pt = EOS_Leptonic::press_eps_yle_ymu__beta_eq__rho_temp(eps_pt, yle_pt,
  //                                                                         ymu_pt, rho_pt, temp_pt, error1);
  //    //myfile << rho_pt <<" "<< eps_pt<<" "<<yle_pt<<" "<<ymu_pt << "\n";
  //    myfile << rho_pt <<" "<< press_pt<<" "<<yle_pt<<" "<<ymu_pt << "\n";

  //    //std::cout << "\n press eps yle ymu after beta eqm muonic=\n"<<
  //    //              press_pt<<" "<< eps_pt<<" "<<yle_pt<<" "<<ymu_pt<<" \n";

  //    press_pt = EOS_Leptonic::press_eps_ye__beta_eq__rho_temp(eps_pt, yp_pt,
  //                                                              rho_pt, temp_pt, error2);
  //    //std::cout << "\n press eps yp after beta eqm new=\n"<<
  //    //              rho_pt<<" "<< press_pt<<" "<<yp_pt<<" \n";

  //    //myfile2 << rho_pt <<" "<< eps_pt<<" "<<yp_pt<<" "<<0.0 << "\n";
  //    myfile2 << rho_pt <<" "<< press_pt<<" "<<yp_pt<<" "<<0.0 << "\n";
  //
  //    // code test passed, but for T < 1 MeV, for high Ymu, the delta of inversion of eps --> T could be ~ 10%
  //    // Physically, at these pts are not likely to have so much muons!!
  //    // even tolerance of rootfinding does not change
  //    /*
  //    for (int j = 0; j < ntemp; ++j) {
  //      for (int k = 0; k < nyle; ++k) {
  //        yle_pt = yle[k];
  //        temp_pt = exp(logtemp[j]);
  //        eps_pt = EOS_Leptonic::eps__temp_rho_yle_ymu(temp_pt, rho_pt, yle_pt, ymu_pt,
  //                                                     error1);
  //        press_pt = EOS_Leptonic::press_temp__eps_rho_yle_ymu(temp_new_pt, eps_pt,
  //                                                       rho_pt, yle_pt, ymu_pt,
  //                                                       error1);
  //        delta = (temp_new_pt - temp_pt)/temp_pt;
  //        if (std::abs(delta) > 1e-10) {
  //          std::cout << "\n delta T=\n"<<
  //                        rho_pt<<" "<< temp_pt<<" "<<yle_pt<<" "<<delta<<" \n";
  //         }
  //      }
  //    }
  //    */
  //    
  //  }
  //myfile.close();
  //myfile2.close();
  //MPI_Abort(MPI_COMM_WORLD, 911);
  // -------end of code test //

  if(EOS_Leptonic::generate_cold_table_local){
   EOS_Leptonic::generate_cold_table(logrho, nrho);
  } 
  free(logrho);
  free(logtemp);
  free(yle);
  free(logymu);
};

#endif
