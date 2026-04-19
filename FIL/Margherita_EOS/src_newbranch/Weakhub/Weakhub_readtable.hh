//
//  Copyright (C) 2021, Harry Ho-Yin Ng
//  Based on routines by Elias Roland Most
//                       Ludwig Jens Papenfort
//

// logrho_IVtable, logT_IVtable, logymu_IVtable are using natural log
// rho in code unit, T in MeV, kappa in cm^{-1}
// kappa_a_st is kappa_a with stimulated absorption factor
// ***** Weakhub microphysics has 3 or 6 species of neutrino 
// However, FIL_M1 uses 3 or 5 species of neutrino to speed up the code
#include "Weakhub.hh"
#define H5_USE_16_API 1
#include <hdf5.h>

#ifndef M1_WEAKHUB_READTABLE_HH
#define M1_WEAKHUB_READTABLE_HH

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                          \
  do {                                                               \
    int _error_code = fn_call;                                       \
    if (_error_code < 0) {                                           \
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,              \
                  "HDF5 call '%s' returned error code %d", #fn_call, \
                  _error_code);                                      \
    }                                                                \
  } while (0)

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

#ifndef EOS_TABULATED_READTABLE_COMPOSE_HH
static inline int file_is_readable(const char *filename) {
  FILE *fp = NULL;
  fp = fopen(filename, "r");
  if (fp != NULL) {
    fclose(fp);
    return 1;
  }
  return 0;
}
#endif

void M1_Weakhub::readtable_Weakhub(const char *Weakhub_table_name) {
  //using namespace Margherita_constants;

//#ifndef STANDALONE
  CCTK_VInfo(CCTK_THORNSTRING, "*******************************");
  CCTK_VInfo(CCTK_THORNSTRING, "Reading Weakhub table file:");
  CCTK_VInfo(CCTK_THORNSTRING, "%s", Weakhub_table_name);
  CCTK_VInfo(CCTK_THORNSTRING, "*******************************");
//#endif

  hid_t file;
  if (!file_is_readable(Weakhub_table_name)) {
//#ifndef STANDALONE
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read Weakhub_table_name %s \n", Weakhub_table_name);
//#else
//    std::cout << "Cannot open table" << std::endl;
//#endif
  }

  HDF5_ERROR(file = H5Fopen(Weakhub_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

  READ_EOS_HDF5("n_species", &n_spec, H5T_NATIVE_INT, H5S_ALL);

  int Weakhub_table_type = 0;

  READ_EOS_HDF5("weakhub_table_type", &Weakhub_table_type, H5T_NATIVE_INT, H5S_ALL);

  if (Weakhub_table_type != 100) {
        std::cout << "Weakhub_table is not a grey table!\n" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 911);
  }

  // Easier reader
  // Read size of tables
  READ_EOS_HDF5("IVrho", &IVrho, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("IVtemp", &IVtemp, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("IVye", &IVye, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("IVymu", &IVymu, H5T_NATIVE_INT, H5S_ALL);

  std::cout << "\n IVrho in the Weakhub table =\n"<< IVrho;
  std::cout << "\n IVtemp in the Weakhub table =\n"<< IVtemp;
  std::cout << "\n IVye in the Weakhub table =\n"<< IVye;
  std::cout << "\n IVymu in the Weakhub table =\n"<< IVymu;

  // Allocate memory for tables
  double *logrho_IVtable  = new double[IVrho];
  double *logtemp_IVtable = new double[IVtemp];
  double *ye_IVtable      = new double[IVye];
  double *logymu_IVtable  = new double[IVymu]; // Ymu is in log scale just like rho, temp


  // Read additional tables and variables
  READ_EOS_HDF5("logrho_IVtable", logrho_IVtable, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("logtemp_IVtable", logtemp_IVtable, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("ye_IVtable", ye_IVtable, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("logymu_IVtable", logymu_IVtable, H5T_NATIVE_DOUBLE, H5S_ALL);

  logrho_max_IV  = logrho_IVtable[IVrho-1];
  logrho_min_IV  = logrho_IVtable[0];
  logtemp_max_IV = logtemp_IVtable[IVtemp-1];
  logtemp_min_IV = logtemp_IVtable[0];
  ye_max_IV      = ye_IVtable[IVye-1];
  ye_min_IV      = ye_IVtable[0];
  logymu_max_IV  = logymu_IVtable[IVymu-1];
  logymu_min_IV  = logymu_IVtable[0];

  std::cout << "\n logrho_max_IV\n" << logrho_max_IV;
  std::cout << "\n logrho_min_IV\n" << logrho_min_IV;
  std::cout << "\n logtemp_max_IV\n" << logtemp_max_IV;
  std::cout << "\n logtemp_min_IV\n" << logtemp_min_IV;
  std::cout << "\n ye_max_IV\n" << ye_max_IV;
  std::cout << "\n ye_min_IV\n" << ye_min_IV;
  std::cout << "\n logymu_max_IV\n" << logymu_max_IV;
  std::cout << "\n logymu_min_IV\n" << logymu_min_IV;

 // Read the opacity tables

  // Allocate memory of opacity tables 
  double *kappa_a_en_grey_tables = new double[IVrho * IVtemp * IVye * IVymu * n_spec];
  double *kappa_a_num_grey_tables = new double[IVrho * IVtemp * IVye * IVymu * n_spec];
  double *kappa_s_grey_tables = new double[IVrho * IVtemp * IVye * IVymu * n_spec];

  // Read 4-D arrays from h5table to 1D arrays

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[4] = {n_spec, (hsize_t)IVrho * IVtemp * IVye * IVymu * n_spec};
  hsize_t var4[4] = {1, (hsize_t)IVrho * IVtemp * IVye * IVymu * n_spec};
  hid_t mem4 = H5Screate_simple(4, table_dims, NULL);
  READ_EOS_HDF5("kappa_a_en_grey_table", kappa_a_en_grey_tables, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("kappa_a_num_grey_table", kappa_a_num_grey_tables, H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("kappa_s_grey_table", kappa_s_grey_tables, H5T_NATIVE_DOUBLE, H5S_ALL);

  HDF5_ERROR(H5Fclose(file));
  HDF5_ERROR(H5Sclose(mem4));
  // End reader

  auto alltables_kappa_a_en =
      std::unique_ptr<double[]>(new double[IVrho * IVtemp * IVye * IVymu * n_spec]);
  auto alltables_kappa_a_num =
      std::unique_ptr<double[]>(new double[IVrho * IVtemp * IVye * IVymu * n_spec]);
  auto alltables_kappa_s =
      std::unique_ptr<double[]>(new double[IVrho * IVtemp * IVye * IVymu * n_spec]);

  for (int iv = NEU_NUM::i_nu_e; iv < n_spec; iv++)
    for (int l = 0; l < IVymu; l++)
      for (int k = 0; k < IVye; k++)
        for (int j = 0; j < IVtemp; j++)
          for (int i = 0; i < IVrho; i++) {
            int indold = i + IVrho * (j + IVtemp * (k + IVye * (l + IVymu * iv)));
            //int indold = i + IVrho * (j + IVtemp * (k + IVye * iv));
            int indnew = iv + n_spec * (i + IVrho * (j + IVtemp * (k + IVye * l)  ) );
            //int indnew = iv + n_spec * (i + IVrho * (j + IVtemp * k));
            alltables_kappa_a_en[indnew]  = kappa_a_en_grey_tables[indold];
            alltables_kappa_a_num[indnew] = kappa_a_num_grey_tables[indold];
            alltables_kappa_s[indnew]     = kappa_s_grey_tables[indold];
          }

  delete[] kappa_a_en_grey_tables;
  delete[] kappa_a_num_grey_tables;
  delete[] kappa_s_grey_tables;

  double kappa_min_local[n_spec]  = {1.0e99};
  double kappa_max_local[n_spec];

  // reader checkers
  for (int i = 0; i < IVrho * IVtemp * IVye * IVymu; i++) {
      int idx;
      for (int i_spec = NEU_NUM::i_nu_e; i_spec < n_spec; i_spec++) {
        idx = i_spec + n_spec * i;
        kappa_max_local[i_spec] = std::max(alltables_kappa_a_en[idx], kappa_max_local[i_spec]);
        kappa_min_local[i_spec] = std::min(alltables_kappa_a_en[idx], kappa_min_local[i_spec]);
      }
  }
  
  for (int i_spec = i_nu_e; i_spec < n_spec; i_spec++) {
    std::cout << "\n kappa_max_local[i_spec]\n" << std::scientific <<kappa_max_local[i_spec];
    std::cout << "\n kappa_min_local[i_spec]\n" << std::scientific <<kappa_min_local[i_spec];
  }
  //std::cout << "Passed checker of kappa_a_en\n";
  //MPI_Abort(MPI_COMM_WORLD, 911);
  
  auto logrho_IVtable_ptr1 = std::unique_ptr<double[]>(new double[IVrho]);
  auto logtemp_IVtable_ptr1 = std::unique_ptr<double[]>(new double[IVtemp]);
  auto ye_IVtable_ptr1 = std::unique_ptr<double[]>(new double[IVye]);

  for (int i = 0; i < IVrho; ++i)  logrho_IVtable_ptr1[i] = logrho_IVtable[i];
  for (int i = 0; i < IVtemp; ++i) logtemp_IVtable_ptr1[i] = logtemp_IVtable[i];
  for (int i = 0; i < IVye; ++i)   ye_IVtable_ptr1[i] = ye_IVtable[i];

  auto logrho_IVtable_ptr2 = std::unique_ptr<double[]>(new double[IVrho]);
  auto logtemp_IVtable_ptr2 = std::unique_ptr<double[]>(new double[IVtemp]);
  auto ye_IVtable_ptr2 = std::unique_ptr<double[]>(new double[IVye]);

  for (int i = 0; i < IVrho; ++i)  logrho_IVtable_ptr2[i] = logrho_IVtable[i];
  for (int i = 0; i < IVtemp; ++i) logtemp_IVtable_ptr2[i] = logtemp_IVtable[i];
  for (int i = 0; i < IVye; ++i)   ye_IVtable_ptr2[i] = ye_IVtable[i];

  auto logrho_IVtable_ptr3 = std::unique_ptr<double[]>(new double[IVrho]);
  auto logtemp_IVtable_ptr3 = std::unique_ptr<double[]>(new double[IVtemp]);
  auto ye_IVtable_ptr3 = std::unique_ptr<double[]>(new double[IVye]);

  for (int i = 0; i < IVrho; ++i)  logrho_IVtable_ptr3[i] = logrho_IVtable[i];
  for (int i = 0; i < IVtemp; ++i) logtemp_IVtable_ptr3[i] = logtemp_IVtable[i];
  for (int i = 0; i < IVye; ++i)   ye_IVtable_ptr3[i] = ye_IVtable[i];


  // Storing alltables to global alltables

  // Warning: everytime after allocate the linear_interp_ND --> your ptr will disappear, you have to make ptr1, 2, 3 for this
  if (n_spec == 6) {
    auto num_points1 =
        std::array<size_t, 4>{size_t(IVrho), size_t(IVtemp), size_t(IVye), size_t(IVymu)};
    auto num_points2 =
        std::array<size_t, 4>{size_t(IVrho), size_t(IVtemp), size_t(IVye), size_t(IVymu)};
    auto num_points3 =
        std::array<size_t, 4>{size_t(IVrho), size_t(IVtemp), size_t(IVye), size_t(IVymu)};
    auto logymu_IVtable_ptr1 = std::unique_ptr<double[]>(new double[IVymu]);
    auto logymu_IVtable_ptr2 = std::unique_ptr<double[]>(new double[IVymu]);
    auto logymu_IVtable_ptr3 = std::unique_ptr<double[]>(new double[IVymu]);
    for (int i = 0; i < IVymu; ++i)  {
         logymu_IVtable_ptr1[i] = logymu_IVtable[i];
         logymu_IVtable_ptr2[i] = logymu_IVtable[i];
         logymu_IVtable_ptr3[i] = logymu_IVtable[i];
    }
         M1_Weakhub::kappa_a_st_en_table_6spec = linear_interp_uniform_ND_t<double, 4, 6>(
             std::move(alltables_kappa_a_en), std::move(num_points1), std::move(logrho_IVtable_ptr1),
             std::move(logtemp_IVtable_ptr1), std::move(ye_IVtable_ptr1), std::move(logymu_IVtable_ptr1));
       
         M1_Weakhub::kappa_a_st_num_table_6spec = linear_interp_uniform_ND_t<double, 4, 6>(
             std::move(alltables_kappa_a_num), std::move(num_points2), std::move(logrho_IVtable_ptr2),
             std::move(logtemp_IVtable_ptr2), std::move(ye_IVtable_ptr2), std::move(logymu_IVtable_ptr2));
       
         M1_Weakhub::kappa_s_table_6spec = linear_interp_uniform_ND_t<double, 4, 6>(
             std::move(alltables_kappa_s), std::move(num_points3), std::move(logrho_IVtable_ptr3),
             std::move(logtemp_IVtable_ptr3), std::move(ye_IVtable_ptr3), std::move(logymu_IVtable_ptr3));
  } else if (n_spec == 3){
    auto num_points1 =
        std::array<size_t, 3>{size_t(IVrho), size_t(IVtemp), size_t(IVye)};
    auto num_points2 =
        std::array<size_t, 3>{size_t(IVrho), size_t(IVtemp), size_t(IVye)};
    auto num_points3 =
        std::array<size_t, 3>{size_t(IVrho), size_t(IVtemp), size_t(IVye)};

         M1_Weakhub::kappa_a_st_en_table_3spec = linear_interp_uniform_ND_t<double, 3, 3>(
             std::move(alltables_kappa_a_en), std::move(num_points1), std::move(logrho_IVtable_ptr1),
             std::move(logtemp_IVtable_ptr1), std::move(ye_IVtable_ptr1));

         M1_Weakhub::kappa_a_st_num_table_3spec = linear_interp_uniform_ND_t<double, 3, 3>(
             std::move(alltables_kappa_a_num), std::move(num_points2), std::move(logrho_IVtable_ptr2),
             std::move(logtemp_IVtable_ptr2), std::move(ye_IVtable_ptr2));

         M1_Weakhub::kappa_s_table_3spec = linear_interp_uniform_ND_t<double, 3, 3>(
             std::move(alltables_kappa_s), std::move(num_points3), std::move(logrho_IVtable_ptr3),
             std::move(logtemp_IVtable_ptr3), std::move(ye_IVtable_ptr3));
  } else {
      std::cout << "n_spec not equal to 6 or 3, Wrong n_species inside the table" << std::endl;
  }

   // code test   passed num_spec = 3 and also 5!
   //int num_spec = 5;
   //double *kappa_a_st_en  = new double[num_spec];
   //double *kappa_a_st_num  = new double[num_spec];
   //double *kappa_s  = new double[num_spec];
   //double *eta  = new double[num_spec];
   //typename M1_Weakhub::error_type error1;
   //double temp_pt = 1.0;
   //double rho_pt  = 8.9850285983761747e-4;
   //double yle_pt = 0.3, ymu_pt = 0.01;

   // //for (int i = 0; i < IVye; ++i) {
   // //  yle_pt = ye_IVtable[i];
   // for (int i = 0; i < IVrho; ++i) {
   //   rho_pt = exp(logrho_IVtable[i]);
   //   //std::cout << "\n rho after beta eqm muonic=\n"<< rho;
   //   //std::cout << "\n yle after beta eqm muonic=\n"<< yle;
   //   //std::cout << "\n ymu after beta eqm muonic=\n"<< ymu;
   //   kappa_ast_kappa_s__temp_rho_yle_ymu(
   //       kappa_a_st_en, kappa_a_st_num, kappa_s, num_spec,  temp_pt,
   //       rho_pt,  yle_pt, ymu_pt, error1);

   //   if (error1[M1_Weakhub::errors::YEIV_TOO_LOW] == true || error1[M1_Weakhub::errors::YEIV_TOO_HIGH] == true 
   //       || error1[M1_Weakhub::errors::YMUIV_TOO_LOW] == true || error1[M1_Weakhub::errors::YMUIV_TOO_HIGH] == true) {
   //      std::cout << error1 << "\n";
   //      std::cout << "\n ye/ymu out of bounds in Weakhub table \n" << std::endl;
   //      MPI_Abort(MPI_COMM_WORLD, 911);
   //   }
   //   for (int j = 0; j < num_spec; ++j) {
   //     std::cout << "\n kappa =\n"<< kappa_a_st_en[j]<<" "<<kappa_a_st_num[j]<<" "<<kappa_s[j]<<" "<<i<<" "<<j<<" \n";
   //   }
   // }
   //MPI_Abort(MPI_COMM_WORLD, 911);
   // end of code test


  free(logrho_IVtable);
  free(logtemp_IVtable);
  free(ye_IVtable);
  free(logymu_IVtable);

};

#endif
