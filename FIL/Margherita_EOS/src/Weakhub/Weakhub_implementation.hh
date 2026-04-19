//
//  Copyright (C) 2021, Harry Ho-Yin Ng
//  Based on routines by Elias Roland Most
//                       Ludwig Jens Papenfort
//
#include <array>
#include <bitset>
#include <functional>
#include <iostream>

#include "../utils/brent.hh"
#include "../utils/spline.hh"
#include "Weakhub.hh"

#ifndef M1_WEAKHUB_IMPL_HH
#define M1_WEAKHUB_IMPL_HH

// Initialize storage for static members
// Harry: If you do not initialize storage --> compile undefine
int M1_Weakhub::n_spec;
int M1_Weakhub::IVrho;
int M1_Weakhub::IVtemp;
int M1_Weakhub::IVye;
int M1_Weakhub::IVymu;
double M1_Weakhub::logrho_max_IV;
double M1_Weakhub::logrho_min_IV;
double M1_Weakhub::logtemp_max_IV;
double M1_Weakhub::logtemp_min_IV;
double M1_Weakhub::logymu_max_IV;
double M1_Weakhub::logymu_min_IV;
double M1_Weakhub::ye_max_IV;
double M1_Weakhub::ye_min_IV;
bool M1_Weakhub::use_Weakhub_eas = false;

linear_interp_uniform_ND_t<double, 4, 6>
    M1_Weakhub::kappa_a_st_en_table_6spec;
linear_interp_uniform_ND_t<double, 4, 6>
    M1_Weakhub::kappa_a_st_num_table_6spec;
linear_interp_uniform_ND_t<double, 4, 6>
    M1_Weakhub::kappa_s_table_6spec;
linear_interp_uniform_ND_t<double, 3, 3>
    M1_Weakhub::kappa_a_st_en_table_3spec;
linear_interp_uniform_ND_t<double, 3, 3>
    M1_Weakhub::kappa_a_st_num_table_3spec;
linear_interp_uniform_ND_t<double, 3, 3>
    M1_Weakhub::kappa_s_table_3spec;

// ###################################################################################
// 	Private functions
// ###################################################################################
template <bool check_temp>
M1_Weakhub::error_type M1_Weakhub::checkbounds(double &xrho,
                                               double &xtemp,
                                               double &xye,
                                               double &xymu, int &num_spec) {
  // initialize error codes to zero
  error_type error_codes;

  const auto lrho = log(xrho);
  const auto ltemp = log(xtemp);

  if (lrho > logrho_max_IV) {
    // fix
      error_codes[errors::RHOIV_TOO_HIGH] = true;
      xrho = exp(logrho_max_IV);
  } else if (lrho < logrho_min_IV) {
      // fix
      error_codes[errors::RHOIV_TOO_LOW] = true;
      xrho = exp(logrho_min_IV);
  }
  if (xye > ye_max_IV) {
      error_codes[errors::YEIV_TOO_HIGH] = true;
      //std::cout << "\n yle > ye_max_IV in Weakhub table \n" << std::endl;
      //MPI_Abort(MPI_COMM_WORLD, 911);
  } else if (xye < ye_min_IV) {
      error_codes[errors::YEIV_TOO_LOW] = true;
      //std::cout << "\n yle < ye_min_IV in Weakhub table \n" << std::endl;
      //MPI_Abort(MPI_COMM_WORLD, 911);
  }

  if (num_spec == 5 || num_spec == 6) {
    if (log(xymu) > logymu_max_IV) {
      error_codes[errors::YMUIV_TOO_HIGH] = true;
      //std::cout << "\n log(ymu) > logymu_max_IV in Weakhub table \n" << std::endl;
      //MPI_Abort(MPI_COMM_WORLD, 911);
    } else if (log(xymu) < logymu_min_IV) {
      error_codes[errors::YMUIV_TOO_LOW] = true;
      //std::cout << "\n log(ymu) < logymu_min_IV in Weakhub table \n" << std::endl;
      //MPI_Abort(MPI_COMM_WORLD, 911);
    }
  }

  if (ltemp > logtemp_max_IV) {
      // fix
      error_codes[errors::TEMPIV_TOO_HIGH] = true;
      xtemp = exp(logtemp_max_IV);
  } else if (ltemp < logtemp_min_IV) {
      // fix
      error_codes[errors::TEMPIV_TOO_LOW] = true;
      xtemp = exp(logtemp_min_IV);
  }
  return error_codes;
}
// ###################################################################################
//      Public functions
// ###################################################################################
void M1_Weakhub::kappa_ast_kappa_s__temp_rho_yle_ymu(
      double kappa_a_st_en[], double kappa_a_st_num[], double kappa_s[], int num_spec, double &temp,
      double &rho, double &yle, double &ymu, error_type &error){

  error = checkbounds<true>(rho, temp, yle, ymu, num_spec);
  
  if (error[errors::TEMPIV_TOO_LOW] == true || error[errors::RHOIV_TOO_LOW] == true) {
     for (int i = 0; i < num_spec; i++) {
        kappa_a_st_en[i] = 0.0; 
        kappa_a_st_num[i] = 0.0; 
        kappa_s[i] = 0.0; 
        return;
     }
  }

  const double lrho = log(rho);
  const double ltemp = log(temp);

  if (num_spec == 3) {
     auto const kappa_a_st_en_local = 
        kappa_a_st_en_table_3spec.interpolate<NEU_NUM::i_nu_e, NEU_NUM::i_nu_e_bar, NEU_NUM::i_nu_mu>(lrho, ltemp, yle);
     auto const kappa_a_st_num_local = 
        kappa_a_st_num_table_3spec.interpolate<NEU_NUM::i_nu_e, NEU_NUM::i_nu_e_bar, NEU_NUM::i_nu_mu>(lrho, ltemp, yle);
     auto const kappa_s_local = 
        kappa_s_table_3spec.interpolate<NEU_NUM::i_nu_e, NEU_NUM::i_nu_e_bar, NEU_NUM::i_nu_mu>(lrho, ltemp, yle);
     for (int i = 0; i < num_spec; i++) {
       kappa_a_st_en[i]  = kappa_a_st_en_local[i];
       kappa_a_st_num[i] = kappa_a_st_num_local[i];
       kappa_s[i]        = kappa_s_local[i];
     }
  }
  else if (num_spec == 5) {
    // this option treats kappa[i_nu_tau] = kappa[i_nu_tau_bar]
    const double lymu = log(ymu);

    auto const kappa_a_st_en_local =
       kappa_a_st_en_table_6spec.interpolate
         <NEU_NUM::i_nu_e, NEU_NUM::i_nu_e_bar, NEU_NUM::i_nu_mu, NEU_NUM::i_nu_mu_bar, NEU_NUM::i_nu_tau>(lrho, ltemp, yle, lymu);
    auto const kappa_a_st_num_local =
       kappa_a_st_num_table_6spec.interpolate
         <NEU_NUM::i_nu_e, NEU_NUM::i_nu_e_bar, NEU_NUM::i_nu_mu, NEU_NUM::i_nu_mu_bar, NEU_NUM::i_nu_tau>(lrho, ltemp, yle, lymu);
    auto const kappa_s_local =
       kappa_s_table_6spec.interpolate
         <NEU_NUM::i_nu_e, NEU_NUM::i_nu_e_bar, NEU_NUM::i_nu_mu, NEU_NUM::i_nu_mu_bar, NEU_NUM::i_nu_tau>(lrho, ltemp, yle, lymu);
     for (int i = 0; i < num_spec; i++) {
       kappa_a_st_en[i]  = kappa_a_st_en_local[i];
       kappa_a_st_num[i] = kappa_a_st_num_local[i];
       kappa_s[i]        = kappa_s_local[i];
     }
  }
  else {
     std::cout << "\n Unknown num_spec inside kappa_ast_kappa_s \n" << std::endl;
     //MPI_Abort(MPI_COMM_WORLD, 911);
     abort();
  }
  // why I cannot put there?
     //for (int i = 0; i < num_spec; i++) {
     //  kappa_a_st_en[i]  = kappa_a_st_en_local[i];
     //  kappa_a_st_num[i] = kappa_a_st_num_local[i];
     //  kappa_s[i]        = kappa_s_local[i];
     //}


  return;

}

#include "Weakhub_readtable.hh"

#endif
