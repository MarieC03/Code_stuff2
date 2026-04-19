//
//  Copyright (C) 2021, Harry Ho-Yin Ng
//  Based on routines by Elias Roland Most
//                       Ludwig Jens Papenfort
//

//#include "cctk.h"
#include <array>
#include <bitset>
#include <iostream>
#include "../margherita.hh"
#include "../utils/brent.hh"
#include "../utils/linear_interp_ND.hh"

#ifndef M1_WEAKHUB_HH
#define M1_WEAKHUB_HH

class M1_Weakhub {
 public:
  // error definitions
  enum errors {
    NO_ERRORS = 0,
    RHOIV_TOO_HIGH,
    RHOIV_TOO_LOW,
    TEMPIV_TOO_HIGH,
    TEMPIV_TOO_LOW,
    YEIV_TOO_HIGH,
    YEIV_TOO_LOW,
    YMUIV_TOO_HIGH,
    YMUIV_TOO_LOW,
    num_errors
  };
  // enum checked_vars { RHO, TEMP, YE, YMU };
  // return type of errors
  typedef std::bitset<errors::num_errors> error_type;

 // neutrino number
  enum NEU_NUM {
    i_nu_e = 0,
    i_nu_e_bar,
    i_nu_mu,
    i_nu_mu_bar,
    i_nu_tau,
    i_nu_tau_bar
  };

 private:
  template <bool check_temp>
  static inline error_type checkbounds(double &xrho, double &xtemp,
                                       double &xye, double &xymu, int &num_spec);

 public:

  static void readtable_Weakhub(const char *Weakhub_table_name);
  static void test_Weakhub();

  // For three species case
  //static int    i_nu_x, i_nu_tautaubar;
  static int    n_spec;
  static int    IVrho, IVtemp, IVye, IVymu;
  static double logrho_max_IV, logrho_min_IV;
  static double logtemp_max_IV, logtemp_min_IV;
  static double ye_max_IV, ye_min_IV;
  static double logymu_max_IV, logymu_min_IV;

  static bool use_Weakhub_eas;

  // General accessors
  // yle could be yp when num_spec = 3 (no muons):
  static void kappa_ast_kappa_s__temp_rho_yle_ymu(
      double kappa_a_st_en[], double kappa_a_st_num[], double kappa_s[], int num_spec, double &temp,
      double &rho, double &yle, double &ymu, error_type &error);

  //static linear_interp_uniform_ND_t<double, 5, 3> kappa_table;
    // as function of (rho, T, Yp)
  static linear_interp_uniform_ND_t<double, 3, 3> kappa_a_st_en_table_3spec;
  static linear_interp_uniform_ND_t<double, 3, 3> kappa_a_st_num_table_3spec;
  static linear_interp_uniform_ND_t<double, 3, 3> kappa_s_table_3spec;
    // as function of (rho, T, Yle, Ymu)
  static linear_interp_uniform_ND_t<double, 4, 6> kappa_a_st_en_table_6spec;
  static linear_interp_uniform_ND_t<double, 4, 6> kappa_a_st_num_table_6spec;
  static linear_interp_uniform_ND_t<double, 4, 6> kappa_s_table_6spec;
};

#endif
