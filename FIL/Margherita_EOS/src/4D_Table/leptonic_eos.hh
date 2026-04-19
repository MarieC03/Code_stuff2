// fixme: cs2 does not have leptons contribution
//  Copyright (C) 2021, Harry Ho-Yin Ng
//  Based on routines by Elias Roland Most
//  			 Ludwig Jens Papenfort
//
// TODO: Check table bounds!
//        When I need to extend the table like EOS_tabulated::

//#include "cctk.h"
#include <array>
#include <bitset>
#include <iostream>
#include "../margherita.hh"
#include "../utils/brent.hh"
#include "../utils/linear_interp_ND.hh"

#ifndef EOS_LEPTONIC_HH
#define EOS_LEPTONIC_HH

class EOS_Leptonic {
 public:
  // error definitions
  enum errors {
    NO_ERRORS = 0,
    RHOLEP_TOO_HIGH,
    RHOLEP_TOO_LOW,
    TEMPLEP_TOO_HIGH,
    TEMPLEP_TOO_LOW,
    YELEP_TOO_HIGH,
    YELEP_TOO_LOW,
    YMULEP_TOO_HIGH,
    YMULEP_TOO_LOW,
    num_errors
  };
  // enum checked_vars { RHO, TEMP, YE, YMU };
  // return type of errors
  typedef std::bitset<errors::num_errors> error_type;

 //
  enum EMUON {
    MUMU = 0,
    YMU_MINUS,
    YMU_PLUS,
    PRESS_MU_MINUS,
    PRESS_MU_PLUS,
    EPS_MU_MINUS,
    EPS_MU_PLUS,
    S_MU_MINUS,
    S_MU_PLUS,
    NUM_VARS_MUON
  };

  enum EELE {
    MUELE = 0,
    YLE_MINUS,
    YLE_PLUS,
    PRESS_E_MINUS,
    PRESS_E_PLUS,
    EPS_E_MINUS,
    EPS_E_PLUS,
    S_E_MINUS,
    S_E_PLUS,
    NUM_VARS_ELE
  };

  enum EBARYON {
    PRESS = 0,
    EPS,
    S,
    CS2,
    MUE,
    MUP,
    MUN,
    XA,
    XH,
    XN,
    XP,
    ABAR,
    ZBAR,
    NUM_VARS_BARYON
  } ; 

  enum Lepton_type {
    i_ele = 0,
    i_muon,
    NUM_LEPTON
  };
  static constexpr bool has_ye = true;
  static constexpr bool has_ymu = true;
 private:
  template <bool check_temp, bool check_ymu>
  static inline error_type checkbounds(double &xrho, double &xtemp,
                                       double &xyle, double &xymu);

  static inline double find_logtemp_from_eps(const double &lrho, double &eps,
                                             const double yle, const double ymu);

  static inline double find_logtemp_from_entropy(const double &lrho,
                                                 double &entropy,
                                                 const double yle, const double ymu);

  static double press_eps_ye__beta_eq__rho_temp(double& eps, double& ye, double& rho, double& temp, error_type& err) ;
  static double press_eps_ye__nu_beta_eq__rho_temp(double& eps, double& ye, double& rho, double& temp, const double& mu_nu_e_bar, error_type& err) ;
  
 public:
  // Read in leptonic table
  static void readtable_leptonic(const char *leptonic_eos_table_name);
  static void readtable_scollapse(const char *nuceos_table_name,
                                  bool do_energy_shift, bool recompute_mu_nu);
  // Read in compose table
  static void readtable_compose(const char *nuceos_table_name);
  static void generate_cold_table(double logrho[], int &nrho);
  static void test_leptonic_table(double logrho[], int &nrho);

  // General accessors
  static double press__eps_rho_yle_ymu(double &eps, double &rho, double &yle, double &ymu, 
                                       error_type &error);
  static double press_temp__eps_rho_yle_ymu(double &temp, double &eps, double &rho,
                                            double &yle, double &ymu, error_type &error);
  static double press__temp_rho_yle_ymu(double &temp, double &rho, double &yle, double &ymu, 
                                        error_type &error);
  static double eps__temp_rho_yle_ymu(double &temp, double &rho, double &yle, double &ymu, 
                                      error_type &error);
  static double press_cold__rho_yle_ymu(double &rho, double &yle, double &ymu, error_type &error);
  static double eps_cold__rho_yle_ymu(double &rho, double &yle, double &ymu, error_type &error);
  static double temp_cold__rho_yle_ymu(const double &rho, const double &yle, const double &ymu, 
                                       error_type &error);

  static std::array<double, 2> eps_range__rho_yle_ymu(double &rho, double &yle, double &ymu, 
                                                      error_type &error);

  static std::array<double, 2> entropy_range__rho_yle_ymu(double &rho, double &yle, double &ymu, 
                                                      error_type &error);

  // Economical interface. We might want to include also the entropy
  static double press_h_csnd2__eps_rho_yle_ymu(double &h, double &csnd2, double &eps,
                                               double &rho, double &yle, double &ymu, 
                                               error_type &error);
  static double press_h_csnd2__temp_rho_yle_ymu(double &h, double &csnd2,
                                                double &temp, double &rho,
                                                double &yle, double &ymu, error_type &error);
  static double eps_h_csnd2__press_rho_yle_ymu(double &h, double &csnd2,
                                               const double &press, double &rho,
                                               double &yle, double &ymu, error_type &error);

  static double press_eps_csnd2__temp_rho_yle_ymu(double &eps, double &csnd2,
                                                  double &temp, double &rho,
                                                  double &yle, double &ymu, error_type &error);


  // Get everything
  static double press_h_csnd2_temp_entropy__eps_rho_yle_ymu(
      double &h, double &csnd2, double &temp, double &entropy, double &eps,
      double &rho, double &yle, double &ymu, error_type &error);

  static double eps_csnd2_entropy__temp_rho_yle_ymu(double &csnd2, double &entropy,
                                                    double &temp, double &rho,
                                                    double &yle, double &ymu, error_type &error);

  static double press_h_csnd2_temp_eps__entropy_rho_yle_ymu(
      double &h, double &csnd2, double &temp, double &eps, double &entropy,
      double &rho, double &yle, double &ymu, error_type &error);
  static double eps_h_csnd2_temp_entropy__press_rho_yle_ymu(
      double &eps, double &h, double &csnd2, double &temp, double &entropy,
      const double &press, double &rho, double &yle, double &ymu, error_type &error);

  static double press_eps_yle_ymu__beta_eq__rho_temp(double &eps, double &yle, double &ymu,
                                                     double &rho, double &temp,
                                                     error_type &error);

  static double press_eps_yle_ymu__nu_beta_eq__rho_temp(double &eps, double &yle, double &ymu,
                                                        double &rho, double &temp, const double &mu_nu_e_bar, const double &mu_nu_mu_bar, 
                                                        error_type &error);

  static double mue_mumu_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_yle_ymu(
      double &mumu, double &mup, double &mun, double &Xa, double &Xh, double &Xn, double &Xp,
      double &Abar, double &Zbar, double &temp, double &rho, double &yle, double &ymu, 
      error_type &error);

  static double yl_plus_minus_press_l_press_b__temp_rho_yle_ymu(
      double &yle_plus, double &ymu_minus, double &ymu_plus, 
      double &press_e_minus, double &press_e_plus, double &press_mu_minus, double &press_mu_plus, double &press_b, 
      double &eps_e_minus, double &eps_e_plus, double &eps_mu_minus, double &eps_mu_plus, double &eps_b, 
      double &temp, double &rho, double &yle, double &ymu,
      error_type &error);


  static inline double find_yle_of_mumu(const double &lrho, const double &ltemp, const double &mu_mu);
  static inline double find_yle_of_mumu_nubetaeq(const double &lrho, const double &ltemp,
                                                 const double &mu_nu_e_bar, const double &mu_nu_mu_bar, const double &mu_mu);


  static double eos_rhomax, eos_rhomin;
  static double eos_tempmin, eos_tempmax;
  static double eos_ylemax, eos_ylemin;
  // but in hydro, we actually use eos_yemin/max instead of eos_ylemin/max
  static double eos_yemax, eos_yemin;
  static double eos_ymumax, eos_ymumin;
  static double eos_mumumax, eos_mumumin;
  static double eos_muelemax, eos_muelemin;
  static double eos_prs_mu_max, eos_prs_mu_min;
  static double eos_eps_mu_max, eos_eps_mu_min;
  static double eos_prs_ele_max, eos_prs_ele_min;
  static double eos_eps_ele_max, eos_eps_ele_min;
  // need to call beta eqm with muon to find them
  static double energy_shift ;
  static double baryon_mass ;
  static double temp0,temp1 ; 
  // Con2prim limits:
  static double c2p_rho_atm;
  static double c2p_temp_atm;
  static double c2p_eps_atm;
  static double c2p_eps_min;
  static double c2p_h_min;
  static double c2p_h_max;
  static double c2p_eps_max;
  static double c2p_press_max;
  static double c2p_ymu_atm, c2p_ye_atm;
  static bool atm_muonic_beta_eq;
  static bool atm_beta_eq; // FIXME:do I sill need this?

  static bool use_muonic_eos;
  static bool generate_cold_table_local;
  static bool extend_table_high;
  static bool fix_ymu_for_too_high_yp;
  static bool add_ele_contribution;
  static bool presubtract_e_contri;

  // temperature for initial data passing Margherita_ID and ID_converter
  static double temp_ID;

  static linear_interp_uniform_ND_t<double, 3, EOS_Leptonic::EELE::NUM_VARS_ELE> alltables_ele;
  static linear_interp_uniform_ND_t<double, 3, EOS_Leptonic::EMUON::NUM_VARS_MUON> alltables_muon;
  static linear_interp_uniform_ND_t<double, 3, EOS_Leptonic::EBARYON::NUM_VARS_BARYON> alltables_baryon;

#ifdef DEBUG
  static int c2p_roofinding_call_iteration_masfunconly;
#endif

};

#endif
