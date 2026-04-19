// Copyright (C) 2021, Harry Ho-Yin Ng
//    tested : eps to logtemp; beta eqm; error_bounds

#include <array>
#include <bitset>
#include <functional>
#include <iostream>
#include "../M1/fermi.hh"


#include "../utils/brent.hh"
#include "../utils/spline.hh"
#include "leptonic_eos.hh"
//#include "tabulated.hh"

#ifndef EOS_LEPTONIC_IMPL_HH
#define EOS_LEPTONIC_IMPL_HH

// Initialize storage for static members
double EOS_Leptonic::eos_rhomax;
double EOS_Leptonic::eos_rhomin;
double EOS_Leptonic::eos_tempmax;
double EOS_Leptonic::eos_tempmin;
double EOS_Leptonic::eos_ylemax;
double EOS_Leptonic::eos_ylemin;
double EOS_Leptonic::eos_yemax;
double EOS_Leptonic::eos_yemin;
double EOS_Leptonic::eos_ymumin;
double EOS_Leptonic::eos_ymumax;
double EOS_Leptonic::eos_muelemin;
double EOS_Leptonic::eos_muelemax;
double EOS_Leptonic::eos_mumumin;
double EOS_Leptonic::eos_mumumax;
double EOS_Leptonic::eos_prs_mu_min;
double EOS_Leptonic::eos_prs_mu_max;
double EOS_Leptonic::eos_prs_ele_min;
double EOS_Leptonic::eos_prs_ele_max;
double EOS_Leptonic::eos_eps_mu_min;
double EOS_Leptonic::eos_eps_mu_max;
double EOS_Leptonic::eos_eps_ele_min;
double EOS_Leptonic::eos_eps_ele_max;

double EOS_Leptonic::c2p_ymu_atm;
double EOS_Leptonic::c2p_ye_atm;
double EOS_Leptonic::c2p_rho_atm;
double EOS_Leptonic::c2p_temp_atm;
double EOS_Leptonic::c2p_eps_atm;
double EOS_Leptonic::c2p_eps_min;
double EOS_Leptonic::c2p_eps_max;
double EOS_Leptonic::c2p_press_max;
double EOS_Leptonic::c2p_h_min;
double EOS_Leptonic::c2p_h_max;

double EOS_Leptonic::energy_shift ;
double EOS_Leptonic::baryon_mass;
double EOS_Leptonic::temp0; 
double EOS_Leptonic::temp1; 
double EOS_Leptonic::temp_ID; 

bool EOS_Leptonic::atm_muonic_beta_eq = true;
bool EOS_Leptonic::atm_beta_eq = true;
bool EOS_Leptonic::use_muonic_eos = false;
bool EOS_Leptonic::generate_cold_table_local = false;
bool EOS_Leptonic::extend_table_high = false;
bool EOS_Leptonic::fix_ymu_for_too_high_yp = true;
bool EOS_Leptonic::add_ele_contribution = false;
bool EOS_Leptonic::presubtract_e_contri = false;

#ifdef DEBUG
int EOS_Leptonic::c2p_roofinding_call_iteration_masfunconly = 0;
#endif

linear_interp_uniform_ND_t<double, 3, EOS_Leptonic::EELE::NUM_VARS_ELE>
    EOS_Leptonic::alltables_ele;
linear_interp_uniform_ND_t<double, 3, EOS_Leptonic::EMUON::NUM_VARS_MUON>
    EOS_Leptonic::alltables_muon;
linear_interp_uniform_ND_t<double, 3, EOS_Leptonic::EBARYON::NUM_VARS_BARYON>
    EOS_Leptonic::alltables_baryon;

// ###################################################################################
// 	Private functions
// ###################################################################################

//  error and bounds checkers
template <bool check_temp, bool check_ymu>
EOS_Leptonic::error_type EOS_Leptonic::checkbounds(double &xrho,
                                                   double &xtemp,
                                                   double &xyle, 
                                                   double &xymu) {
  // initialize error codes to zero
  error_type error_codes;
  if (xrho > eos_rhomax && !extend_table_high) {
    // This happens only inside the AH.
    // (assuming that the table has a sane range)
    error_codes[errors::RHOLEP_TOO_HIGH] = true;
    xrho = eos_rhomax;
  } else if (xrho < eos_rhomin) {
    xrho = eos_rhomin;
    error_codes[errors::RHOLEP_TOO_LOW] = true;
  }
  if (xyle > eos_ylemax) {
    xyle = eos_ylemax;
    error_codes[errors::YELEP_TOO_HIGH] = true;
  } else if (xyle < eos_ylemin) {
    xyle = eos_ylemin;
    error_codes[errors::YELEP_TOO_LOW] = true;
  }
  if ( check_ymu) {
    if (xymu > eos_ymumax) {
      xymu = eos_ymumax;
      error_codes[errors::YMULEP_TOO_HIGH] = true;
      //abort();
      //CCTK_ERROR("eos: xymu > eos_ymumax, stop");
    } else if (xymu < eos_ymumin) {
      xymu = eos_ymumin;
      error_codes[errors::YMULEP_TOO_LOW] = true;
    }

    // these bounds --> NaN for Ymu in c2p
    //if (log(xymu) > log(eos_ymumax)) {
    //  xymu = eos_ymumax;
    //  error_codes[errors::YMULEP_TOO_HIGH] = true;
    //} else if (log(xymu) < log(eos_ymumin)) {
    //  xymu = eos_ymumin;
    //  error_codes[errors::YMULEP_TOO_LOW] = true;
    //}
  }
  if (check_temp) {
    if (xtemp > eos_tempmax) {
      xtemp = eos_tempmax;
      error_codes[errors::TEMPLEP_TOO_HIGH] = true;
    } else if (xtemp < eos_tempmin) {
      xtemp = eos_tempmin;
      error_codes[errors::TEMPLEP_TOO_LOW] = true;
    }
  }
  return error_codes;
}
//Dont know how to do this yet.
inline double EOS_Leptonic::find_logtemp_from_eps(const double &lrho,
                                                  double &eps,
                                                  const double yle, const double ymu) {
  //auto const leps = log(eps + energy_shift);

  double eps_local = eps;

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  const double lymu = log(ymu);

  // Get eps ranges
  auto const varsmin =
      alltables_baryon.interpolate<EBARYON::EPS>(lrho, alltables_baryon.xmin<1>(), yp);
  auto const varsmax =
      alltables_baryon.interpolate<EBARYON::EPS>(lrho, alltables_baryon.xmax<1>(), yp);

  auto const varsmin_muon =
      alltables_muon.interpolate<EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, alltables_baryon.xmin<1>(), lymu);
  auto const varsmax_muon =
      alltables_muon.interpolate<EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, alltables_baryon.xmax<1>(), lymu);

  double eps_lower = 0.0;
  double eps_upper = 0.0;
 
  if (add_ele_contribution) {
	  auto const varsmin_ele =
	      alltables_ele.interpolate<EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, alltables_baryon.xmin<1>(), yle);
	  auto const varsmax_ele =
	      alltables_ele.interpolate<EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, alltables_baryon.xmax<1>(), yle);

          eps_lower = exp(varsmin[0]) - energy_shift +
                  varsmin_ele[0] + varsmin_ele[1] + varsmin_muon[0] + varsmin_muon[1];

          eps_upper = exp(varsmax[0]) - energy_shift +
                  varsmax_ele[0] + varsmax_ele[1] + varsmax_muon[0] + varsmax_muon[1];
  } else {
          eps_lower = exp(varsmin[0]) - energy_shift +
                   varsmin_muon[0] + varsmin_muon[1];

          eps_upper = exp(varsmax[0]) - energy_shift +
                   varsmax_muon[0] + varsmax_muon[1];
  }


  if (eps_local <= eps_lower) {
    eps = eps_lower;
    return alltables_baryon.xmin<1>();
  }
  if (eps_local >= eps_upper ) {
    eps = eps_upper;
    return alltables_baryon.xmax<1>();
  }

  // Root finding interface closure
  auto const func = [&](double &lt) {
    auto const vars = alltables_baryon.interpolate<EBARYON::EPS>(lrho, lt, yp);

    auto const vars_muon = alltables_muon.interpolate<EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, lt, lymu);
    if (add_ele_contribution) {
   	 auto const vars_ele = alltables_ele.interpolate<EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, lt, yle);
    	 return eps_local - (exp(vars[0]) - energy_shift + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1]);
    } else {
    	 return eps_local - (exp(vars[0]) - energy_shift + vars_muon[0] + vars_muon[1]);
    }
  };

  //return zero_brent(alltables_baryon.xmin<1>(), alltables_baryon.xmax<1>(), 1.e-12, func);
  return zero_brent(alltables_baryon.xmin<1>(), alltables_baryon.xmax<1>(), 1.e-14, func);
}

//Dont know how to do this yet.
inline double EOS_Leptonic::find_logtemp_from_entropy(const double &lrho,
                                                       double &entropy,
                                                       const double yle, const double ymu) {

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  const double lymu = log(ymu);

  // Get eps ranges
  auto const varsmin =
      alltables_baryon.interpolate<EBARYON::S>(lrho, alltables_baryon.xmin<1>(), yp);
  auto const varsmax =
      alltables_baryon.interpolate<EBARYON::S>(lrho, alltables_baryon.xmax<1>(), yp);

  auto const varsmin_muon =
      alltables_muon.interpolate<EMUON::S_MU_MINUS, EMUON::S_MU_PLUS>(lrho, alltables_baryon.xmin<1>(), lymu);
  auto const varsmax_muon =
      alltables_muon.interpolate<EMUON::S_MU_MINUS, EMUON::S_MU_PLUS>(lrho, alltables_baryon.xmax<1>(), lymu);

  double entropy_lower = 0.0;
  double entropy_upper = 0.0;

  if (add_ele_contribution) {
	auto const varsmin_ele =
	    alltables_ele.interpolate<EELE::S_E_MINUS, EELE::S_E_PLUS>(lrho, alltables_baryon.xmin<1>(), yle);
	auto const varsmax_ele =
	    alltables_ele.interpolate<EELE::S_E_MINUS, EELE::S_E_PLUS>(lrho, alltables_baryon.xmax<1>(), yle);
  	entropy_lower = varsmin[0] + 
  	        varsmin_ele[0] + varsmin_ele[1] + varsmin_muon[0] + varsmin_muon[1];

  	entropy_upper = varsmax[0] +
  	        varsmax_ele[0] + varsmax_ele[1] + varsmax_muon[0] + varsmax_muon[1];
  } else {
  	entropy_lower = varsmin[0] + 
  	        varsmin_muon[0] + varsmin_muon[1];

  	entropy_upper = varsmax[0] +
  	        varsmax_muon[0] + varsmax_muon[1];
  }

  if (entropy <= entropy_lower) {
    entropy = entropy_lower;
    return alltables_baryon.xmin<1>();
  }
  if (entropy >= entropy_upper) {
    entropy = entropy_upper;
    return alltables_baryon.xmax<1>();
  }

  // Root finding interface closure
  auto const func = [&](double &lt) {
    auto const vars = alltables_baryon.interpolate<EBARYON::S>(lrho, lt, yp);

    auto const vars_muon = alltables_muon.interpolate<EMUON::S_MU_MINUS, EMUON::S_MU_PLUS>(lrho, lt, lymu);
    if (add_ele_contribution) {
   	 auto const vars_ele = alltables_ele.interpolate<EELE::S_E_MINUS, EELE::S_E_PLUS>(lrho, lt, yle);
   	 return entropy - (vars[0] + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1]);
    } else {
   	 return entropy - (vars[0] + vars_muon[0] + vars_muon[1]);
    }
  };

  return zero_brent(alltables_baryon.xmin<1>(), alltables_baryon.xmax<1>(), 1.e-14, func);
}

// ###################################################################################
// 	Public functions
// ###################################################################################
double EOS_Leptonic::press__eps_rho_yle_ymu(double &eps, double &rho, double &yle, double &ymu, 
                                        error_type &error) {
  // Check if rho and Y_e lie inside the table, otherwise abort!
  double temp_tmp = 0;
  error = checkbounds<false,true>(rho, temp_tmp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) { 
        ymu = eos_ymumin;
     }
  }

  const double lrho = log(rho);
  const double ltemp = find_logtemp_from_eps(lrho, eps, yle, ymu);
  const double lymu = log(ymu);

  double press;

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars = alltables_baryon.interpolate<EBARYON::PRESS>(lrho, ltemp, yp);
  auto const vars_muon = alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS>(lrho, ltemp, lymu);

  if (add_ele_contribution) {
  	auto const vars_ele  = alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS>(lrho, ltemp, yle);
  	return exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  } else {
  	return exp(vars[0]) + vars_muon[0] + vars_muon[1];
  }
}

double EOS_Leptonic::press_temp__eps_rho_yle_ymu(double &temp, double &eps,
                                                 double &rho, double &yle, double &ymu,
                                                 error_type &error) {
  // Check if rho and Y_e lie inside the table, otherwise abort!
  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }


  const double lrho = log(rho);
  const double ltemp = find_logtemp_from_eps(lrho, eps, yle, ymu);
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  temp = exp(ltemp);

  auto const vars = alltables_baryon.interpolate<EBARYON::PRESS>(lrho, ltemp, yp);
  auto const vars_muon = alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS>(lrho, ltemp, lymu);

  if (add_ele_contribution) {
  	auto const vars_ele  = alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS>(lrho, ltemp, yle);
  	return exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  } else {
  	return exp(vars[0]) + vars_muon[0] + vars_muon[1];
  }

}

double EOS_Leptonic::press__temp_rho_yle_ymu(double &temp, double &rho, double &yle, double &ymu, 
                                                 error_type &error) {
  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const auto lrho = log(rho);
  const auto ltemp = log(temp);
  const auto lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars = alltables_baryon.interpolate<EBARYON::PRESS>(lrho, ltemp, yp);
  auto const vars_muon = alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS>(lrho, ltemp, lymu);

  if (add_ele_contribution) {
  	auto const vars_ele  = alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS>(lrho, ltemp, yle);
  	return exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  } else {
  	return exp(vars[0]) + vars_muon[0] + vars_muon[1];
  }
}

double EOS_Leptonic::eps__temp_rho_yle_ymu(double &temp, double &rho, double &yle, double &ymu, 
                                           error_type &error) {
  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }


  const double lrho = log(rho);
  const double ltemp = log(temp);
  const auto lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars = alltables_baryon.interpolate<EBARYON::EPS>(lrho, ltemp, yp);
  auto const vars_muon = alltables_muon.interpolate<EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);

  if (add_ele_contribution) {
  	auto const vars_ele  = alltables_ele.interpolate<EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  	return exp(vars[0]) - energy_shift + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  } else {
  	return exp(vars[0]) - energy_shift + vars_muon[0] + vars_muon[1];
  }
}

double EOS_Leptonic::press_cold__rho_yle_ymu(double &rho, double &yle, double &ymu,
                                         error_type &error) {
  double temp_tmp = 0;
  error = checkbounds<false,true>(rho, temp_tmp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double lrho = log(rho);
  const double ltemp = log(eos_tempmin);  
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars = alltables_baryon.interpolate<EBARYON::PRESS>(lrho, ltemp, yp);
  auto const vars_muon = alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS>(lrho, ltemp, lymu);

  if (add_ele_contribution) {
  	auto const vars_ele  = alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS>(lrho, ltemp, yle);
  	return exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  } else {
  	return exp(vars[0]) + vars_muon[0] + vars_muon[1];
  }
}

double EOS_Leptonic::eps_cold__rho_yle_ymu(double &rho, double &yle, double &ymu,
                                           error_type &error) {
  double temp_tmp = 0;
  error = checkbounds<false,true>(rho, temp_tmp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double lrho = log(rho);
  const double ltemp = log(eos_tempmin);
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars = alltables_baryon.interpolate<EBARYON::EPS>(lrho, ltemp, yp);
  auto const vars_muon = alltables_muon.interpolate<EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);

  if (add_ele_contribution) {
  	auto const vars_ele  = alltables_ele.interpolate<EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  	return exp(vars[0]) - energy_shift + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  } else {
  	return exp(vars[0]) - energy_shift + vars_muon[0] + vars_muon[1];
  }
}

double EOS_Leptonic::temp_cold__rho_yle_ymu(const double &rho, const double &yle, const double &ymu, 
                                            error_type &error) {
  //  for (auto i = 0; i < error.size(); ++i)
  //    error[i] = errors::NO_ERRORS;
  //return temp0;
  return eos_tempmin;
}
///////////

std::array<double, 2> EOS_Leptonic::eps_range__rho_yle_ymu(double &rho, double &yle,double &ymu, 
                                                           error_type &error) {
  double temp_tmp = 0;
  error = checkbounds<false,true>(rho, temp_tmp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double lrho = log(rho);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  double lymu = log(ymu);

  auto const varsmin =
      alltables_baryon.interpolate<EBARYON::EPS>(lrho, alltables_baryon.xmin<1>(), yp);
  auto const varsmax =
      alltables_baryon.interpolate<EBARYON::EPS>(lrho, alltables_baryon.xmax<1>(), yp);

  auto const varsmin_muon =
      alltables_muon.interpolate<EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, alltables_baryon.xmin<1>(), lymu);
  auto const varsmax_muon =
      alltables_muon.interpolate<EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, alltables_baryon.xmax<1>(), lymu);

  if (add_ele_contribution) {
  	auto const varsmin_ele =
  	    alltables_ele.interpolate<EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, alltables_baryon.xmin<1>(), yle);
  	auto const varsmax_ele =
  	    alltables_ele.interpolate<EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, alltables_baryon.xmax<1>(), yle);
  	return std::array<double, 2>{exp(varsmin[0]) - energy_shift + 
  	                               varsmin_ele[0] + varsmin_ele[1] + varsmin_muon[0] + varsmin_muon[1],
  	                             exp(varsmax[0]) - energy_shift + 
  	                               varsmax_ele[0] + varsmax_ele[1] + varsmax_muon[0] + varsmax_muon[1]};
  } else {
        return std::array<double, 2>{exp(varsmin[0]) - energy_shift +
                                       varsmin_muon[0] + varsmin_muon[1],
                                     exp(varsmax[0]) - energy_shift +
                                        varsmax_muon[0] + varsmax_muon[1]};
  }

}

std::array<double, 2> EOS_Leptonic::entropy_range__rho_yle_ymu(double &rho,
                                                               double &yle,
                                                               double &ymu,
                                                               error_type &error) {
  double temp_tmp = 0;
  error = checkbounds<false,true>(rho, temp_tmp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double lrho = log(rho);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  double lymu = log(ymu);

  auto const varsmin =
      alltables_baryon.interpolate<EBARYON::S>(lrho, alltables_baryon.xmin<1>(), yp);
  auto const varsmax =
      alltables_baryon.interpolate<EBARYON::S>(lrho, alltables_baryon.xmax<1>(), yp);

  auto const varsmin_muon =
      alltables_muon.interpolate<EMUON::S_MU_MINUS, EMUON::S_MU_PLUS>(lrho, alltables_baryon.xmin<1>(), lymu);
  auto const varsmax_muon =
      alltables_muon.interpolate<EMUON::S_MU_MINUS, EMUON::S_MU_PLUS>(lrho, alltables_baryon.xmax<1>(), lymu);

  if (add_ele_contribution) {
  	auto const varsmin_ele =
  	    alltables_ele.interpolate<EELE::S_E_MINUS, EELE::S_E_PLUS>(lrho, alltables_baryon.xmin<1>(), yle);
  	auto const varsmax_ele =
  	    alltables_ele.interpolate<EELE::S_E_MINUS, EELE::S_E_PLUS>(lrho, alltables_baryon.xmax<1>(), yle);
  	return std::array<double, 2>{varsmin[0] +
  	                               varsmin_ele[0] + varsmin_ele[1] + varsmin_muon[0] + varsmin_muon[1],
  	                             varsmax[0] +
  	                               varsmax_ele[0] + varsmax_ele[1] + varsmax_muon[0] + varsmax_muon[1]};
  } else {
  	return std::array<double, 2>{varsmin[0] +
  	                                varsmin_muon[0] + varsmin_muon[1],
  	                             varsmax[0] +
  	                                varsmax_muon[0] + varsmax_muon[1]};
  }


}
/////////////////
// the eps here is total eps
double EOS_Leptonic::press_h_csnd2__eps_rho_yle_ymu(double &h, double &csnd2,
                                                    double &eps, double &rho,
                                                    double &yle, double &ymu, error_type &error) {
  using namespace Margherita_helpers;

  double temp_tmp = 0;
  error = checkbounds<false,true>(rho, temp_tmp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const auto lrho = log(rho);
  const auto ltemp = find_logtemp_from_eps(lrho, eps, yle, ymu);
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars = alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::CS2>(lrho, ltemp, yp);

  auto const vars_muon = alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS>(lrho, ltemp, lymu);

  double press = 0.0;
  if (add_ele_contribution) {
  	auto const vars_ele  = alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS>(lrho, ltemp, yle);
  	press = exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  } else {
  	press = exp(vars[0]) + vars_muon[0] + vars_muon[1];
  }

  h = 1. + eps + press / rho;

  csnd2 = vars[1];  // min(vars[1] / h, max_csnd2);

  return press;
};

double EOS_Leptonic::press_h_csnd2__temp_rho_yle_ymu(double &h, double &csnd2,
                                                     double &temp, double &rho,
                                                     double &yle, double &ymu,
                                                     error_type &error) {
  using namespace Margherita_helpers;

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double ltemp = log(temp);
  const double lrho = log(rho);
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars =
    alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::CS2, EBARYON::EPS>(lrho, ltemp, yp);

  auto const vars_muon = 
    alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS, EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);
  double press = 0.0;
  double eps   = 0.0;

  if (add_ele_contribution) {
  	auto const vars_ele  = 
  	  alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS, EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  	press = exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[2]) - energy_shift + vars_ele[2] + vars_ele[3] + vars_muon[2] + vars_muon[3];
  } else {
  	press = exp(vars[0]) + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[2]) - energy_shift + vars_muon[2] + vars_muon[3];
  }

  h = 1. + eps + press / rho;

  csnd2 = vars[1];  // min(vars[1] / h, max_csnd2);

  return press;
}

double EOS_Leptonic::eps_h_csnd2__press_rho_yle_ymu(double &h, double &csnd2,
                                                    const double &press,
                                                    double &rho, double &yle, double &ymu, 
                                                    error_type &error) {
  using namespace Margherita_helpers;

  assert(!"This routine should not be used. There is no monotonicity condition "
          "to enforce a succesfull inversion from eps(press). So you better "
          "rewrite your code to not require this call... (EOS_Leptonic)");

  return 0;
}

double EOS_Leptonic::press_eps_csnd2__temp_rho_yle_ymu(double &eps, double &csnd2,
                                                       double &temp, double &rho,
                                                       double &yle, double &ymu, 
                                                       error_type &error) {
  using namespace Margherita_helpers;

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double ltemp = log(temp);
  const double lrho = log(rho);
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars =
    alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::CS2, EBARYON::EPS>(lrho, ltemp, yp);

  auto const vars_muon =
    alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS, EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);
  double press = 0.0;
  if (add_ele_contribution) {
  	auto const vars_ele  =
  	  alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS, EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  	press = exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[2]) - energy_shift + vars_ele[2] + vars_ele[3] + vars_muon[2] + vars_muon[3];
  } else {
  	press = exp(vars[0]) + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[2]) - energy_shift + vars_muon[2] + vars_muon[3];
  }

  csnd2 = vars[1];

  return press;
}

double EOS_Leptonic::press_h_csnd2_temp_entropy__eps_rho_yle_ymu(
                       double &h, double &csnd2, double &temp, double &entropy, double &eps,
                       double &rho, double &yle, double &ymu, error_type &error) {

  using namespace Margherita_helpers;

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double lrho = log(rho);
  const auto ltemp = find_logtemp_from_eps(lrho, eps, yle, ymu);
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars =
    alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::CS2, EBARYON::S>(lrho, ltemp, yp);

  auto const vars_muon =
    alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS, EMUON::S_MU_MINUS, EMUON::S_MU_PLUS>(lrho, ltemp, lymu);

  double press = 0.0;
  if (add_ele_contribution) {
  	auto const vars_ele  =
  	  alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS, EELE::S_E_MINUS, EELE::S_E_PLUS>(lrho, ltemp, yle);
  	press = exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  	entropy = vars[2] + vars_ele[2] + vars_ele[3] + vars_muon[2] + vars_muon[3];
  } else {
  	press = exp(vars[0]) + vars_muon[0] + vars_muon[1];
  	entropy = vars[2] + vars_muon[2] + vars_muon[3];
  }

  temp = exp(ltemp);

  h = 1. + eps + press / rho;

  csnd2 = vars[1];  // min(vars[1] / h, max_csnd2);

  return press;
}

double EOS_Leptonic::eps_csnd2_entropy__temp_rho_yle_ymu(double &csnd2,
                                                         double &entropy,
                                                         double &temp, double &rho,
                                                         double &yle, double &ymu, 
                                                         error_type &error) {
  using namespace Margherita_helpers;

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double lrho = log(rho);
  const double ltemp = log(temp);
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars =
    alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::CS2, EBARYON::S, EBARYON::EPS>(lrho, ltemp, yp);

  auto const vars_muon =
    alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS, EMUON::S_MU_MINUS, EMUON::S_MU_PLUS, 
                               EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);

  double press = 0.0;
  double eps   = 0.0;

  if (add_ele_contribution) {
  	auto const vars_ele  =
  	  alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS, EELE::S_E_MINUS, EELE::S_E_PLUS, 
                              		EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  	press = exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  	entropy = vars[2] + vars_ele[2] + vars_ele[3] + vars_muon[2] + vars_muon[3];
  	eps = exp(vars[3]) - energy_shift + vars_ele[4] + vars_ele[5] + vars_muon[4] + vars_muon[5];
  } else {
  	press = exp(vars[0]) + vars_muon[0] + vars_muon[1];
  	entropy = vars[2] + vars_muon[2] + vars_muon[3];
  	eps = exp(vars[3]) - energy_shift + vars_muon[4] + vars_muon[5];
  }

  const double h = 1. + eps + press / rho;

  csnd2 = vars[1];  // min(vars[1] / h, max_csnd2);

  return eps;
}

double EOS_Leptonic::press_h_csnd2_temp_eps__entropy_rho_yle_ymu(
    double &h, double &csnd2, double &temp, double &eps, double &entropy,
    double &rho, double &yle, double &ymu, error_type &error) {

  using namespace Margherita_helpers;

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const auto lrho = log(rho);
  const auto ltemp = find_logtemp_from_entropy(lrho, entropy, yle, ymu);
  const double lymu = log(ymu);

  double yp;
  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars =
    alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::CS2, EBARYON::EPS>(lrho, ltemp, yp);

  auto const vars_muon =
    alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS,
                               EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);

  double press = 0.0;
  if (add_ele_contribution) {
  	auto const vars_ele  =
  	  alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS,
  	                            EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  	press = exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[2]) - energy_shift + vars_ele[2] + vars_ele[3] + vars_muon[2] + vars_muon[3];
  } else {
  	press = exp(vars[0]) + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[2]) - energy_shift + vars_muon[2] + vars_muon[3];
  }

  temp = exp(ltemp);

  h = 1. + eps + press / rho;

  csnd2 = vars[1];  // min(vars[1] / h, max_csnd2);

  return press;
}

double EOS_Leptonic::eps_h_csnd2_temp_entropy__press_rho_yle_ymu(
    double &eps, double &h, double &csnd2, double &temp, double &entropy,
    const double &press, double &rho, double &yle, double &ymu, error_type &error) {
  using namespace Margherita_helpers;

  assert(!"This routine should not be used. There is no monotonicity condition "
          "to enforce a succesfull inversion from eps(press). So you better "
          "rewrite your code to not require this call...(Leptonic)");

  return 0;
}
// leptonic mu included rest mass, nucleonic mu do not.
// We do not use the mue inside nucleonic table anymore
double EOS_Leptonic::mue_mumu_mup_mun_Xa_Xh_Xn_Xp_Abar_Zbar__temp_rho_yle_ymu(
    double &mumu, double &mup, double &mun, double &Xa, double &Xh, double &Xn, double &Xp,
    double &Abar, double &Zbar, double &temp, double &rho, double &yle, double &ymu, 
    error_type &error) {

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }

  const double lrho = log(rho);
  const double ltemp = log(temp);
  const double lymu = log(ymu);

  double yp, yp_local, ymu_local, yle_local;
  yp = yle + ymu;
  yp_local = yp;
  ymu_local = ymu;
  yle_local = yle;
 
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars =
      alltables_baryon.interpolate<EBARYON::MUE, EBARYON::MUP, 
                            EBARYON::MUN, EBARYON::XA, EBARYON::XH, 
                            EBARYON::XN, EBARYON::XP, EBARYON::ABAR, 
                            EBARYON::ZBAR>(lrho, ltemp, yp);
  
  mup = vars[1];
  mun = vars[2];

  Xa = vars[3];
  Xh = vars[4];
  Xn = vars[5];
  Xp = vars[6];

  Abar = vars[7];
  Zbar = vars[8];

  // Prevent corrupted mu_e inside leptonic table, we use baryon table (rho, yle, T)
  if (yp_local >= eos_yemax) {
        // in this case, we use baryonic table mu_e and force mu_mu = 0
  	mumu = 105.6583755; // mu_mu at atmosphere value is not 0 mev
  	return vars[0];  // mue
  } else {
  	auto const vars_muon = alltables_muon.interpolate<EMUON::MUMU>(lrho, ltemp, lymu);
        mumu = vars_muon[0];
  	auto const vars_ele = alltables_ele.interpolate<EELE::MUELE>(lrho, ltemp, yle);
  	return vars_ele[0];  // mue
  }
  // Original leptonic tables
  //auto const vars_muon = alltables_muon.interpolate<EMUON::MUMU>(lrho, ltemp, lymu);
  //auto const vars_ele = alltables_ele.interpolate<EELE::MUELE>(lrho, ltemp, yle);
  //mumu = vars_muon[0];
  //return vars_ele[0];
};

double EOS_Leptonic::yl_plus_minus_press_l_press_b__temp_rho_yle_ymu(
      double &yle_plus, double &ymu_minus, double &ymu_plus,
      double &press_e_minus, double &press_e_plus, double &press_mu_minus, double &press_mu_plus, double &press_b,
      double &eps_e_minus, double &eps_e_plus, double &eps_mu_minus, double &eps_mu_plus, double &eps_b,
      double &temp, double &rho, double &yle, double &ymu,
      error_type &error) {

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  if (fix_ymu_for_too_high_yp) {
     if (yle + ymu > eos_yemax) {
        ymu = eos_ymumin;
     }
  }


  const double lrho = log(rho);
  const double ltemp = log(temp);
  const double lymu = log(ymu);

  double yp, yp_local, ymu_local, yle_local;
  yp = yle + ymu;
  yp_local = yp;
  ymu_local = ymu;
  yle_local = yle;
 
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  auto const vars =
      alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::EPS>(lrho, ltemp, yp);
  
  press_b = exp(vars[0]);
  eps_b   = exp(vars[1]) - energy_shift;

  auto const vars_ele = alltables_ele.interpolate<EELE::YLE_MINUS, EELE::YLE_PLUS, EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS, 
                                                  EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  yle_plus      = vars_ele[1];
  press_e_minus = vars_ele[2];
  press_e_plus  = vars_ele[3];
  eps_e_minus   = vars_ele[4];
  eps_e_plus    = vars_ele[5];

  auto const vars_muon = alltables_muon.interpolate<EMUON::YMU_MINUS, EMUON::YMU_PLUS, EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS,
                                                    EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);

  ymu_minus      = vars_muon[0];
  ymu_plus       = vars_muon[1];
  press_mu_minus = vars_muon[2];
  press_mu_plus  = vars_muon[3];
  eps_mu_minus   = vars_muon[4];
  eps_mu_plus    = vars_muon[5];


  //if (yp_local >= eos_yemax) {
  //      press_mu_minus = 0.0;
  //      press_mu_plus  = 0.0;
  //      eps_mu_minus   = 0.0;
  //      eps_mu_plus    = 0.0;
  //}
  return vars_ele[0]; // return yle_minus
};


// Two constraints: 1. mu_e = mu_mu
//                  2. mu_n - mu_p - mu_mu = 0
double EOS_Leptonic::press_eps_yle_ymu__beta_eq__rho_temp(double &eps, double &yle,
                                              double &ymu, double &rho, double &temp, 
                                              error_type &error){

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  const auto lrho = log(rho);
  const auto ltemp = log(temp);

  // Rootfinding f(Y_mu)
  auto const func_y = [&](double &lymu_in) {

    auto const mu_mu = alltables_muon.interpolate<EMUON::MUMU>(lrho, ltemp, lymu_in);

    // Here input mu_mu to find yle_local
    auto const yle_local = find_yle_of_mumu(lrho,ltemp,mu_mu[0]);
    
    // Charge Neutrality:
    double yp = exp(lymu_in) + yle_local;
    yp = std::min(eos_yemax,std::max(yp, eos_yemin));

    // get mu_p and mu_n as function of (rho, T, yp = yle+ymu)
    auto const mu_nucleons = alltables_baryon.interpolate<EBARYON::MUP, EBARYON::MUN>(lrho, ltemp, yp);

    yle = yle_local;
 
    // Satisfying mu_n - mu_p - mu_mu = 0
    return mu_mu[0] + mu_nucleons[0] - mu_nucleons[1] - M1_Constants::Qnp;

  };

  // ln scale of ymu and then change it back to fraction value
  ymu = zero_brent(log(eos_ymumin), log(eos_ymumax), 1.e-14, func_y);

  double yp;
  // Now ymu is log(Y_mu)
  // log(Y_mu) <= log(eos_ymumin) is more accurate than Y_mu <= eos_ymumin
  // The floating error with 5e-5 is good to capture correct Yp when T = 0.1 MeV
  // Becoz there should not be muon gas in atmosphere and stellar medium
  if (ymu <= log(eos_ymumin + 0.00001)) {
  //if (ymu <= log(eos_ymumin)) {
    //std::cout << "passing <= eos_ymumin" << rho<<" "<<yle<<" "<<ymu<<" \n";
    // Fix ymu to be lower bound, since log and exp induced floating pt error
    ymu = eos_ymumin;

    double dummy;
    error_type error2;
    // use pure npe-eqm instead
    press_eps_ye__beta_eq__rho_temp(dummy, yp, rho, temp, error2);
    yle = yp - ymu;
  } else {  
    ymu = exp(ymu);
  }

  yp = yle + ymu;
  // bound yp for baryonic table:
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));

  const auto lymu = log(ymu);

  auto const vars =
    alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::EPS>(lrho, ltemp, yp);

  auto const vars_muon =
    alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS,
                               EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);

  double press = 0.0;
  if (add_ele_contribution) {
  	auto const vars_ele  =
  	  alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS,
  	                            EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  	press = exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[1]) - energy_shift + vars_ele[2] + vars_ele[3] + vars_muon[2] + vars_muon[3];
  } else {
  	press = exp(vars[0]) + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[1]) - energy_shift + vars_muon[2] + vars_muon[3];
  }
  return press;
};

inline double EOS_Leptonic::find_yle_of_mumu(const double &lrho, const double &ltemp, const double &mu_mu){
  double yle_local2;
  auto const func_yle_of_mumu = [&](double &yle_local2) {
   auto const mu_e =
//FIXME should I change here for without add_ele_contribution
     alltables_ele.interpolate<EELE::MUELE>(lrho, ltemp, yle_local2);
   return mu_e[0] - mu_mu;
  };

  yle_local2 = zero_brent(eos_ylemin, eos_ylemax, 1.e-14, func_yle_of_mumu);
  return yle_local2;
}

double EOS_Leptonic::press_eps_ye__beta_eq__rho_temp(double &eps, double &ye,
						     double &rho, double &temp,
						     error_type &error)
{
  double dummy{} ;
  // Check if rho and Y_e lie inside the table, otherwise abort!
  error = checkbounds<true,false>(rho, temp, ye, dummy);

  const double lrho = log(rho);
  const double ltemp = log(temp);

  // Beta equilibrium requires that \mu_n =  \mu_p +\mu_e -\mu_nu
  // Since we assume that \mu_nu should be negligible we demand
  // \mu_n-\mu_p-\mu_e = \mu_hat =0

  // Root finding interface closure
  auto const func = [&](double &ye) {
    auto const vars =
        alltables_baryon.interpolate<EBARYON::MUE, EBARYON::MUP, EBARYON::MUN>(lrho, ltemp, ye);
    return vars[0] + vars[1] - vars[2] - M1_Constants::Qnp;
    // Qnp corrections make a significant difference
    //return vars[0] + vars[1] - vars[2];
  };
  

  ye = zero_brent(alltables_baryon.xmin<2>(), alltables_baryon.xmax<2>(), 1.e-14, func);

  auto const vars = alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::EPS>(lrho, ltemp, ye);

  eps = exp(vars[1]) - energy_shift;
  return exp(vars[0]);
}

double EOS_Leptonic::press_eps_ye__nu_beta_eq__rho_temp(double &eps, double &ye,
							double &rho, double &temp, const double &mu_nu_e_bar, 
							error_type &error) {
  double dummy{} ;
  // Check if rho and Y_e lie inside the table, otherwise abort!
  error = checkbounds<true,false>(rho, temp, ye,dummy);

  const double lrho = log(rho);
  const double ltemp = log(temp);

  // Beta equilibrium requires that \mu_n =  \mu_p +\mu_e -\mu_nu
  // Since we assume that \mu_nu should be negligible we demand
  // \mu_n-\mu_p-\mu_e = \mu_hat =0

  // Root finding interface closure
  auto const func = [&](double &ye) {
    auto const vars =
        alltables_baryon.interpolate<EBARYON::MUE, EBARYON::MUP, EBARYON::MUN>(lrho, ltemp, ye);
    return vars[0] + vars[1] - vars[2] - M1_Constants::Qnp - mu_nu_e_bar;
  };

  ye = zero_brent(alltables_baryon.xmin<2>(), alltables_baryon.xmax<2>(), 1.e-14, func);

  auto const vars = alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::EPS>(lrho, ltemp, ye);

  eps = exp(vars[1]) - energy_shift;
  return exp(vars[0]);
}


// STILL NEED TO TEST FIXME:
// Neutrino muonic beta eqm:
// Two constraints: 1. mu_e + mu_nu_e_bar = mu_mu + mu_nu_mu_bar
//                  2. mu_n - mu_p - mu_mu - mu_nu_e_bar = 0
double EOS_Leptonic::press_eps_yle_ymu__nu_beta_eq__rho_temp(double &eps, double &yle,
                                                             double &ymu, double &rho, double &temp, const double &mu_nu_e_bar, const double &mu_nu_mu_bar, 
                                                             error_type &error){

  error = checkbounds<true,true>(rho, temp, yle, ymu);

  const auto lrho = log(rho);
  const auto ltemp = log(temp);

  // Should not apply neutrino muonic beta eq in these regions
  if (rho/Margherita_constants::RHOGF < 1.0e11) {
     assert(!"rho/RHOGF < 1e11 joining to muonic nu beta eq, wrong");
     return 0;
  }

  // Rootfinding f(Y_mu)
  auto const func_y = [&](double &lymu_in) {

    auto const mu_mu = alltables_muon.interpolate<EMUON::MUMU>(lrho, ltemp, lymu_in);

    // Here input mu_mu to find yle_local
    auto const yle_local = find_yle_of_mumu_nubetaeq(lrho,ltemp,mu_nu_e_bar,mu_nu_mu_bar,mu_mu[0]);

    // Charge Neutrality:
    double yp = exp(lymu_in) + yle_local;
    yp = std::min(eos_yemax,std::max(yp, eos_yemin));

    // get mu_p and mu_n as function of (rho, T, yp = yle+ymu)
    auto const mu_nucleons = alltables_baryon.interpolate<EBARYON::MUP, EBARYON::MUN>(lrho, ltemp, yp);

    yle = yle_local;

    // Satisfying mu_n - mu_p - mu_mu = 0
    return mu_mu[0] + mu_nucleons[0] - mu_nucleons[1] - M1_Constants::Qnp - mu_nu_mu_bar;

  };

  // ln scale of ymu and then change it back to fraction value
  ymu = zero_brent(log(eos_ymumin), log(eos_ymumax), 1.e-14, func_y);

  double yp;
  // Now ymu is log(Y_mu)
  // log(Y_mu) <= log(eos_ymumin) is more accurate than Y_mu <= eos_ymumin
  // The floating error with 5e-5 is good to capture correct Yp when T = 0.1 MeV
  if (ymu <= log(eos_ymumin + 0.00001)) {
  //if (ymu <= log(eos_ymumin)) {
    //std::cout << "passing <= eos_ymumin" << rho<<" "<<yle<<" "<<ymu<<" \n";
    double dummy;
    error_type error2;
    // Fix ymu to be lower bound, since log and exp induced floating pt error
    ymu = eos_ymumin;
    yp = yle + ymu;
    yp = std::min(eos_yemax,std::max(yp, eos_yemin));
    press_eps_ye__nu_beta_eq__rho_temp(dummy, yp, rho, temp, mu_nu_e_bar, error2);
  } else {
    ymu = exp(ymu);
  }

  yp = yle + ymu;
  yp = std::min(eos_yemax,std::max(yp, eos_yemin));


  const auto lymu = log(ymu);

  auto const vars =
    alltables_baryon.interpolate<EBARYON::PRESS, EBARYON::CS2, EBARYON::EPS>(lrho, ltemp, yp);

  auto const vars_muon =
    alltables_muon.interpolate<EMUON::PRESS_MU_MINUS, EMUON::PRESS_MU_PLUS,
                               EMUON::EPS_MU_MINUS, EMUON::EPS_MU_PLUS>(lrho, ltemp, lymu);

  double press = 0.0;
  if (add_ele_contribution) {
  	auto const vars_ele  =
  	  alltables_ele.interpolate<EELE::PRESS_E_MINUS, EELE::PRESS_E_PLUS,
  	                            EELE::EPS_E_MINUS, EELE::EPS_E_PLUS>(lrho, ltemp, yle);
  	press = exp(vars[0]) + vars_ele[0] + vars_ele[1] + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[2]) - energy_shift + vars_ele[2] + vars_ele[3] + vars_muon[2] + vars_muon[3];
  } else {
  	press = exp(vars[0]) + vars_muon[0] + vars_muon[1];
  	eps = exp(vars[2]) - energy_shift + vars_muon[2] + vars_muon[3];
  }


  return press;
};

inline double EOS_Leptonic::find_yle_of_mumu_nubetaeq(const double &lrho, const double &ltemp, 
                                                      const double &mu_nu_e_bar, const double &mu_nu_mu_bar, const double &mu_mu){
  double yle_local2;
  auto const func_yle_of_mumu = [&](double &yle_local2) {
   auto const mu_e =
     alltables_ele.interpolate<EELE::MUELE>(lrho, ltemp, yle_local2);
   return mu_e[0] - mu_mu + mu_nu_mu_bar - mu_nu_e_bar;
  };

  yle_local2 = zero_brent(eos_ylemin, eos_ylemax, 1.e-14, func_yle_of_mumu);
  return yle_local2;
}



#include "leptonic_eos_readtable_leptonic.hh"
#include "leptonic_eos_readtable_scollapse.hh"
#include "leptonic_eos_readtable_compose.hh"
#include "leptonic_eos_generate_Cold_table.hh"
#endif
