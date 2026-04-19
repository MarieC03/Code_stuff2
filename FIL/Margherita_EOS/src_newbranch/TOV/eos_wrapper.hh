/*
 * =====================================================================================
 *
 *       Filename:  eos_wrapper.hh
 *
 *    Description:  EOS Wrapper for Margherita EOS
 *
 *        Version:  1.0
 *        Created:  12/01/2018 11:10:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#ifndef TOV_EOS_WRAPPER_H
#define TOV_EOS_WRAPPER_H


template<typename EOS>
class Cold_EOS_Wrapper{

public:
 
  using error_t  = typename EOS::error_t;  

  using EOS_t = EOS;

  static inline double rho__press_cold(double &press, error_t & error){

    return EOS::rho__press_cold(press,error);

  };

  static inline double rho_energy_dedp__press_cold(double &energy, double &dedp, double &press,
      error_t &error){

    auto rho = rho__press_cold(press,error);

    double eps;
    EOS::press_cold_eps_cold__rho(eps,rho,error);

    auto const dpdrho = EOS::dpress_cold_drho__rho(rho,error);

    energy = rho*(1.+eps);
    auto const rhoh = energy + press;

     dedp = rhoh/(dpdrho*rho);

     return rho;

  };

};

#endif
