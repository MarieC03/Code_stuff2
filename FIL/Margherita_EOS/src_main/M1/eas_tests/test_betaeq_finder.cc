#include <array>
#include <iostream>
#include <iomanip>

#define DEBUG
#define STANDALONE
//#define _NO_GSL__
#define DEBUG_ROOTFINDING
#define COLDTABLE_SETUP
//#define DEBUG
#define DEBUG_BETAEQ

using CCTK_REAL = double;
using CCTK_INT = int ;

#include "../../3D_Table/readtable.cc"
#include "../../Margherita_M1.h"

const std::string table_name = "/home/relastro-shared/EOS/scollapse/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5" ; //TODO

int main() {
    using namespace Margherita_M1_EAS; 
    constexpr bool use_energy_shift = true ;
    constexpr bool recompute_mu_nu = false ;
    EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, recompute_mu_nu);

    std::cout << "Read Table" << std::endl;


    double rho = 7.582629253507815e14*Margherita_constants::RHOGF;//7.582629253507815e14*Margherita_constants::RHOGF;
    double ye = 0.234693877551020;
    double temp = 16.852590447507527;

    //rho = 0.00111869;
    //temp = 0.01;//29.63349846;
    //ye = 0.1;

    CCTK_REAL _T,_ye ;

    std::array<CCTK_REAL,REQ_FLUID_VARS> U;
    U[RHO]  = rho;
    U[YE]   = ye;
    U[TEMP] = temp;

    std::cout << "rho: " << rho << std::endl;
    std::cout << "ye: " << ye << std::endl;
    std::cout << "temp: " << temp << std::endl;

    const auto tau_init= EAS<double>::compute_analytic_opacity(rho);

    std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
    auto tau_e = tau_n;

    auto F= Fugacities<double>(rho,temp,ye, std::move(tau_n));

    std::array<bool, NINTERACTIONS> which_interactions { true, true, true, true } ;
    auto eas = EAS<double>(F, which_interactions);
    eas.calc_eas() ; 

}