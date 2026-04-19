/*
 * =====================================================================================
 *
 *       Filename:  leak_constants.hh
 *
 *    Description:  Constants for Leakage
 *
 *        Version:  1.0
 *        Created:  13/05/2017 21:17:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#ifndef _HH_LEAK_CONSTANTS
#define _HH_LEAK_CONSTANTS

namespace Leak_Constants {
constexpr double me_mev = 0.510998910;  // mass of electron in MeV
constexpr double sigma_0 = 1.76e-44;    // in units of cm^2
constexpr double alpha = 1.23;          // dimensionless
constexpr double Qnp = 1.293333;        // neutron proton mass difference in MeV
constexpr double hc_mevcm = 1.23984172e-10;  // hc in units of MeV*cm
//  constexpr double rho_min_leak = 1.0e4 //in g/cm^3
constexpr double Cv = 0.5 + 2.0 * 0.23;  //  dimensionless
constexpr double Ca = 0.5;               // dimensionless
constexpr double gamma_0 = 5.565e-2;     // dimensionless
constexpr double fsc = 1.0 / 137.036;  // fine structure constant, dimensionless
constexpr double pi = 3.14159265358979;
constexpr double clite = 2.99792458e10;
constexpr double a_diff = 6.;  // O'Connor says 6., Sekiguchi 3. ...

constexpr double mev_to_erg = 1.60217733e-6;
constexpr double erg_to_mev = 6.24150636e5;
constexpr double amu_cgs = 1.66053873e-24;
constexpr double massn_cgs = 1.674927211e-24;
constexpr double amu_mev = 931.49432e0;
constexpr double kb_erg = 1.380658e-16;
constexpr double kb_mev = 8.61738568e-11;
constexpr double temp_mev_to_kelvin = 1.1604447522806e10;
constexpr double planck = 6.626176e-27;
constexpr double avo = 6.0221367e23;
constexpr double hbarc_mevcm = 1.97326966e-11;

constexpr double beta = pi * clite * (1.0 + 3.0 * alpha * alpha) * sigma_0 /
                        (hc_mevcm * hc_mevcm * hc_mevcm * me_mev * me_mev);

constexpr double normfact = 1.;

}  // namespace Leak_Constants

#endif
