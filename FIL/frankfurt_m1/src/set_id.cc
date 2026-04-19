#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "frankfurt_m1.h"
#include "Margherita_M1.h"

#define VERR_DEF_PARAMS __LINE__, __FILE__, CCTK_THORNSTRING

extern "C" void m1_set_initial(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    using namespace Margherita_M1_EAS ;


    //CCTK_VInfo(CCTK_THORNSTRING,"Setting up initial data for M1 standard test: %s", M1_WhichTest);

    CCTK_REAL * restrict velx = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,0)]);
    CCTK_REAL * restrict vely = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,1)]);
    CCTK_REAL * restrict velz = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,2)]);

    #pragma omp parallel for
    for(int ijk=0;ijk<cctk_lsh[2]*cctk_lsh[1]*cctk_lsh[0];++ijk) {
        //Opacity
        kappa_nue_a[ijk] = 0.;
        kappa_nue_bar_a[ijk] = 0.;
        kappa_numu_a[ijk] = 0.;
        kappa_numu_bar_a[ijk] = 0.;
        kappa_nux_a[ijk] = 0.;

        Qnue[ijk] = kappa_nue_a[ijk];
        Qnue_bar[ijk] = kappa_nue_bar_a[ijk];
        Qnumu[ijk] = kappa_numu_a[ijk];
        Qnumu_bar[ijk] = kappa_numu_bar_a[ijk];
        Qnux[ijk] = kappa_nux_a[ijk];
  
        Rnue[ijk] = 0.;
        Rnue_bar[ijk] = 0.;
        Rnumu[ijk] = 0.;
        Rnumu_bar[ijk] = 0.;
        Rnux[ijk] = 0.;
  
        eta_nue[ijk] = 0.;
        eta_nue_bar[ijk] = 0.;
        eta_numu[ijk] = 0.;
        eta_numu_bar[ijk] = 0.;
        eta_nux[ijk] = 0.;

        kappa_nue_s[ijk] = 0.;
        kappa_nue_bar_s[ijk] = 0.;
        kappa_numu_s[ijk] = 0.;
        kappa_numu_bar_s[ijk] = 0.;
        kappa_nux_s[ijk] = 0.;

        eps_nue[ijk] = 0.;
        eps_nue_bar[ijk] = 0.;
        eps_numu[ijk] = 0.;
        eps_numu_bar[ijk] = 0.;
        eps_nux[ijk] = 0.;

        //Radiation
        CCTK_REAL const Eatm = E_atmo ;
        CCTK_REAL const Natm = N_atmo ;

        Enue[ijk] = Eatm;
        Enue_bar[ijk] = Eatm;
        Enumu[ijk] = Eatm;
        Enumu_bar[ijk] = Eatm;
        Enux[ijk] = Eatm;

        Nnue[ijk] = Natm;
        Nnue_bar[ijk] = Natm;
        Nnumu[ijk] = Natm;
        Nnumu_bar[ijk] = Natm;
        Nnux[ijk] = Natm;

        Fnue_x[ijk] = 0.; Fnue_y[ijk]=0.; Fnue_z[ijk]=0.;
        Fnue_bar_x[ijk] = 0.; Fnue_bar_y[ijk]=0.; Fnue_bar_z[ijk]=0.;
        Fnumu_x[ijk] = 0.; Fnumu_y[ijk]=0.; Fnumu_z[ijk]=0.;
        Fnumu_bar_x[ijk] = 0.; Fnumu_bar_y[ijk]=0.; Fnumu_bar_z[ijk]=0.;
        Fnux_x[ijk] = 0.; Fnux_y[ijk]=0.; Fnux_z[ijk]=0.;

//CCTK_ERROR("Code test, passing m1_set_initial");

        //const auto tau_init = EAS<double>::compute_analytic_opacity(rho_b[ijk]);
//    CCTK_INFO("passing ser tau = 0") ;
//
//	tau_nue[ijk] = 0.;
//    CCTK_INFO("passed ser tau = 0") ;

    }
//CCTK_ERROR("Code test, passed m1 initialization");
} //




