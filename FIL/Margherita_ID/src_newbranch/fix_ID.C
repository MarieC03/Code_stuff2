/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 * 
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed 
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <limits>
#include <math.h>
#include <sys/time.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "Margherita_EOS.h"
#include "margherita_interface.hh"


template<typename eos>
void fix_ID_type(CCTK_ARGUMENTS){

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  double ymu{} ; // This does not exist in hydrobase and it will be set by ID_converter.. anyway
  ymu = 0.0;
  
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
	int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
	double xL = x[index];
	double yL = y[index];
	double zL = z[index];

	rho[index] = std::max(rho[index],eos::c2p_rho_atm);

        //Atmosphere needs to be treated differently...
        const double tempL = (rho[index] <= 1.001*eos::c2p_rho_atm) ? eos::c2p_temp_atm : eos::eos_tempmin;

        //temperature[index] = ID_is_cold ? tempL : temperature[index];
        //temp[index] = ID_is_cold ? tempL : temperature[index];

          // Harry: it doesnot matter what value of temperature it takes here??? tried Nan or 0.1 as well, same profile outputed
          // or it just gives ye, p, eps from beta eqm with T=0.01mev --> construct ID...
          // becoz it all depends on hotslice.cc in Margherita_EOS
        if (ID_is_cold) {
           temperature[index] = tempL;
        } else {
           temperature[index] = eos::temp_ID;
        }

	if(ID_uses_beta_equilibrium){
		
          typename eos::error_type error; // TODO How do we handle this?


          //CCTK_ERROR("stop first") ;


	  press[index] = Pressure_Interface<eos>::press_eps_yle_ymu__beta_eq__rho_temp(
										       eps[index]
										       , Y_e[index]
										       , ymu
										       , rho[index]
										       , temperature[index]
										       ,error);
	}

	double dummy_csnd2;
	typename eos::error_type error;
	eps[index]=eos::eps_csnd2_entropy__temp_rho_yle_ymu(
							    dummy_csnd2, entropy[index], temperature[index], rho[index], Y_e[index], ymu, error );
	press[index] = eos::press__temp_rho_yle_ymu(temperature[index], rho[index], Y_e[index], ymu, error);


	double ETvx = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
	double ETvy = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
	double ETvz = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];

	if(add_velocity_pert){
   	    const double varpi = r[index]/ r_approx;
	    const double vp = vel_pert/r_approx; //0.5* vel_pert* fabs((SQR(varpi) - 3.) * varpi ); 	   	
	    //Need to decompose this into x,y,z components
	    ETvx += vp*x[index]/std::max(r[index],1.0e-10);
	    ETvy += vp*y[index]/std::max(r[index],1.0e-10);
	    ETvz += vp*z[index]/std::max(r[index],1.0e-10);
	    
	}
	if(rho[index]<=1.000001*eos::c2p_rho_atm){
	       ETvx = 0.;
	       ETvy = 0.;
	       ETvz = 0.;
	    }


	vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = ETvx;
	vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = ETvy;
	vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = ETvz;


	const CCTK_REAL v2= gxx[index]*ETvx*ETvx
			   +gyy[index]*ETvy*ETvy
			   +gzz[index]*ETvz*ETvz
			   +2.0*(
				gxy[index]*ETvx*ETvy +
				gxz[index]*ETvx*ETvz +
				gyz[index]*ETvz*ETvy 
				);
	w_lorentz[index]= 1./sqrt(1.-v2);
	assert(std::isfinite(w_lorentz[index]));

      }
}

#define VERR_DEF_PARAMS __LINE__, __FILE__, CCTK_THORNSTRING
extern "C" void fix_ID(CCTK_ARGUMENTS){
  DECLARE_CCTK_PARAMETERS;

 if(CCTK_Equals(eos_type,"tabulated"))
  fix_ID_type<EOS_Tabulated>(CCTK_PASS_CTOC);  
 else if(CCTK_Equals(eos_type,"tabulated_leptonic"))
  fix_ID_type<EOS_Leptonic>(CCTK_PASS_CTOC);
 else if(CCTK_Equals(eos_type,"polytropic"))
  fix_ID_type<EOS_Polytropic>(CCTK_PASS_CTOC);  
 else if(CCTK_Equals(eos_type,"tabulated_cold"))
  fix_ID_type<EOS_Hybrid_ColdTable>(CCTK_PASS_CTOC);  
 else
    CCTK_VError(VERR_DEF_PARAMS,"EOS unknown..");
}
