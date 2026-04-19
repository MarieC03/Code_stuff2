/**
 * @file driver_evaluate_m1_rhs.cc
 * 
 * Compute explicit part of M1 evolution
 * i.e. fluxes and geometric sources 
 *
 * Everything is computed one species 
 * of neutrinos at a time
 *
 * The branch optimize contains a version
 * that does all species at a time,
 * it's just as fast (slow?).
 * This version is cleaner in my
 * opinion
 *
 * The structure follows the one from 
 * the Frankfurt-IllinoisGRMHD code 
 * by Elias Most and that of the 
 * original M1 scheme by Lukas Weih and
 * E. Most.
 *
 * @author: Carlo Musolino
 */

#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <tuple>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "frankfurt_m1.h" /* Generic #define's and function prototypes */
#include "driver_evaluate_M1_rhs.h" /* Function prototypes for this file only */
#include "m1_closure.hh"

#define velx_p (&vel_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p (&vel_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p (&vel_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

extern "C" void frankfurt_m1_calc_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  using namespace frankfurt_m1 ;

  if( ! evolve_radiation )
    return ;

  int const substep = 2 - (*RK_counter);

  using recon_function_t = void (*) (const cGH*, const int*, const int, const int, const int*,
				    gf_and_gz_struct *, gf_and_gz_struct *, gf_and_gz_struct *, CCTK_REAL* ) ;


  recon_function_t reconstruct_set_of_prims ;
  
  if ( CCTK_Equals(m1_recon_type,"minmod" ) ) {
    reconstruct_set_of_prims = (&reconstruct_set_of_prims_MINMOD) ;
  } else if ( CCTK_Equals(m1_recon_type,"mc2") ) {
    reconstruct_set_of_prims = (&reconstruct_set_of_prims_MC2) ;
  } else if ( CCTK_Equals(m1_recon_type,"weno5") ) {
    reconstruct_set_of_prims = (&reconstruct_set_of_prims_WENO5)  ;
  } 
  
  int species = (*frankfurt_m1_rhs_species) - 1; 

  if( verbosity>1)
    CCTK_VInfo(CCTK_THORNSTRING, "RHS species index: %d",species ) ; 
  
  CCTK_REAL dX[3] = { CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2) };

  // in_prims,out_prims_r, and out_prims_l are arrays of pointers to the actual gridfunctions.
  // for gf_and_gz_struct definition see frankfurt_m1.h -- or IllinoisGRMHD_headers.h
  gf_and_gz_struct in_prims[NUM_RECON_PRIM],out_prims_r[NUM_RECON_PRIM],out_prims_l[NUM_RECON_PRIM];
  int which_prims_to_reconstruct[NUM_RECON_PRIM],num_prims_to_reconstruct;


  // We read these in just because we need their reconstructed values
  #pragma omp simd
  for(int i=0;i<cctk_lsh[0];i++)for(int j=0;j<cctk_lsh[1];j++)for(int k=0;k<cctk_lsh[2];k++){
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /*
    v_x[index] = zvecx[index];
    v_y[index] = zvecy[index];
    v_z[index] = zvecz[index];
    */
	  v_x[index] = (substep==0) * velx_p[index] * w_lorentz_p[index] + (substep==1) * zvecx[index];
	  v_y[index] = (substep==0) * vely_p[index] * w_lorentz_p[index] + (substep==1) * zvecy[index];
	  v_z[index] = (substep==0) * velz_p[index] * w_lorentz_p[index] + (substep==1) * zvecz[index];
  }

  // This routine deals with just one neutrino species
  // the pointers below will be set to the relevant GFs 
  CCTK_REAL *Enu, *Fnux, *Fnuy, *Fnuz, *Nnu ;
  CCTK_REAL *k_a, *k_s;
  CCTK_REAL *zeta ; 
  CCTK_REAL *Enu_rhs, *Fx_rhs, *Fy_rhs, *Fz_rhs, *Nnu_rhs;
  CCTK_REAL *Pnu[6] ;
  
  // Set pointers to correct species 
  if( species == 0 ) {
    zeta  = zeta_e   ;
    Enu = Enue; Fnux=Fnue_x; Fnuy = Fnue_y; Fnuz=Fnue_z; Nnu=Nnue;
    Enu_rhs = Enue_rhs; Fx_rhs = Fnue_x_rhs; Fy_rhs=Fnue_y_rhs; Fz_rhs=Fnue_z_rhs ;
    Nnu_rhs = Nnue_rhs ;
    k_a=kappa_nue_a; k_s=kappa_nue_s;
    Pnu[0] = PUPe_xx;
    Pnu[1] = PUPe_xy;
    Pnu[2] = PUPe_xz;
    Pnu[3] = PUPe_yy;
    Pnu[4] = PUPe_yz;
    Pnu[5] = PUPe_zz;    
    //Electron neutrinos
  } else if ( species==1 ) {
    zeta  = zeta_a ;
    Enu = Enua; Fnux=Fnua_x; Fnuy = Fnua_y; Fnuz=Fnua_z; Nnu=Nnua;
    Enu_rhs = Enua_rhs; Fx_rhs = Fnua_x_rhs; Fy_rhs=Fnua_y_rhs; Fz_rhs=Fnua_z_rhs ;
    Nnu_rhs = Nnua_rhs ;
    k_a=kappa_nua_a; k_s=kappa_nua_s;
    Pnu[0] = PUPa_xx;
    Pnu[1] = PUPa_xy;
    Pnu[2] = PUPa_xz;
    Pnu[3] = PUPa_yy;
    Pnu[4] = PUPa_yz;
    Pnu[5] = PUPa_zz;
    // Electron antineutrinos
  } else if ( species == 2 ) {
    zeta  = zeta_x   ;
    Enu = Enux; Fnux=Fnux_x; Fnuy = Fnux_y; Fnuz=Fnux_z; Nnu=Nnux;
    Enu_rhs = Enux_rhs; Fx_rhs = Fnux_x_rhs; Fy_rhs=Fnux_y_rhs; Fz_rhs=Fnux_z_rhs ;
    Nnu_rhs = Nnux_rhs ;
    k_a=kappa_nux_a; k_s=kappa_nux_s;
    Pnu[0] = PUPx_xx;
    Pnu[1] = PUPx_xy;
    Pnu[2] = PUPx_xz;
    Pnu[3] = PUPx_yy;
    Pnu[4] = PUPx_yz;
    Pnu[5] = PUPx_zz;
    // Effective heavy neutrino
  } else {
    CCTK_ERROR("Code should not be here!") ;
  }

  int ww=0;
  in_prims[ww].gf=Nnu; out_prims_r[ww].gf=Nnu_r; out_prims_l[ww].gf=Nnu_l; ww++;
  in_prims[ww].gf=Enu; out_prims_r[ww].gf=Enu_r; out_prims_l[ww].gf=Enu_l; ww++;
  in_prims[ww].gf=Fnux; out_prims_r[ww].gf=Fnux_r; out_prims_l[ww].gf=Fnux_l; ww++;
  in_prims[ww].gf=Fnuy; out_prims_r[ww].gf=Fnuy_r; out_prims_l[ww].gf=Fnuy_l; ww++;
  in_prims[ww].gf=Fnuz; out_prims_r[ww].gf=Fnuz_r; out_prims_l[ww].gf=Fnuz_l; ww++;
  
  in_prims[ww].gf=v_x; out_prims_r[ww].gf=v_x_r; out_prims_l[ww].gf=v_x_l; ww++;
  in_prims[ww].gf=v_y; out_prims_r[ww].gf=v_y_r; out_prims_l[ww].gf=v_y_l; ww++;
  in_prims[ww].gf=v_z; out_prims_r[ww].gf=v_z_r; out_prims_l[ww].gf=v_z_l; ww++;

  // these won't be reconstructed -- they are just used as initial
  // guess for the NR in the closure solver
  in_prims[ww].gf=zeta; out_prims_r[ww].gf=nullptr; out_prims_l[ww].gf=nullptr; ww++;

  int numvars = ww ;
  
  // Prims are defined AT ALL GRIDPOINTS, so we set the # of ghostzones to zero:
  for(int i=0;i<numvars;i++) for(int j=1;j<=3;j++) { in_prims[i].gz_lo[j]=0; in_prims[i].gz_hi[j]=0; }
  // Left/right variables are not yet defined, yet we set the # of gz's to zero by default:
  for(int i=0;i<numvars;i++) for(int j=1;j<=3;j++) { out_prims_r[i].gz_lo[j]=0; out_prims_r[i].gz_hi[j]=0; }
  for(int i=0;i<numvars;i++) for(int j=1;j<=3;j++) { out_prims_l[i].gz_lo[j]=0; out_prims_l[i].gz_hi[j]=0; }
  
  /* SET POINTERS TO METRIC GRIDFUNCTIONS */
  CCTK_REAL *metric[NUMVARS_METRIC]; // "metric" here is array of pointers to the actual gridfunctions.
  CCTK_REAL *k_xx, *k_xy, *k_xz, *k_yy, *k_yz, *k_zz ;
  if( substep == 0 ) {
    ww=0;
    metric[ww]=gxx_p;    ww++;
    metric[ww]=gxy_p;    ww++;
    metric[ww]=gxz_p;    ww++;
    metric[ww]=gyy_p;    ww++;
    metric[ww]=gyz_p;    ww++;
    metric[ww]=gzz_p;    ww++;
    metric[ww]=alp_p;    ww++;
    metric[ww]=betax_p;   ww++;
    metric[ww]=betay_p;   ww++;
    metric[ww]=betaz_p;   ww++;

    k_xx = kxx_p ; k_xy = kxy_p ;
    k_xz = kxz_p ; k_yy = kyy_p ;
    k_yz = kyz_p ; k_zz = kzz_p ;
    
  } else {
    ww=0;
    metric[ww]=gxx;    ww++;
    metric[ww]=gxy;    ww++;
    metric[ww]=gxz;    ww++;
    metric[ww]=gyy;    ww++;
    metric[ww]=gyz;    ww++;
    metric[ww]=gzz;    ww++;
    metric[ww]=alp;    ww++;
    metric[ww]=betax;   ww++;
    metric[ww]=betay;   ww++;
    metric[ww]=betaz;   ww++;

    k_xx = kxx ; k_xy = kxy ;
    k_xz = kxz ; k_yy = kyy ;
    k_yz = kyz ; k_zz = kzz ;
    
  }
  
  
  // Add K_(ij)P^(ij) to E_rhs and store P^(ij) in its GF 
  compute_E_rhs_extrinsic_curvature_terms_and_Pmunu(cctkGH,cctk_lsh,cctk_nghostzones,metric,
						    Pnu,
						    in_prims,
						    k_xx,k_xy,k_xz,k_yy,k_yz,k_zz,
						    Enu_rhs, zeta) ;

  // reconstruction is performed on F/E N/E
  // See Foucart+ but the reason is that this
  // helps avoiding unphysical fluxes where F_sq > E_sq
  #pragma omp parallel for simd
 for(int ijk=0;ijk<cctk_lsh[2]*cctk_lsh[1]*cctk_lsh[0];++ijk) {
    Fnux[ijk] /= Enu[ijk];
    Fnuy[ijk] /= Enu[ijk];
    Fnuz[ijk] /= Enu[ijk];
    Nnu[ijk] /= Enu[ijk];
  }


  ww=0;
  which_prims_to_reconstruct[ww]=NNU;          ww++;
  which_prims_to_reconstruct[ww]=ENU;	       ww++;
  which_prims_to_reconstruct[ww]=FNUX;	       ww++;
  which_prims_to_reconstruct[ww]=FNUY;	       ww++;
  which_prims_to_reconstruct[ww]=FNUZ;	       ww++;

  which_prims_to_reconstruct[ww]=VELX;	       ww++;
  which_prims_to_reconstruct[ww]=VELY;	       ww++;
  which_prims_to_reconstruct[ww]=VELZ;	       ww++;
  
  num_prims_to_reconstruct=ww;

  /* Once per direction:
   * 1) Reconstruct
   * 2) Compute fluxes and ( vector ) source terms
   *    and add to RHSs
   */
  
  int flux_dirn ;
  flux_dirn = 1 ; // X direction
  // This function is housed in the file: "reconstruct_set_of_prims_MC2.C"
  reconstruct_set_of_prims(cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
				  in_prims,out_prims_r,out_prims_l,temporary);
  // Then add fluxes to RHS for M1 variables:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_M1_rhss.C"
  add_fluxes_and_source_terms_to_M1_rhss<0>(cctkGH,cctk_lsh,cctk_nghostzones,dX,metric,in_prims,
					    num_prims_to_reconstruct,out_prims_r,out_prims_l,
					    k_a, k_s,
					    Pnu,
					    Nflux, Eflux, Fxflux, Fyflux, Fzflux,
					    Nnu_rhs, Enu_rhs, Fx_rhs, Fy_rhs, Fz_rhs,
					    Nflux_LF, Eflux_LF, Fxflux_LF, Fyflux_LF, Fzflux_LF);
  
  flux_dirn=2; // Y direction
  // This function is housed in the file: "reconstruct_set_of_prims_M.C"
  reconstruct_set_of_prims(cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
				  in_prims,out_prims_r,out_prims_l,temporary);
  // Then add fluxes to RHS for M1 variables:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_M1_rhss.C"
  add_fluxes_and_source_terms_to_M1_rhss<1>(cctkGH,cctk_lsh,cctk_nghostzones,dX,metric,in_prims,
					    num_prims_to_reconstruct,out_prims_r,out_prims_l,
					    k_a, k_s,
					    Pnu,
					    Nflux, Eflux, Fxflux, Fyflux, Fzflux,
					    Nnu_rhs, Enu_rhs, Fx_rhs, Fy_rhs, Fz_rhs,
					    Nflux_LF, Eflux_LF, Fxflux_LF, Fyflux_LF, Fzflux_LF);
  
  flux_dirn=3; // Z direction
  // This function is housed in the file: "reconstruct_set_of_prims_MC2.C"
  reconstruct_set_of_prims(cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
				  in_prims,out_prims_r,out_prims_l,temporary);
  // Then add fluxes to RHS for M1 variables:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_M1_rhss.C"
  add_fluxes_and_source_terms_to_M1_rhss<2>(cctkGH,cctk_lsh,cctk_nghostzones,dX,metric,in_prims,
					    num_prims_to_reconstruct,out_prims_r,out_prims_l,
					    k_a, k_s,
					    Pnu,
					    Nflux, Eflux, Fxflux, Fyflux, Fzflux,
					    Nnu_rhs, Enu_rhs, Fx_rhs, Fy_rhs, Fz_rhs,
					    Nflux_LF, Eflux_LF, Fxflux_LF, Fyflux_LF, Fzflux_LF);
  // Remove normalisation 
#pragma omp parallel for simd
  for(int ijk=0;ijk<cctk_lsh[2]*cctk_lsh[1]*cctk_lsh[0];++ijk) {
    Fnux[ijk] *= Enu[ijk];
    Fnuy[ijk] *= Enu[ijk];
    Fnuz[ijk] *= Enu[ijk];
    Nnu[ijk]  *= Enu[ijk];
  }
  
  return;
}


// We add #include's here instead of compiling these separately to help ensure that functions are properly inlined.
//    These files only include about 800 lines of code in total (~1200 lines in total), but it's arguably more
//    convenient to edit a 600 line file than an 1800 line file, so I'd prefer to leave this unconventional structure
//    alone.
#include "reconstruct_set_of_prims_MINMOD.cc"
#include "reconstruct_set_of_prims_MC2.cc"
#include "reconstruct_set_of_prims_WENO5.cc" 
#include "compute_E_rhs_extrinsic_curvature_terms_and_Pmunu.cc"
#include "add_fluxes_and_source_terms_to_M1_rhss.cc"
