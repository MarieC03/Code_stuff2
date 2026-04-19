/**
 * @file add_fluxes_and_source_terms_to_M1_rhss.cc
 *
 * This file is part of frankfurt_m1
 *
 *  Compute fluxes and RHSs of M1 equations
 *
 *  Based on previous version by L. Weih and E. R. Most.
 *
 * @author  Carlo Musolino
 */ 

// Side note: the following values could be used for cell averaged gfs: 
//     am2=-1.0/12.0, am1=7.0/12.0, a0=7.0/12.0, a1=-1.0/12.0
// However, since the metric gfs store the grid point values instead of the cell average, 
//     the following coefficients should be used: 
//     am2 = -1/16, am1 = 9/16, a0 = 9/16, a1 = -1/16
// This will yield the third-order-accurate face values at m-1/2, 
//      using values specified at {m-2,m-1,m,m+1}
#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625

#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

#define COMPUTE_DERIV_4(METRICm2,METRICm1,METRICp1,METRICp2) ((1./12.)*(METRICm2) + (-8./12.)*(METRICm1) + (8./12.)*(METRICp1) + (-1./12.)*(METRICp2))

#include "ECHO_helper_functions.hh"


#define COMPUTE_FOURMETRIC(g4tt,g4tx,g4ty,g4tz,g4xx,g4xy,g4xz,g4yy,g4yz,g4zz,METRIC,METRIC_AUX)  ( { \
      /* g_{0i} = beta_i */                                             \
      g4tx = METRIC_AUX[PSI4]*(METRIC[GXX]*METRIC[SHIFTX] + METRIC[GXY]*METRIC[SHIFTY] + METRIC[GXZ]*METRIC[SHIFTZ]); \
      g4ty = METRIC_AUX[PSI4]*(METRIC[GXY]*METRIC[SHIFTX] + METRIC[GYY]*METRIC[SHIFTY] + METRIC[GYZ]*METRIC[SHIFTZ]); \
      g4tz = METRIC_AUX[PSI4]*(METRIC[GXZ]*METRIC[SHIFTX] + METRIC[GYZ]*METRIC[SHIFTY] + METRIC[GZZ]*METRIC[SHIFTZ]); \
      /* g_{00} = -alpha^2 + beta^i beta^j gamma_{ij} = -alpha^2 + beta^i beta_i = -alpha^2 + beta^i g_{0i} */ \
      g4tt = -SQR(METRIC_AUX[LAPSE]) + g4tx*METRIC[SHIFTX] + g4ty*METRIC[SHIFTY] + g4tz*METRIC[SHIFTZ]; \
      g4xx = METRIC_AUX[PSI4]*METRIC[GXX];                              \
      g4xy = METRIC_AUX[PSI4]*METRIC[GXY];                              \
      g4xz = METRIC_AUX[PSI4]*METRIC[GXZ];                              \
      g4yy = METRIC_AUX[PSI4]*METRIC[GYY];                              \
      g4yz = METRIC_AUX[PSI4]*METRIC[GYZ];                              \
      g4zz = METRIC_AUX[PSI4]*METRIC[GZZ];                              \
    } )


// Harry: this is reading one of the nu species
template <size_t flux_dirn>
static void add_fluxes_and_source_terms_to_M1_rhss(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,CCTK_REAL *dX,
						   CCTK_REAL **metric,gf_and_gz_struct *in_prims,
						   int numvars_reconstructed,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l,
						   CCTK_REAL *kappa_nu_a, CCTK_REAL* kappa_nu_s,
						   CCTK_REAL **Pnu,
						   CCTK_REAL *Nnu_flux,
						   CCTK_REAL *Enu_flux, CCTK_REAL *Fnu_x_flux, CCTK_REAL *Fnu_y_flux, CCTK_REAL *Fnu_z_flux,
						   CCTK_REAL *Nnu_rhs, 
						   CCTK_REAL *Enu_rhs, CCTK_REAL *Fnu_x_rhs, CCTK_REAL *Fnu_y_rhs, CCTK_REAL *Fnu_z_rhs,
						   CCTK_REAL *Nnu_flux_LF ,
						   CCTK_REAL *Enu_flux_LF, CCTK_REAL *Fnu_x_flux_LF, CCTK_REAL *Fnu_y_flux_LF, CCTK_REAL *Fnu_z_flux_LF)  {
  using namespace frankfurt_m1 ;

  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL dxi[3] = { 1.0/dX[0],1.0/dX[1],1.0/dX[2] };

  // Notice in the loop below that we go from 3 to cctk_lsh-2 for i, j, AND k, even though
  //   we are only computing the flux in one direction at a time. This is because in the end,
  //   we only need the rhs's from 3 to cctk_lsh-3 for i, j, and k.

#pragma omp parallel for collapse(3)
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2] - cctk_nghostzones[2]+1;k++) 
    for(int j=cctk_nghostzones[1];j<cctk_lsh[1]- cctk_nghostzones[1]+1;j++) 
      for(int i=cctk_nghostzones[0];i<cctk_lsh[0] - cctk_nghostzones[0]+1;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	// Set metric and associated variables
	CCTK_REAL METRIC[NUMVARS_METRIC]; for(int ii=0;ii<NUMVARS_METRIC;ii++) METRIC[ii] = metric[ii][index];	
	// Next set metric on the faces, applying a 3rd-order lopsided stencil.
	int indexm2 = CCTK_GFINDEX3D(cctkGH,i-2*kronecker_delta(flux_dirn,0),j-2*kronecker_delta(flux_dirn,1),k-2*kronecker_delta(flux_dirn,2));
	int indexm1 = CCTK_GFINDEX3D(cctkGH,i-  kronecker_delta(flux_dirn,0),j-  kronecker_delta(flux_dirn,1),k-  kronecker_delta(flux_dirn,2));
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+  kronecker_delta(flux_dirn,0),j+  kronecker_delta(flux_dirn,1),k+  kronecker_delta(flux_dirn,2));
	int indexp2 = CCTK_GFINDEX3D(cctkGH,i+2*kronecker_delta(flux_dirn,0),j+2*kronecker_delta(flux_dirn,1),k+2*kronecker_delta(flux_dirn,2));
	// The "vector" METRIC stores needed metric-related quantities.
	CCTK_REAL METRICm2[NUMVARS_METRIC]; for(int ii=0;ii<NUMVARS_METRIC;ii++) METRICm2[ii] = metric[ii][indexm2];
	CCTK_REAL METRICm1[NUMVARS_METRIC]; for(int ii=0;ii<NUMVARS_METRIC;ii++) METRICm1[ii] = metric[ii][indexm1];
	CCTK_REAL METRICp1[NUMVARS_METRIC]; for(int ii=0;ii<NUMVARS_METRIC;ii++) METRICp1[ii] = metric[ii][indexp1];
	CCTK_REAL METRICp2[NUMVARS_METRIC]; for(int ii=0;ii<NUMVARS_METRIC;ii++) METRICp2[ii] = metric[ii][indexp2];
	
	// Next compute the metric values at the {i,j,k} +/- 1/2 faces (i.e., the "face values" of the metric)
	CCTK_REAL FACEVAL[NUMVARS_METRIC],FACEVALp1[NUMVARS_METRIC];
	for(int w=0;w<NUMVARS_METRIC;w++) FACEVAL[w]   = COMPUTE_FCVAL(METRICm2[w],METRICm1[w],METRIC[w],METRICp1[w]);
	for(int w=0;w<NUMVARS_METRIC;w++) FACEVALp1[w] = COMPUTE_FCVAL(METRICm1[w],METRIC[w],METRICp1[w],METRICp2[w]);
	
	// initialise metric object at cell interface i - 1/2 
	metric_c gamma_face( {FACEVAL[GXX], FACEVAL[GXY], FACEVAL[GXZ], FACEVAL[GYY], FACEVAL[GYZ], FACEVAL[GZZ]},
			  {FACEVAL[SHIFTX], FACEVAL[SHIFTY], FACEVAL[SHIFTZ]}, FACEVAL[LAPSE]) ;
	// Initial guess for NR solver 
	auto const zeta_old = 0.5 * ( in_prims[ZETA].gf[index] + in_prims[ZETA].gf[indexm1] ) ;
	// Read reconstructed primitives at i-1/2 +- eps 
	auto El  = out_prims_l[ENU].gf[index] ;
	auto Nl  = out_prims_l[NNU].gf[index]  * El ;
	auto Fxl = out_prims_l[FNUX].gf[index] * El ;
	auto Fyl = out_prims_l[FNUY].gf[index] * El ;
	auto Fzl = out_prims_l[FNUZ].gf[index] * El ;
	
	auto Er  = out_prims_r[ENU].gf[index] ;
	auto Nr  = out_prims_r[NNU].gf[index]  * Er ;
	auto Fxr = out_prims_r[FNUX].gf[index] * Er ;
	auto Fyr = out_prims_r[FNUY].gf[index] * Er ;
	auto Fzr = out_prims_r[FNUZ].gf[index] * Er ;

	// Initialise closure for left and right face values
	closure_t cl( {out_prims_l[VELX].gf[index], out_prims_l[VELY].gf[index], out_prims_l[VELZ].gf[index]}, &gamma_face ) ;
	closure_t cr( {out_prims_r[VELX].gf[index], out_prims_r[VELY].gf[index], out_prims_r[VELZ].gf[index]}, &gamma_face ) ;

	// compute closures	
	cl.update_closure<m1_closure_f>( El,
					 { Fxl,
					   Fyl,
					   Fzl }, zeta_old, true ) ;
	
	
	cr.update_closure<m1_closure_f>( Er,
					 { Fxr,
					   Fyr,
					   Fzr }, zeta_old, true ) ;
       
	

	//We need the factor for the flux-corrections
	auto const kappa_bar = std::sqrt(  ( kappa_nu_a[index] + kappa_nu_s[index] ) * ( kappa_nu_a[indexm1] + kappa_nu_s[indexm1])) + 1.e-45 ; 
	auto A  = std::min(1., 1./(kappa_bar) * dxi[flux_dirn]);

	
	auto const sqrtg = gamma_face.sqrtdet;
	
	// get wavespeeds
	double cp_l, cm_l, cp_r, cm_r ;

	cl.get_wavespeeds<flux_dirn>(cp_l,cm_l) ;
	cr.get_wavespeeds<flux_dirn>(cp_r,cm_r) ;

	double cmax = std::max(std::max(cp_l,cp_r), 0. ) ;
	double cmin = -std::min(std::min(cm_l,cm_r), 0. ) ;

	if ( cmax < speed_eps && cmin < speed_eps) {
	  cmax = 1.;
	  cmin = 1.;
	}

	// We need to correct the fluxes to recover the correct diffusion
	// limit in the asymptotically thick regime.
	// Here we follow Skinner+ ( https://collaborate.princeton.edu/en/publications/fornax-a-flexible-code-for-multiphysics-astrophysical-simulations)
	// And Kuroda+ (https://iopscience.iop.org/article/10.3847/0067-0049/222/2/20)
	auto const hlle_E_flux = [&] (double const& R, double const& L, double const& FR, double const& FL) {
				     return ( cmax*FL + cmin*FR - A*cmax*cmin*(R-L) ) / ( cmax + cmin ) ;
				 };

	auto const hlle_F_flux = [&] (double const& R, double const& L, double const& FR, double const& FL) {
				     return ( A*A*(cmax*FL + cmin*FR) - A*cmax*cmin*(R-L) ) / ( cmax + cmin ) + 0.5*(1.-A*A)*(FR+FL) ;
				 };

	// Compute left and right fluxes 
	double Nfl, Efl, Fxfl, Fyfl, Fzfl;
	cl.get_fluxes<flux_dirn>(Nl, Nfl, Efl, Fxfl, Fyfl, Fzfl ) ;
	
	double Nfr, Efr, Fxfr, Fyfr, Fzfr;
	cr.get_fluxes<flux_dirn>(Nr, Nfr, Efr, Fxfr, Fyfr, Fzfr ) ;
	
	// Solve riemann problem approximately 
	Enu_flux[index] = hlle_E_flux( Er*sqrtg, El*sqrtg, Efr, Efl ) ;
	Nnu_flux[index] = hlle_E_flux( Nr*sqrtg*cr.Gamma, Nl*sqrtg*cl.Gamma, Nfr, Nfl ) ;

	Fnu_x_flux[index] = hlle_F_flux(Fxr*sqrtg, Fxl*sqrtg, Fxfr, Fxfl ) ;
	Fnu_y_flux[index] = hlle_F_flux(Fyr*sqrtg, Fyl*sqrtg, Fyfr, Fyfl ) ;
	Fnu_z_flux[index] = hlle_F_flux(Fzr*sqrtg, Fzl*sqrtg, Fzfr, Fzfl ) ;
	
	#ifdef FRANKFURT_M1_USE_PPLIM
	// Mix in some diffusive fluxes 
	// to ensure positivity
	// Here the left -- right
	// states are the cell centered
	// values at i and i-1
	// the diffusive flux is thus computed
	// according to a Godunov like LF prescription
	Er  = in_prims[ENU].gf[index] ;
	Nr  = in_prims[NNU].gf[index]  * Er ;
	Fxr = in_prims[FNUX].gf[index] * Er ;
	Fyr = in_prims[FNUY].gf[index] * Er ;
	Fzr = in_prims[FNUZ].gf[index] * Er ;

	El  = in_prims[ENU].gf[indexm1] ;
	Nl  = in_prims[NNU].gf[indexm1]  * El ;
	Fxl = in_prims[FNUX].gf[indexm1] * El ;
	Fyl = in_prims[FNUY].gf[indexm1] * El ;
	Fzl = in_prims[FNUZ].gf[indexm1] * El ;

	metric_c gamma_cR( { METRIC[GXX],
			     METRIC[GXY],
			     METRIC[GXZ],
			     METRIC[GYY],
			     METRIC[GYZ],
			     METRIC[GZZ]}, {METRIC[SHIFTX], METRIC[SHIFTY], METRIC[SHIFTZ]}, METRIC[LAPSE] ) ;
	
	metric_c gamma_cL( { METRICm1[GXX],
			     METRICm1[GXY],
			     METRICm1[GXZ],
			     METRICm1[GYY],
			     METRICm1[GYZ],
			     METRICm1[GZZ]}, {METRICm1[SHIFTX], METRICm1[SHIFTY], METRICm1[SHIFTZ]}, METRICm1[LAPSE] ) ;
	
	closure_t crc( { in_prims[VELX].gf[index], in_prims[VELY].gf[index], in_prims[VELZ].gf[index]}, &(gamma_cR) ) ;
	closure_t clc( { in_prims[VELX].gf[indexm1], in_prims[VELY].gf[indexm1], in_prims[VELZ].gf[indexm1]}, &(gamma_cL) ) ; 
	
	crc.update_closure<m1_closure_f>( Er,
					  { Fxr,
					    Fyr,
					    Fzr }, in_prims[ZETA].gf[index], true ) ;
	
	clc.update_closure<m1_closure_f>( El,
					  { Fxl,
					    Fyl,
					    Fzl }, in_prims[ZETA].gf[indexm1], true ) ;
	
	
	cmax = 1.;
	cmin = 1.;
	A    = 1.;
	
	clc.get_fluxes<flux_dirn>(Nl, Nfl, Efl, Fxfl, Fyfl, Fzfl ) ;
	crc.get_fluxes<flux_dirn>(Nr, Nfr, Efr, Fxfr, Fyfr, Fzfr ) ;

	auto const sqrtgamma_cR = gamma_cR.sqrtdet ;
	auto const sqrtgamma_cL = gamma_cL.sqrtdet ;
	
	Enu_flux_LF[index] = hlle_E_flux( Er*sqrtgamma_cR, El*sqrtgamma_cL, Efr, Efl ) ;
	Nnu_flux_LF[index] = hlle_E_flux( Nr*sqrtgamma_cR*crc.Gamma, Nl*sqrtgamma_cL*clc.Gamma, Nfr, Nfl ) ;

	Fnu_x_flux_LF[index] = hlle_F_flux(Fxr*sqrtgamma_cR, Fxl*sqrtgamma_cL, Fxfr, Fxfl ) ;
	Fnu_y_flux_LF[index] = hlle_F_flux(Fyr*sqrtgamma_cR, Fyl*sqrtgamma_cL, Fyfr, Fyfl ) ;
	Fnu_z_flux_LF[index] = hlle_F_flux(Fzr*sqrtgamma_cR, Fzl*sqrtgamma_cL, Fzfr, Fzfl ) ;

	#endif
	
	// ==============================================================================================================
	// Add source terms ( if we're not in the ghostzones ). This routine is called once per direction,
	// therefore here we only add directional source terms and we do this one direction at a time.
	// In particular the F_k source is only computed for k == flux_dirn and the Enu contraction term is
	// added one index at a time 
	// ==============================================================================================================
	if(k<cctk_lsh[2]-cctk_nghostzones[2] && j<cctk_lsh[1]-cctk_nghostzones[1] && i<cctk_lsh[0]-cctk_nghostzones[0] //){
	 && k>=cctk_nghostzones[2] && j>=cctk_nghostzones[1] && i>=cctk_nghostzones[0]) {

	  auto Ec = in_prims[ENU].gf[index] ;
	  
	  std::array<CCTK_REAL,3> F = { in_prims[FNUX].gf[index] * Ec ,
	                                in_prims[FNUY].gf[index] * Ec ,
	                                in_prims[FNUZ].gf[index] * Ec };
	  
	  // cell-centered metric 
	  metric_c gamma_c( {METRIC[GXX], METRIC[GXY], METRIC[GXZ], METRIC[GYY], METRIC[GYZ], METRIC[GZZ]},
			    {METRIC[SHIFTX], METRIC[SHIFTY], METRIC[SHIFTZ]}, METRIC[LAPSE] ) ;
	   
	  CCTK_REAL sqrt_gamma = gamma_c.sqrtdet;
	  
	  // Compute \partial_i g_{\mu \nu} at m+1/2
	  std::array<CCTK_REAL,6> partial_i_gmunu;
	  
	  partial_i_gmunu[GXX] = (FACEVALp1[GXX] - FACEVAL[GXX])*dxi[flux_dirn];
	  partial_i_gmunu[GXY] = (FACEVALp1[GXY] - FACEVAL[GXY])*dxi[flux_dirn];
	  partial_i_gmunu[GXZ] = (FACEVALp1[GXZ] - FACEVAL[GXZ])*dxi[flux_dirn];
	  partial_i_gmunu[GYY] = (FACEVALp1[GYY] - FACEVAL[GYY])*dxi[flux_dirn];
	  partial_i_gmunu[GYZ] = (FACEVALp1[GYZ] - FACEVAL[GYZ])*dxi[flux_dirn];
	  partial_i_gmunu[GZZ] = (FACEVALp1[GZZ] - FACEVAL[GZZ])*dxi[flux_dirn];
	  
	 // Partial shift derivative
	  std::array<CCTK_REAL,3> partial_i_shift;

	  partial_i_shift[0] = (FACEVALp1[SHIFTX] - FACEVAL[SHIFTX])*dxi[flux_dirn];
	  partial_i_shift[1] = (FACEVALp1[SHIFTY] - FACEVAL[SHIFTY])*dxi[flux_dirn];
	  partial_i_shift[2] = (FACEVALp1[SHIFTZ] - FACEVAL[SHIFTZ])*dxi[flux_dirn];

	  // Needed for tau_rhs computation:
	  std::array<CCTK_REAL,3> lapse_deriv = { 0,0,0 };

	  lapse_deriv[flux_dirn] = (FACEVALp1[LAPSE] - FACEVAL[LAPSE])*dxi[flux_dirn] ; // CM 
	  
	  // ----------------------------------------------------------------------------
	  // ----------------------------- FNU source term -----------------------------
	  // ----------------------------------------------------------------------------
	  auto Fdbeta = gamma_c.scalar_product(F,partial_i_shift) ;
	  
	  std::array<double,6> P ;
	  for( int ii=0; ii<6; ii++)
	    P[ii]=Pnu[ii][index] ;
	  auto Pdgmunu = gamma_c.contract_two_tensor_symm(P, partial_i_gmunu ) ;
	  
	  //Source sqrt_gamma * [ -E*d_k(alpha) + alpha/2 * P^ij d_k(gamma_ij) + F_i \partial_k \beta^i ]
	  std::array<CCTK_REAL,3> F_i_curvature_terms = { 0,0,0 };
	  F_i_curvature_terms[flux_dirn] = -Ec*lapse_deriv[flux_dirn]
	    + Fdbeta
	    + 0.5 * gamma_c.lapse * Pdgmunu ;

	  // Notice that Fnue_i_curvature_terms[N]=0 for N!=flux_dirn.
	  // these terms are computed in a directionally split way 
	  Fnu_x_rhs[index] += F_i_curvature_terms[0] * sqrt_gamma;
	  Fnu_y_rhs[index] += F_i_curvature_terms[1] * sqrt_gamma;
	  Fnu_z_rhs[index] += F_i_curvature_terms[2] * sqrt_gamma;

	  // ----------------------------------------------------------------------------
	  // ----------------------------- ENU source term -----------------------------
	  // ----------------------------------------------------------------------------
	  auto F_up = gamma_c.raise_index<0,3>(F) ;
	  // - \alpha sqrt_gamma F^i \partial_i ln( \alpha )
	  Enu_rhs[index]  -= sqrt_gamma * F_up[flux_dirn]*lapse_deriv[flux_dirn] ;
	}

      }


#ifdef FRANKFURT_M1_USE_PPLIM
  // ==============================================================================================================
  // Here we mix the previosly computed LF fluxes into the high order HLLE fluxes.
  // See e.g.  Hu Adams Shu (https://arxiv.org/abs/1203.1540)
  // For M1 the only additional caveat is that we don't mix the diffusive flux in the optically thick regime
  // as this would mess up the asymptotically preserving flux. Therefore our
  //                                         \phi = (1. - \theta)A
  // And 
  //                                          F = \phi F_LO + (1.-\phi) ( F_HO - F_LO ) 
  // ==============================================================================================================
  const CCTK_REAL a2CFL=  6.*(cctkGH->cctk_levfac[flux_dirn-1]/cctkGH->cctk_delta_space[flux_dirn-1])*cctkGH->cctk_delta_time/cctkGH->cctk_timefac;
  #pragma omp parallel for collapse(3)
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2] - cctk_nghostzones[2]+ (flux_dirn==2);k++) 
    for(int j=cctk_nghostzones[1];j<cctk_lsh[1]- cctk_nghostzones[1]+ (flux_dirn==1);j++) 
      for(int i=cctk_nghostzones[0];i<cctk_lsh[0] - cctk_nghostzones[0]+ (flux_dirn==0);i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	int indexm1 = CCTK_GFINDEX3D(cctkGH,i-  kronecker_delta(flux_dirn,0),j-  kronecker_delta(flux_dirn,1),k-  kronecker_delta(flux_dirn,2));
	
	double theta = 1. ;

	auto const kappa_bar = std::sqrt( ( kappa_nu_a[index] + kappa_nu_s[index] ) * ( kappa_nu_a[indexm1] + kappa_nu_s[indexm1]) ) ; 
	auto A  = std::min(1., 1./(kappa_bar+1.e-45) * dxi[flux_dirn]);

	metric_c gamma( { 0.5*(metric[GXX][index]+metric[GXX][indexm1]),
			  0.5*(metric[GXY][index]+metric[GXY][indexm1]),
			  0.5*(metric[GXZ][index]+metric[GXZ][indexm1]),
			  0.5*(metric[GYY][index]+metric[GYY][indexm1]),
			  0.5*(metric[GYZ][index]+metric[GYZ][indexm1]),
			  0.5*(metric[GZZ][index]+metric[GZZ][indexm1])},
	  { 0.5*(metric[SHIFTX][index]+metric[SHIFTX][indexm1]),
	    0.5*(metric[SHIFTY][index]+metric[SHIFTY][indexm1]),
	    0.5*(metric[SHIFTZ][index]+metric[SHIFTZ][indexm1]) }, 0.5*(metric[LAPSE][index]+metric[LAPSE][indexm1]) ) ;

	metric_c gamma_r( { metric[GXX][index],
			    metric[GXY][index],
			    metric[GXZ][index],
			    metric[GYY][index],
			    metric[GYZ][index],
			    metric[GZZ][index],},
	  { metric[SHIFTX][index],
	    metric[SHIFTY][index],
	    metric[SHIFTZ][index],}, metric[LAPSE][index]) ;

	metric_c gamma_l( { metric[GXX][indexm1],
			    metric[GXY][indexm1],
			    metric[GXZ][indexm1],
			    metric[GYY][indexm1],
			    metric[GYZ][indexm1],
			    metric[GZZ][indexm1],},
	  { metric[SHIFTX][indexm1],
	    metric[SHIFTY][indexm1],
	    metric[SHIFTZ][indexm1],}, metric[LAPSE][indexm1]) ;



	CCTK_REAL Er  = in_prims[ENU].gf[index] ;
	
	CCTK_REAL El  = in_prims[ENU].gf[indexm1] ;
		
	CCTK_REAL const sqrtg = gamma_r.sqrtdet ;
	CCTK_REAL const sqrtgm1 = gamma_l.sqrtdet ;

	CCTK_REAL Enu_min = sqrtg*E_atmo;
	CCTK_REAL Enu_minm1 = sqrtgm1 * E_atmo ;

	CCTK_REAL const Enu_LLF_m = sqrtg*Er + a2CFL*Enu_flux_LF[index] ;
	CCTK_REAL const Enu_LLF_p = sqrtgm1*El - a2CFL*Enu_flux_LF[index] ;

	CCTK_REAL const Enu_m = sqrtg*Er + a2CFL*Enu_flux[index] ;
	CCTK_REAL const Enu_p = sqrtgm1*El - a2CFL*Enu_flux[index] ;
	
	
	CCTK_REAL theta_m=1., theta_p=1.;
	
	if( Enu_m < Enu_min )
	  theta_m = MIN(theta, MAX( 0.,(Enu_min-Enu_LLF_m)/(a2CFL*(Enu_flux[index]-Enu_flux_LF[index])))) ;
	if( Enu_p < Enu_minm1 )
	  theta_p = MIN(theta, MAX( 0.,-(Enu_min-Enu_LLF_p)/(a2CFL*(Enu_flux[index]-Enu_flux_LF[index]))));

	theta = MIN(theta_m,theta_p) ;

	CCTK_REAL const phi=  (1.-theta)*A ;

	Nnu_flux[index] = phi*Nnu_flux_LF[index] + (1.-phi)*Nnu_flux[index] ;
	Enu_flux[index] = phi*Enu_flux_LF[index] + (1.-phi)*Enu_flux[index] ;
	Fnu_x_flux[index] = phi*Fnu_x_flux_LF[index] + (1.-phi)*Fnu_x_flux[index] ;
	Fnu_y_flux[index] = phi*Fnu_y_flux_LF[index] + (1.-phi)*Fnu_y_flux[index] ;
	Fnu_z_flux[index] = phi*Fnu_z_flux_LF[index] + (1.-phi)*Fnu_z_flux[index] ;
	  
      }

#endif
  
  
  // ==============================================================================================================
  // Finally add flux derivatives to RHSs. This is done with a simple first order 2 point stencil
  // Notice that we use  F_{i} - F_{i+1} and therefore the above loop was extended one point into the ghostzones
  // in the direction of the flux computation 
  // ==============================================================================================================
#pragma omp parallel for collapse(3)
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++)
    for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++)
      for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
	
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+kronecker_delta(flux_dirn,0),j+kronecker_delta(flux_dirn,1),k+kronecker_delta(flux_dirn,2));
	
	Nnu_rhs[index] += (Nnu_flux[index] - Nnu_flux[indexp1]) * dxi[flux_dirn];

	Enu_rhs[index] += (Enu_flux[index] - Enu_flux[indexp1]) * dxi[flux_dirn];
	
	Fnu_x_rhs[index] += (Fnu_x_flux[index] - Fnu_x_flux[indexp1]) * dxi[flux_dirn];	
	Fnu_y_rhs[index] += (Fnu_y_flux[index] - Fnu_y_flux[indexp1]) * dxi[flux_dirn];
	Fnu_z_rhs[index] += (Fnu_z_flux[index] - Fnu_z_flux[indexp1]) * dxi[flux_dirn];
	
	/*
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	int index_minus[2*(DER_order+1)];
	
	for(int nn=-DER_order;nn<DER_order+2;++nn){
            index_minus[DER_order+nn] = CCTK_GFINDEX3D(cctkGH,i+nn*kronecker_delta[flux_dirn][0],j+nn*kronecker_delta[flux_dirn][1],k+nn*kronecker_delta[flux_dirn][2]);
	}
	#ifdef _M1_N_EVOL__
	Nnue_rhs[index] -= DER_comb<DER_order>(Nnue_flux,index_minus)*dxi[flux_dirn];
	Nnua_rhs[index] -= DER_comb<DER_order>(Nnua_flux,index_minus)*dxi[flux_dirn];
	Nnux_rhs[index] -= DER_comb<DER_order>(Nnux_flux,index_minus)*dxi[flux_dirn];
	#endif
	Enue_rhs[index] -= DER_comb<DER_order>(Enue_flux,index_minus)*dxi[flux_dirn];
	Enua_rhs[index] -= DER_comb<DER_order>(Enua_flux,index_minus)*dxi[flux_dirn];
	Enux_rhs[index] -= DER_comb<DER_order>(Enux_flux,index_minus)*dxi[flux_dirn];

	Fnue_x_rhs[index] -= DER_comb<DER_order>(Fnue_x_flux,index_minus)*dxi[flux_dirn];
	Fnue_y_rhs[index] -= DER_comb<DER_order>(Fnue_y_flux,index_minus)*dxi[flux_dirn];
	Fnue_z_rhs[index] -= DER_comb<DER_order>(Fnue_z_flux,index_minus)*dxi[flux_dirn];
	
	Fnua_x_rhs[index] -= DER_comb<DER_order>(Fnua_x_flux,index_minus)*dxi[flux_dirn];
	Fnua_y_rhs[index] -= DER_comb<DER_order>(Fnua_y_flux,index_minus)*dxi[flux_dirn];
	Fnua_z_rhs[index] -= DER_comb<DER_order>(Fnua_z_flux,index_minus)*dxi[flux_dirn];

	Fnux_x_rhs[index] -= DER_comb<DER_order>(Fnux_x_flux,index_minus)*dxi[flux_dirn];
	Fnux_y_rhs[index] -= DER_comb<DER_order>(Fnux_y_flux,index_minus)*dxi[flux_dirn];
	Fnux_z_rhs[index] -= DER_comb<DER_order>(Fnux_z_flux,index_minus)*dxi[flux_dirn];
	*/
      }
}
