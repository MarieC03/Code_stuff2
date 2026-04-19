/**
 * @file compute_E_rhs_extrinsic_curvature_terms_and_Pmunu.cc
 *
 * This file is part of frankfurt_m1
 *
 *  Compute curvature terms on the RHS
 *  of the energy density equation
 *
 *   \alpha \sqrt{\gamma} (T^{00} \beta^i \beta ^j + 2T^{0i} \beta^j + T^{ij})*K_{ij}
 *
 *  Based on previous version by L. Weih and E. R. Most.
 *
 * @author  Carlo Musolino
 */ 
static void compute_E_rhs_extrinsic_curvature_terms_and_Pmunu
(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones, CCTK_REAL **metric,
 CCTK_REAL **Pnu,
 gf_and_gz_struct *in_prims,
 CCTK_REAL *kxx,CCTK_REAL *kxy,CCTK_REAL *kxz,CCTK_REAL *kyy,CCTK_REAL *kyz,CCTK_REAL *kzz,
 CCTK_REAL *Enu_rhs, CCTK_REAL* zeta) {
  using namespace frankfurt_m1;
  
  // These loop extents must be consistent with add_fluxes_and_source_terms_to_M1_rhs(), 
  // since we use PUPmunu there as well. On the lower end, the loop is extended by one 
  // index in each direction, because the diffusion-term in the flux-corrections needs 
  // one more cell-centered Jthick for computing d(Jthick)/dx

  //for(int k=cctk_nghostzones[2]-DER_order-1;k<cctk_lsh[2]-(cctk_nghostzones[2]-1);k++) for(int j=cctk_nghostzones[1]-DER_order-1;j<cctk_lsh[1]-(cctk_nghostzones[1]-1);j++) for(int i=cctk_nghostzones[0]-DER_order-1;i<cctk_lsh[0]-(cctk_nghostzones[0]-1);i++) {

  int const imin{cctk_nghostzones[0]}, imax{cctk_lsh[0]-(cctk_nghostzones[0])};
  int const jmin{cctk_nghostzones[1]}, jmax{cctk_lsh[1]-(cctk_nghostzones[1])};
  int const kmin{cctk_nghostzones[2]}, kmax{cctk_lsh[2]-(cctk_nghostzones[2])};
  #pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {

	int index = CCTK_GFINDEX3D(cctkGH,i,j,k) ;
	
	std::array<CCTK_REAL,6> Kij { kxx[index],
	    kxy[index],
	    kxz[index],
	    kyy[index],
	    kyz[index],
	    kzz[index] } ;

	metric_c gamma( { metric[GXX][index],
			  metric[GXY][index],
			  metric[GXZ][index],
			  metric[GYY][index],
			  metric[GYZ][index],
			  metric[GZZ][index] } , { metric[SHIFTX][index], metric[SHIFTY][index], metric[SHIFTZ][index] }, metric[LAPSE][index] ) ; 
	

	auto const alpha_sqrtgamma = gamma.lapse * gamma.sqrtdet ;

	closure_t cl( { in_prims[VELX].gf[index],
	      in_prims[VELY].gf[index],
	      in_prims[VELZ].gf[index] }, &gamma ) ;

	cl.update_closure<m1_closure_f>(in_prims[ENU].gf[index], { in_prims[FNUX].gf[index],
								   in_prims[FNUY].gf[index],
								   in_prims[FNUZ].gf[index] }, zeta[index] , true) ;

	auto const P  = cl.get_pressure() ;
	auto PK = gamma.contract_two_tensor_symm(P, Kij) ;
	
	//  \alpha \sqrt{\gamma} (T^{00} \beta^i \beta ^j + 2T^{0i} \beta^j + T^{ij})*K_{ij}
	Enu_rhs[index] += alpha_sqrtgamma*PK ;

	for(int ii=0; ii<6; ii++)
	  Pnu[ii][index] = P[ii];
	
	
      }
}
