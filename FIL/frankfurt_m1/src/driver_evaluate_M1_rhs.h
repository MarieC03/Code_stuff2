#ifndef FRANKFURT_M1_RHS_
#define FRANKFURT_M1_RHS_

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
						   CCTK_REAL *Nnu_flux_LO ,
						   CCTK_REAL *Enu_flux_LO, CCTK_REAL *Fnu_x_flux_LO, CCTK_REAL *Fnu_y_flux_LO, CCTK_REAL *Fnu_z_flux_LO) ;


static void compute_E_rhs_extrinsic_curvature_terms_and_Pmunu(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones, CCTK_REAL **metric,
								CCTK_REAL **Pnu,
								gf_and_gz_struct *in_prims,
								CCTK_REAL *kxx,CCTK_REAL *kxy,CCTK_REAL *kxz,CCTK_REAL *kyy,CCTK_REAL *kyz,CCTK_REAL *kzz,
							      CCTK_REAL *Enu_rhs, CCTK_REAL* zeta) ;

static void reconstruct_set_of_prims_MINMOD(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,const int num_prims_to_reconstruct,const int *which_prims_to_reconstruct,
					    gf_and_gz_struct *in_prims,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l, CCTK_REAL *temporary) ;


static void reconstruct_set_of_prims_MC2(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,const int num_prims_to_reconstruct,const int *which_prims_to_reconstruct,
					 gf_and_gz_struct *in_prims,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l, CCTK_REAL *temporary);


static void reconstruct_set_of_prims_WENO5(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,const int num_prims_to_reconstruct,const int *which_prims_to_reconstruct,
					   gf_and_gz_struct *in_prims,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l, CCTK_REAL *temporary);

#endif

