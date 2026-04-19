#ifndef FRANKFURT_M1_HEADER_HH
#define FRANKFURT_M1_HEADER_HH

#include "utils/metric.hh"
#include "constexpr_utils.hh"
#include "m1_implicit_conventions.hh"
#include "m1_closure_formulations.hh"

#ifndef CCTK_REAL
#define CCTK_REAL double
#endif 

#ifndef FRANKFURT_M1_USE_PPLIM
#define FRANKFURT_M1_USE_PPLIM
#endif 

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SQR(x) ((x) * (x))
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780

#define VERR_DEF_PARAMS __LINE__, __FILE__, CCTK_THORNSTRING

//! layout memory size
#define UTILS_GFSIZE(CGH)                                                      \
    (CGH->cctk_ash[0]*CGH->cctk_ash[1]*CGH->cctk_ash[2])

// indices
//Harry: added numu, numu_bar for the indices, any problems? FIXME
enum m1_species_t { NUE=0, NUE_BAR, NUMU, NUMU_BAR, NUX, NUMSPECIES} ;
enum m1_eas_index_t { ETA=0, K_A, K_S, ETA_N, K_N, NUM_EAS} ;
enum recon_index_t { NNU=0, ENU, FNUX, FNUY, FNUZ, VELX, VELY, VELZ, ZETA, NUM_RECON_PRIM };
enum metric_index_t { GXX=0, GXY, GXZ, GYY, GYZ, GZZ, LAPSE, SHIFTX, SHIFTY, SHIFTZ, NUMVARS_METRIC} ;
//Harry: added ymu for the indices, any problems? FIXME
// need to add ymu after fixing Marghertia
enum m1_fluid_index_t { RHO=0, YE, YMU, TEMP, REQ_FLUID_VARS} ;


static constexpr const int MAXNUMVARS = NUM_RECON_PRIM ;

// reconstruction helper class
struct gf_and_gz_struct {
  CCTK_REAL *gf;
  int gz_lo[4],gz_hi[4];
};

// FIXME: For cosmetic purposes, we might want to make everything either zero-offset or one-offset, instead of a mixture.
const int kronecker_delta_arr[4][3] = { { 0,0,0 },
					{ 1,0,0 },
					{ 0,1,0 },
					{ 0,0,1 } };


using m1_closure_f = frankfurt_m1::minerbo_closure_f ; 

#endif
