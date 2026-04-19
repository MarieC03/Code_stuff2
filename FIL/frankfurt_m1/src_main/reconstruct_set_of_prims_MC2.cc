/*****************************************
 * WENO Reconstruction interface
 * written in Feb. 2017 by Elias R. Most
 *
 * based on the
 *
 * PPM Reconstruction Interface.
 * Zachariah B. Etienne (2013)
 *
 * This implements the WENO-Z using
 * Finite Volume stencils as in
 * Borges et al 2008.
 *
 *****************************************/
#ifndef REC_INDEX
#define REC_INDEX
static constexpr int MINUS2=0;
static constexpr int MINUS1=1;
static constexpr int PLUS0=2;
static constexpr int PLUS1=3;
static constexpr int PLUS2=4;
static constexpr int MAXNUMINDICES=5;
#endif
//    ^^^^^^^^^^^^^ Be _sure_ to define MAXNUMINDICES appropriately!


// You'll find the #define's for LOOP_DEFINE and SET_INDEX_ARRAYS inside:
#include "loop_defines_reconstruction.h"
#include "REC_helpers.hh"


static void reconstruct_set_of_prims_MC2(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,const int num_prims_to_reconstruct,const int *which_prims_to_reconstruct,
					    gf_and_gz_struct *in_prims,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l, CCTK_REAL *temporary) {

  CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES], Ur[MAXNUMVARS][MAXNUMINDICES],Ul[MAXNUMVARS][MAXNUMINDICES];
  int ijkgz_lo_hi[4][2];


  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    int whichvar=which_prims_to_reconstruct[ww];

    if(in_prims[whichvar].gz_lo[flux_dirn]!=0 || in_prims[whichvar].gz_hi[flux_dirn]!=0) {
      CCTK_VError(VERR_DEF_PARAMS,"TOO MANY GZ'S! WHICHVAR=%d: %d %d %d : %d %d %d DIRECTION %d",whichvar,
		  in_prims[whichvar].gz_lo[1],in_prims[whichvar].gz_lo[2],in_prims[whichvar].gz_lo[3],
		  in_prims[whichvar].gz_hi[1],in_prims[whichvar].gz_hi[2],in_prims[whichvar].gz_hi[3],flux_dirn);
    }

    // *** LOOP 1: Interpolate to Ur and Ul, which are face values ***
    //  You will find that Ur depends on U at MINUS1,PLUS0, PLUS1,PLUS2, and
    //                     Ul depends on U at MINUS2,MINUS1,PLUS0,PLUS1.
    //  However, we define the below loop from MINUS2 to PLUS2. Why not split
    //     this up and get additional points? The reason is that later on,
    //     Ur and Ul depend on ftilde, which is defined from MINUS2 to PLUS2,
    //     so we would lose those points anyway.
  
    LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      /*
     * Loop over gridpoints starting from 2 to cctk_lsh[flux_dirn] - 2
     * and from 0 to cctk_lsh in other directions
     */
      
      SET_INDEX_ARRAYS(-2,2,flux_dirn);
      /*
       * Create int index_arr[4][5] containing the 3D indices of 
       * the gridpoint , the gridpoint \pm 1 and the gridpoint \pm 2 
       * in the flux_dirn direction 
       */
     
      
      /* *** LOOP 1a: READ INPUT *** */
      // Read in a primitive at all gridpoints between m = MINUS2 & PLUS2, where m's direction is given by flux_dirn. Store to U.
      for(int ii=MINUS2;ii<=PLUS2;ii++) U[whichvar][ii] = in_prims[whichvar].gf[index_arr[flux_dirn][ii]];

      // Ur[PLUS0] represents U(i+1/2)
	Ur[whichvar][PLUS0]= mc2_REC<false,0>(U,whichvar);
      // Ul[PLUS0] represents U(i-1/2)
	Ul[whichvar][PLUS0]= mc2_REC<true,0>(U,whichvar);



      /* *** LOOP 1c: WRITE OUTPUT *** */
      // Store right face values to {rho_br,Pr,vxr,vyr,vzr,Bxr,Byr,Bzr},
      //    and left face values to {rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl}
      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ur[whichvar][PLUS0];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ul[whichvar][PLUS0];
    }

    // Ur depends on ftilde, which depends on points of U between MINUS2 and PLUS2
    out_prims_r[whichvar].gz_lo[flux_dirn]+=2;
    out_prims_r[whichvar].gz_hi[flux_dirn]+=2;
    // Ul depends on ftilde, which depends on points of U between MINUS2 and PLUS2
    out_prims_l[whichvar].gz_lo[flux_dirn]+=2;
    out_prims_l[whichvar].gz_hi[flux_dirn]+=2;
  }

  // *** LOOP 4: SHIFT Ur AND Ul ***
  /* Currently face values are set so that
   *      a) Ur(i) represents U(i+1/2), and
   *      b) Ul(i) represents U(i-1/2)
   *    Here, we shift so that the indices are consistent:
   *      a) U(i-1/2+epsilon) = oldUl(i)   = newUr(i)
   *      b) U(i-1/2-epsilon) = oldUr(i-1) = newUl(i)
   *    Note that this step is not strictly necessary if you keep
   *      track of indices when computing the flux. */
  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    int whichvar=which_prims_to_reconstruct[ww];
    LOOP_DEFINE(3,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(-1,0,flux_dirn);
      temporary[index_arr[flux_dirn][PLUS0]] = out_prims_r[whichvar].gf[index_arr[flux_dirn][MINUS1]];
    }

    LOOP_DEFINE(3,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(0,0,flux_dirn);
      // Then shift so that Ur represents the gridpoint at i-1/2+epsilon,
      //                and Ul represents the gridpoint at i-1/2-epsilon.
      // Ur(i-1/2) = Ul(i-1/2)     = U(i-1/2+epsilon)
      // Ul(i-1/2) = Ur(i+1/2 - 1) = U(i-1/2-epsilon)
      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = temporary[index_arr[flux_dirn][PLUS0]];
    }
    // Ul was just shifted, so we lost another ghostzone.
    out_prims_l[whichvar].gz_lo[flux_dirn]+=1;
    out_prims_l[whichvar].gz_hi[flux_dirn]+=0;
    // As for Ur, we didn't need to get rid of another ghostzone,
    //    but we did ... seems wasteful!
    out_prims_r[whichvar].gz_lo[flux_dirn]+=1;
    out_prims_r[whichvar].gz_hi[flux_dirn]+=0;

  }
}




