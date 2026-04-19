/**
 * @file variable_indices_lep4d_patch.hh
 * @brief  Preprocessor-guarded addendum to variable_indices.hh.
 *         Include this file immediately after the closing brace of
 *         aux_var_idx in variable_indices.hh when GRACE_ENABLE_LEPTONIC_4D
 *         is defined.
 *
 * Usage in variable_indices.hh:
 *
 *   enum aux_var_idx : int {
 *       RHO_ = 0,
 *       ...
 *       N_AUX_VARS
 *   } ;
 *
 *   #ifdef GRACE_ENABLE_LEPTONIC_4D
 *   #include <grace/data_structures/variable_indices_lep4d_patch.hh>
 *   #endif
 *
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifndef GRACE_DATA_STRUCTURES_VAR_IDX_LEP4D_PATCH_HH
#define GRACE_DATA_STRUCTURES_VAR_IDX_LEP4D_PATCH_HH

// ---------------------------------------------------------------------------
//  Additional evolved conserved variable: Y_mu * sqrt(g) * D
// ---------------------------------------------------------------------------
/**
 * @brief Enum extension: HRSC cell-centred evolved variables for leptonic-4D.
 *        YMUSTAR_ must be appended right after ENTROPYSTAR_ and before any
 *        radiation or metric variables.
 *
 * To activate, redefine evol_hrsc_var_cc_idx inserting YMUSTAR_ after
 * ENTROPYSTAR_.  A clean way is to add it in variable_indices.hh:
 *
 *   enum evol_hrsc_var_cc_idx : int {
 *       DENS_ = 0,
 *       SX_, SY_, SZ_,
 *       TAU_,
 *       YESTAR_,
 *       ENTROPYSTAR_,
 *   #ifdef GRACE_ENABLE_LEPTONIC_4D
 *       YMUSTAR_,        // ← new
 *   #endif
 *       ...
 *       N_HRSC_CC
 *   } ;
 */

// ---------------------------------------------------------------------------
//  Additional auxiliary primitive variable: Y_mu
// ---------------------------------------------------------------------------
/**
 * @brief Enum extension: auxiliary variable index for the muon fraction.
 *        YMU_ must be appended right after YE_ in aux_var_idx.
 *
 * In variable_indices.hh, extend aux_var_idx:
 *
 *   enum aux_var_idx : int {
 *       RHO_ = 0,
 *       ...
 *       YE_,
 *   #ifdef GRACE_ENABLE_LEPTONIC_4D
 *       YMU_,            // ← new
 *   #endif
 *       TEMP_,
 *       ...
 *       N_AUX_VARS
 *   } ;
 */

// ---------------------------------------------------------------------------
//  FILL_CONS_ARRAY / FILL_PRIMS_ARRAY macro overrides for leptonic-4D.
//  These macros extend the standard ones in grmhd_helpers.hh to also
//  populate YMUSL / YMUL when the leptonic-4D mode is enabled.
// ---------------------------------------------------------------------------

#define FILL_CONS_ARRAY_LEP4D(consarr, vview, q, ...)     \
consarr[LEP4D_DENSL]  = vview(__VA_ARGS__,DENS_,q);        \
consarr[LEP4D_TAUL]   = vview(__VA_ARGS__,TAU_,q);         \
consarr[LEP4D_STXL]   = vview(__VA_ARGS__,SX_,q);          \
consarr[LEP4D_STYL]   = vview(__VA_ARGS__,SY_,q);          \
consarr[LEP4D_STZL]   = vview(__VA_ARGS__,SZ_,q);          \
consarr[LEP4D_YESL]   = vview(__VA_ARGS__,YESTAR_,q);      \
consarr[LEP4D_ENTSL]  = vview(__VA_ARGS__,ENTROPYSTAR_,q); \
consarr[LEP4D_YMUSL]  = vview(__VA_ARGS__,YMUSTAR_,q)

#define FILL_PRIMS_ARRAY_LEP4D(primsarr, vview, q, ...)   \
primsarr[LEP4D_RHOL]   = vview(__VA_ARGS__,RHO_,q);        \
primsarr[LEP4D_PRESSL] = vview(__VA_ARGS__,PRESS_,q);      \
primsarr[LEP4D_ZXL]    = vview(__VA_ARGS__,ZVECX_,q);      \
primsarr[LEP4D_ZYL]    = vview(__VA_ARGS__,ZVECY_,q);      \
primsarr[LEP4D_ZZL]    = vview(__VA_ARGS__,ZVECZ_,q);      \
primsarr[LEP4D_YEL]    = vview(__VA_ARGS__,YE_,q);         \
primsarr[LEP4D_YMUL]   = vview(__VA_ARGS__,YMU_,q);        \
primsarr[LEP4D_TEMPL]  = vview(__VA_ARGS__,TEMP_,q);       \
primsarr[LEP4D_EPSL]   = vview(__VA_ARGS__,EPS_,q);        \
primsarr[LEP4D_ENTL]   = vview(__VA_ARGS__,ENTROPY_,q)

// ---------------------------------------------------------------------------
//  Scatter-back macro: write primitives back to auxiliary view.
// ---------------------------------------------------------------------------
#define SCATTER_PRIMS_ARRAY_LEP4D(primsarr, vview, q, ...)  \
vview(__VA_ARGS__,RHO_,q)     = primsarr[LEP4D_RHOL];       \
vview(__VA_ARGS__,PRESS_,q)   = primsarr[LEP4D_PRESSL];     \
vview(__VA_ARGS__,ZVECX_,q)   = primsarr[LEP4D_ZXL];        \
vview(__VA_ARGS__,ZVECY_,q)   = primsarr[LEP4D_ZYL];        \
vview(__VA_ARGS__,ZVECZ_,q)   = primsarr[LEP4D_ZZL];        \
vview(__VA_ARGS__,YE_,q)      = primsarr[LEP4D_YEL];        \
vview(__VA_ARGS__,YMU_,q)     = primsarr[LEP4D_YMUL];       \
vview(__VA_ARGS__,TEMP_,q)    = primsarr[LEP4D_TEMPL];      \
vview(__VA_ARGS__,EPS_,q)     = primsarr[LEP4D_EPSL];       \
vview(__VA_ARGS__,ENTROPY_,q) = primsarr[LEP4D_ENTL]

#endif /* GRACE_DATA_STRUCTURES_VAR_IDX_LEP4D_PATCH_HH */
