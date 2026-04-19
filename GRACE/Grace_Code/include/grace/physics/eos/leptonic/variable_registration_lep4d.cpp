/**
 * @file variable_registration_lep4d.cpp
 * @brief  Registers YMUSTAR_ (evolved) and YMU_ (auxiliary) grid functions
 *         when GRACE_ENABLE_LEPTONIC_4D is defined.
 *
 *         Call register_variables_lep4d() from within register_variables()
 *         in src/data_structures/variable_indices.cpp after the base hydro
 *         variables are registered and before the M1 / metric variables.
 *
 *         Example placement in variable_indices.cpp:
 *
 *           register_evolved_scalar(YESTAR_,      "ye_star",  hydro_bc, "second_order") ;
 *           register_evolved_scalar(ENTROPYSTAR_, "s_star",   hydro_bc, "second_order") ;
 *       #ifdef GRACE_ENABLE_LEPTONIC_4D
 *           register_variables_lep4d() ;
 *       #endif
 *
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifdef GRACE_ENABLE_LEPTONIC_4D

#include <grace/data_structures/variable_indices.hh>
#include <grace/utils/grace_utils.hh>

namespace grace { namespace variables {

void register_variables_lep4d()
{
    auto hydro_bc =
        detail::get_bc_type(grace::get_param<std::string>("grmhd","bc_kind")) ;

    // ------------------------------------------------------------------
    //  Evolved: Y_mu * sqrt(g) * W * rho   (muon fraction conserved)
    // ------------------------------------------------------------------
    register_evolved_scalar(YMUSTAR_, "ymu_star", hydro_bc, "second_order") ;

    // ------------------------------------------------------------------
    //  Auxiliary: Y_mu  (primitive muon fraction)
    // ------------------------------------------------------------------
    register_aux_scalar(YMU_, "ymu") ;
}

}} /* namespace grace::variables */

#endif /* GRACE_ENABLE_LEPTONIC_4D */
