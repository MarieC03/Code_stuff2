/**
 * @file grmhd_lep4d_id.cpp
 * @brief  Implementation of set_grmhd_initial_data_lep4d().
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifdef GRACE_ENABLE_LEPTONIC_4D

#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/eos/eos_storage_lep4d.hh>
#include <grace/physics/grmhd.hh>
#include <grace/physics/id/grmhd_lep4d_id.hh>
#include <grace/data_structures/variables.hh>
#include <grace/utils/grace_utils.hh>
#include <grace/system/grace_system.hh>

namespace grace {

void set_grmhd_initial_data_lep4d()
{
    GRACE_INFO("Setting GRMHD initial data for 4D leptonic EOS.") ;
    set_grmhd_initial_data<leptonic_eos_4d_t>() ;

    auto eos = get_leptonic_4d_eos() ;
    GRACE_INFO("Leptonic atmosphere fractions: Y_e={:.4f}  Y_mu={:.4f}",
               eos.ye_atmosphere(), eos.get_c2p_ymu_atm()) ;

    auto& vars = grace::variable_list::get() ;
    set_ymustar_from_beta_eq(eos, vars) ;

    GRACE_INFO("Done setting 4D leptonic initial data.") ;
}

} /* namespace grace */

// ------------------------------------------------------------------
//  Explicit template instantiation for leptonic_eos_4d_t in the
//  standard GRMHD pipeline (conservs_to_prims, set_grmhd_initial_data).
// ------------------------------------------------------------------
#define INSTANTIATE_TEMPLATE(EOS) \
template \
void GRACE_HOST_DEVICE \
grace::conservs_to_prims_lep4d<EOS>( grace::lep4d_cons_array_t&  \
                                   , grace::lep4d_prims_array_t&  \
                                   , grace::metric_array_t const& \
                                   , EOS const&                   \
                                   , grace::atmo_params_t const&  \
                                   , grace::excision_params_t const& \
                                   , grace::c2p_params_t const&   \
                                   , double*                       \
                                   , grace::c2p_err_t& )

INSTANTIATE_TEMPLATE(grace::leptonic_eos_4d_t) ;
#undef INSTANTIATE_TEMPLATE

#endif /* GRACE_ENABLE_LEPTONIC_4D */
