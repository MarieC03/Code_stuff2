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
    // ------------------------------------------------------------------
    //  1. Run the standard GRMHD ID setup (TOV, FUKA, etc.) but with
    //     leptonic_eos_4d_t as the EOS type.  This fills all the
    //     standard variables (DENS_, TAU_, SX_, YESTAR_, ENTROPYSTAR_,
    //     RHO_, YE_, TEMP_, EPS_, PRESS_, ENTROPY_, B, zvec).
    // ------------------------------------------------------------------
    GRACE_INFO("Setting GRMHD initial data for 4D leptonic EOS.") ;
    set_grmhd_initial_data<leptonic_eos_4d_t>() ;

    // ------------------------------------------------------------------
    //  2. Retrieve the device EOS object.
    // ------------------------------------------------------------------
    auto eos = get_leptonic_4d_eos() ;

    // ------------------------------------------------------------------
    //  3. Host-side beta-eq for the atmosphere: determine (ye_atm, ymu_atm)
    //     at (rho_fl, temp_fl) and update the EOS atmosphere values.
    // ------------------------------------------------------------------
    bool const atm_beta_eq =
        grace::get_param<bool>("eos","leptonic_4d","atmosphere_beta_eq", true) ;

    if (atm_beta_eq) {
        double const rho_fl  = grace::get_param<double>("grmhd","atmosphere","rho_fl") ;
        double const temp_fl = grace::get_param<double>("grmhd","atmosphere","temp_fl") ;

        double ye_atm{eos.eos_yemin}, ymu_atm{eos.eos_ymumin} ;
        find_beta_eq_ye_ymu(ye_atm, ymu_atm, rho_fl, temp_fl, eos) ;

        GRACE_INFO("Atmosphere beta-eq: Y_e={:.4f}  Y_mu={:.4f}  "
                   "at rho={:.3e}  T={:.3e}",
                   ye_atm, ymu_atm, rho_fl, temp_fl) ;

        // Update the atmosphere muon fraction stored in the EOS
        // (used by reset_to_atmosphere_lep4d).
        // We do NOT mutate the singleton here; instead the value
        // is captured in the kernel lambda below.
        eos._c2p_ymu_atm = ymu_atm ;
    }

    // ------------------------------------------------------------------
    //  4. Populate YMUSTAR_ and YMU_ from beta equilibrium at each point.
    //     This is done on the GPU with the device-friendly kernel helper.
    // ------------------------------------------------------------------
    auto& vars = grace::variables::get() ;
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
