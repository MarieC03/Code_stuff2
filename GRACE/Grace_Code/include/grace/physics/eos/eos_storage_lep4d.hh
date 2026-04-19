/**
 * @file eos_storage_lep4d.hh
 * @brief  Extension of eos_storage_t for the 4D leptonic EOS.
 *
 *         When GRACE_ENABLE_LEPTONIC_4D is defined, the singleton
 *         eos_storage_t gains an additional member _leptonic_4d of type
 *         leptonic_eos_4d_t and accessor methods consistent with the
 *         existing get_tabulated() / get_hybrid_pwpoly() pattern.
 *
 *         Include this header from eos_storage.hh guarded by the macro.
 *
 * @date   2025
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#ifndef GRACE_PHYSICS_EOS_STORAGE_LEP4D_HH
#define GRACE_PHYSICS_EOS_STORAGE_LEP4D_HH

#ifdef GRACE_ENABLE_LEPTONIC_4D

#include <grace/physics/eos/eos_storage.hh>

namespace grace {

/**
 * @brief Convenience wrapper that returns the 4D leptonic EOS object.
 *        Mirrors grace::eos::get().get_tabulated() for the new type.
 */
inline leptonic_eos_4d_t get_leptonic_4d_eos()
{
    return grace::eos::get().get_eos<grace::leptonic_eos_4d_t>();
}

} /* namespace grace */

#endif /* GRACE_ENABLE_LEPTONIC_4D */
#endif /* GRACE_PHYSICS_EOS_STORAGE_LEP4D_HH */
