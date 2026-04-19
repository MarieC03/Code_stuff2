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

#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/eos/read_leptonic_4d_table.hh>

namespace grace {

// ---------------------------------------------------------------------------
//  Mixin / patch for eos_storage_t.
//  Because eos_storage_t is a singleton with a private constructor, we
//  cannot derive from it.  Instead we provide free functions that the
//  singleton's constructor calls, and an accessor that wraps the
//  storage singleton.
// ---------------------------------------------------------------------------

/**
 * @brief Global singleton-style accessor for the 4D leptonic EOS.
 *
 * The object is default-constructed (empty) until init_leptonic_4d_eos()
 * is called from eos_storage_t's constructor when eos_type == "leptonic_4d".
 */
inline leptonic_eos_4d_t& leptonic_4d_eos_instance()
{
    static leptonic_eos_4d_t _inst ;
    return _inst ;
}

/**
 * @brief Read the 4D table and populate the global EOS instance.
 *        Called once during initialisation.
 */
inline void init_leptonic_4d_eos()
{
    leptonic_4d_eos_instance() = read_leptonic_4d_table() ;
}

/**
 * @brief Convenience wrapper that returns the 4D leptonic EOS object.
 *        Mirrors  grace::eos::get().get_tabulated()  for the new type.
 */
inline leptonic_eos_4d_t get_leptonic_4d_eos()
{
    return leptonic_4d_eos_instance() ;
}

} /* namespace grace */

// ---------------------------------------------------------------------------
//  Additions to eos_storage_t (added via partial-class extension macro).
//  Paste the block below inside the eos_storage_t class body in
//  eos_storage.hh, guarded by GRACE_ENABLE_LEPTONIC_4D.
//
//    #ifdef GRACE_ENABLE_LEPTONIC_4D
//        //! The 4D leptonic EOS object.
//        leptonic_eos_4d_t _leptonic_4d ;
//
//        decltype(auto) GRACE_ALWAYS_INLINE get_leptonic_4d() {
//            return _leptonic_4d ;
//        }
//    #endif
//
//  And in eos_storage_t::get_eos<eos_t>:
//
//    } else if constexpr (std::is_same_v<eos_t, leptonic_eos_4d_t>) {
//        return _leptonic_4d ;
//    }
//
// ---------------------------------------------------------------------------

#endif /* GRACE_ENABLE_LEPTONIC_4D */
#endif /* GRACE_PHYSICS_EOS_STORAGE_LEP4D_HH */
