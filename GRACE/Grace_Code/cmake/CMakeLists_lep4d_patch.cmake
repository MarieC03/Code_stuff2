# CMakeLists_lep4d_patch.cmake
# ============================================================
#  Patch snippet for the GRACE CMakeLists.txt.
#
#  Add these blocks at the appropriate locations in the top-level
#  and src/physics CMakeLists.txt files.
#
#  1. In the top-level CMakeLists.txt:
#     Add an option to enable the 4D leptonic EOS:
#
#       option(GRACE_ENABLE_LEPTONIC_4D
#              "Enable 4D leptonic EOS with muon fraction Y_mu" OFF)
#
#       if(GRACE_ENABLE_LEPTONIC_4D)
#           add_compile_definitions(GRACE_ENABLE_LEPTONIC_4D)
#       endif()
#
#  2. In src/physics/CMakeLists.txt (or src/physics/eos/CMakeLists.txt),
#     inside the target_sources() for the physics library:
#
#       if(GRACE_ENABLE_LEPTONIC_4D)
#           target_sources(grace_physics PRIVATE
#               ${CMAKE_SOURCE_DIR}/src/physics/eos/read_leptonic_4d_table.cpp
#               ${CMAKE_SOURCE_DIR}/src/physics/id/grmhd_lep4d_id.cpp
#               ${CMAKE_SOURCE_DIR}/src/data_structures/variable_registration_lep4d.cpp
#           )
#           target_include_directories(grace_physics PUBLIC
#               ${CMAKE_SOURCE_DIR}/include/grace/physics/eos
#               ${CMAKE_SOURCE_DIR}/include/grace/physics/id
#           )
#       endif()
#
#  3. In test/CMakeLists.txt:
#
#       if(GRACE_ENABLE_LEPTONIC_4D)
#           add_executable(test_leptonic_4d_eos
#               test_leptonic_4d_eos.cpp
#           )
#           target_link_libraries(test_leptonic_4d_eos
#               PRIVATE
#               grace_physics
#               Catch2::Catch2WithMain
#               Kokkos::kokkos
#               ${HDF5_LIBRARIES}
#           )
#           catch_discover_tests(test_leptonic_4d_eos)
#       endif()
#
# ============================================================

# ============================================================
#  Required include directories
# ============================================================
# include/grace/physics/eos/leptonic_eos_4d.hh
# include/grace/physics/eos/read_leptonic_4d_table.hh
# include/grace/physics/eos/kastaun_c2p_lep4d.hh
# include/grace/physics/eos/c2p_lep4d.hh
# include/grace/physics/eos/eos_storage_lep4d.hh
# include/grace/physics/id/grmhd_lep4d_id.hh
# include/grace/data_structures/variable_indices_lep4d_patch.hh
#
# ============================================================
#  Required source files
# ============================================================
# src/physics/eos/read_leptonic_4d_table.cpp
# src/physics/id/grmhd_lep4d_id.cpp
# src/data_structures/variable_registration_lep4d.cpp
#
# ============================================================
#  Modified existing files  (patch instructions)
# ============================================================
#
#  include/grace/data_structures/variable_indices.hh:
#    - In evol_hrsc_var_cc_idx, add after ENTROPYSTAR_:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        YMUSTAR_,
#        #endif
#    - In aux_var_idx, add after YE_:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        YMU_,
#        #endif
#    - At the end, add:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        #include <grace/data_structures/variable_indices_lep4d_patch.hh>
#        #endif
#
#  include/grace/physics/eos/eos_storage.hh:
#    - Add member and accessor per eos_storage_lep4d.hh instructions.
#    - In eos_storage_t constructor add:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        if (eos_type == "leptonic_4d") grace::init_leptonic_4d_eos() ;
#        #endif
#
#  src/data_structures/variable_indices.cpp (register_variables):
#    - After ENTROPYSTAR_ registration, add:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        grace::variables::register_variables_lep4d() ;
#        #endif
#
#  src/physics/grmhd.cpp (set_grmhd_initial_data):
#    - In the id_type dispatch at the end, add:
#        } else if (id_type == "tov" || id_type == "fuka" || ...) {
#            // standard path already handled above
#        }
#        // After all id_type branches, overwrite YMUSTAR_:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        if constexpr (std::is_same_v<eos_t, grace::leptonic_eos_4d_t>) {
#            auto& vars = grace::variables::get() ;
#            grace::set_ymustar_from_beta_eq(
#                grace::get_leptonic_4d_eos(), vars) ;
#        }
#        #endif
#
#  src/physics/grmhd.cpp (flux kernel):
#    - After the YESTAR_ flux line, add:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        fluxes(VEC(i,j,k),YMUSTAR_,idir,q) = ...ymu flux... ;
#        #endif
#    - In the HLL flux section:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        f[YMUSL] = sqrtg * solver(ymul*fdl,ymur*fdr,ymul*densl,ymur*densr,cmin,cmax) ;
#        #endif
#
#  src/physics/eos/c2p.cpp (conservs_to_prims instantiation):
#    - Add:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        INSTANTIATE_TEMPLATE(grace::leptonic_eos_4d_t) ;
#        #endif
#    - And the lep4d-specific instantiation:
#        #ifdef GRACE_ENABLE_LEPTONIC_4D
#        INSTANTIATE_LEP4D_C2P(grace::leptonic_eos_4d_t) ;
#        #endif
#
# ============================================================
#  Example yaml configuration block
# ============================================================
# eos: {
#   eos_type: "leptonic_4d",
#
#   leptonic_4d: {
#     table_filename:      "/path/to/DD2_leptonic4d.h5",
#     cold_table_filename: "/path/to/DD2_cold_lep4d.grace",
#     table_format:        "leptonic_4d_native",
#     do_energy_shift:      true,
#     use_muonic_eos:       true,
#     atmosphere_beta_eq:   true
#   },
#
#   eps_maximum: 1.e+05
# }
