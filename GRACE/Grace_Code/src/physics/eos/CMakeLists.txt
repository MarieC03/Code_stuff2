register_grace_object(grace_c2p c2p.cpp)
#target_compile_options(grace_c2p PRIVATE
#    -fno-fast-math
#    -fno-unsafe-math-optimizations
#    -ffp-model=strict
#)
register_grace_object(grace_eos eos_storage.cpp read_eos_table.cpp)
register_grace_object_conditional(grace_lep4d_eos GRACE_ENABLE_LEPTONIC_4D read_leptonic_4d_table.cpp)
