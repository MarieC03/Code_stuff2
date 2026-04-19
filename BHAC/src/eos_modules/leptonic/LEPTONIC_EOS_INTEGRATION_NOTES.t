!================================================================================
!
!  mod_eos_leptonic_makefile_fragment.t
!
!  This file is NOT a compilable Fortran source.  It documents the changes
!  needed in  src/eos_modules/makefile  and in amrvacpar.t / amrvacdef.t to
!  integrate the 4D leptonic EOS modules into BHAC.
!
!================================================================================

!-----------------------------------------------------------------------
! 1.  src/eos_modules/makefile  -- add new objects
!-----------------------------------------------------------------------

! FOBJECTS += mod_eos_leptonic_parameters$F \
!             mod_eos_leptonic_interpolation$F \
!             mod_eos_readtable_leptonic$F \
!             mod_eos_readtable_leptonic_scollapse$F \
!             mod_eos_leptonic$F \
!             mod_imhd_con2prim_leptonic$F
!
! OBJECTS  += mod_eos_leptonic_parameters.o \
!             mod_eos_leptonic_interpolation.o \
!             mod_eos_readtable_leptonic.o \
!             mod_eos_readtable_leptonic_scollapse.o \
!             mod_eos_leptonic.o \
!             mod_imhd_con2prim_leptonic.o

! Dependency chain:
!
! mod_eos_leptonic_parameters.o    : (none beyond intrinsics)
! mod_eos_leptonic_interpolation.o : mod_eos_leptonic_parameters.o
! mod_eos_readtable_leptonic_scollapse.o : mod_eos_leptonic_parameters.o
! mod_eos_readtable_leptonic.o     : mod_eos_leptonic_parameters.o
! mod_eos_leptonic.o               : mod_eos_leptonic_parameters.o \
!                                    mod_eos_leptonic_interpolation.o \
!                                    mod_eos_readtable_leptonic_scollapse.o \
!                                    mod_eos_readtable_leptonic.o \
!                                    mod_eos.o mod_rootfinding.o
! mod_imhd_con2prim_leptonic.o     : mod_eos_leptonic.o \
!                                    mod_imhd_con2prim.o \
!                                    mod_rootfinding.o

!-----------------------------------------------------------------------
! 2.  amrvacpar.t  -- new variable indices (add inside the existing block)
!-----------------------------------------------------------------------

!   {#IFDEF LEPTONIC_EOS
!   integer :: ymu_          = 1000   ! primitive: muon fraction Y_mu
!   integer :: Dymu_         = 1000   ! conserved: D * Y_mu
!   }

!-----------------------------------------------------------------------
! 3.  mod_variables.t  -- register new variables
!-----------------------------------------------------------------------

! Inside setup_variables(), after the Dye_ / ye_ block, add:
!
!   {#IFDEF LEPTONIC_EOS
!   ymu_  = var_set_wvar(prim,  need_rec=.true.,  fill_bc=.true.)
!   Dymu_ = var_set_wvar(cons,  need_rec=.false., fill_bc=.true.)
!   }

!-----------------------------------------------------------------------
! 4.  amrvacphys.t  -- conserven() and primitive()
!-----------------------------------------------------------------------

! In conserven(), after the Dye_ line, add:
!
!   {#IFDEF LEPTONIC_EOS
!   where(.not.patchw(ixO^S))
!     w(ixO^S, Dymu_) = w(ixO^S, d_) * w(ixO^S, ymu_)
!   endwhere
!   }

!-----------------------------------------------------------------------
! 5.  amrvacio/amrio.t  -- parameter reading
!-----------------------------------------------------------------------

! In the EOS select-case block, add:
!
!   case("leptonic")
!     call eos_leptonic_activate()
!
! And add the leptonic namelist block:
!
!   namelist /leptoniceolist/ leptonic_table_name, baryon_table_name, &
!                             baryon_table_type, lep_use_muons, &
!                             lep_add_ele_contribution, lep_fix_ymu_high_yp

!-----------------------------------------------------------------------
! 6.  mod_imhd_con2prim.t  -- call leptonic setup
!-----------------------------------------------------------------------

! In imhd_con2prim_setup(), add at the end:
!
!   {#IFDEF LEPTONIC_EOS
!   call leptonic_con2prim_setup()
!   }

! In imhd_con2prim(), replace the ye_hat recovery block for tabulated EOS:
!
!   if (eos_type == tabulated) then
!     ye_hat = cons_tmp(ix^D, Dye_) / cons_tmp(ix^D, D_)
!     ye_hat = max(eos_yemin, min(eos_yemax, ye_hat))
!   {#IFDEF LEPTONIC_EOS
!     ymu_hat = cons_tmp(ix^D, Dymu_) / cons_tmp(ix^D, D_)
!     ymu_hat = max(lep_eos_ymumin, min(lep_eos_ymumax, ymu_hat))
!   }
!   end if

! And replace  eos_temp_get_all_one_grid / eos_eps_get_all_one_grid calls
! with the leptonic versions, passing the ymu optional argument:
!
!   {#IFDEF LEPTONIC_EOS
!   call leptonic_eps_get_all_one_grid(rho_hat, eps_hat, ye_hat, ymu_hat, &
!       temp_hat, prs=prs_tmp, cs2=cs2_tmp)
!   }
!   {#IFNDEF LEPTONIC_EOS
!   call eos_eps_get_all_one_grid(rho_hat, eps_hat, ye_hat, temp_hat, &
!       prs=prs_tmp, cs2=cs2_tmp)
!   }

!-----------------------------------------------------------------------
! 7.  definition file flags
!-----------------------------------------------------------------------

! Add to the definition file used by the preprocessor (e.g. definitions.f90):
!
!   #define LEPTONIC_EOS
!   #define TABEOS

!-----------------------------------------------------------------------
! 8.  Parameter file example (.par)
!-----------------------------------------------------------------------

!   &eoslist
!     eos_type_input     = "leptonic"
!     baryon_table_name  = "SFHo_baryon.h5"
!     baryon_table_type  = "scollapse"
!     leptonic_table_name = "SFHo_leptonic.h5"
!     lep_use_muons      = .true.
!   /
!
!   &methodlist
!     typeinversion = "KastaunC2P"
!   /
