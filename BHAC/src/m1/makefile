# Add the objects for M1:
FOBJECTS += mod_m1$F mod_m1_closure$F mod_m1_closure_functions$F mod_m1_metric_interface$F mod_m1_internal$F mod_m1_tests$F mod_m1_collisional$F mod_m1_eas$F mod_m1_eas_param$F mod_m1_eas_test$F mod_m1_beta_eq$F mod_m1_constants$F mod_m1_fermi$F mod_Weakhub_reader$F mod_weakhub_parameters$F mod_m1_eas_microphysical_gray$F mod_m1_tau_optical$F mod_m1_outflow$F mod_m1_backreaction$F mod_m1_diagnostics$F    
OBJECTS  += mod_m1$O mod_m1_closure$O mod_m1_closure_functions$O mod_m1_metric_interface$O mod_m1_internal$O mod_m1_tests$O mod_m1_collisional$O mod_m1_eas$O mod_m1_eas_param$O mod_m1_eas_test$O mod_m1_beta_eq$O mod_m1_constants$O mod_m1_fermi$O mod_Weakhub_reader$O mod_weakhub_parameters$O mod_m1_eas_microphysical_gray$O mod_m1_tau_optical$O mod_m1_outflow$O mod_m1_backreaction$O mod_m1_diagnostics$O
INCLUDES += mod_m1$F mod_m1_closure$F mod_m1_closure_functions$F mod_m1_metric_interface$F mod_m1_internal$F mod_m1_tests$F mod_m1_collisional$F mod_m1_eas$F mod_m1_eas_param$F mod_m1_eas_test$F mod_m1_beta_eq$F mod_m1_constants$F mod_m1_fermi$F mod_Weakhub_reader$F mod_weakhub_parameters$F mod_m1_eas_microphysical_gray$F mod_m1_tau_optical$F mod_m1_outflow$F mod_m1_backreaction$F mod_m1_diagnostics$F

#########################################################################
# Dependencies:
#########################################################################
mod_m1.o: mod_m1_metric_interface.o mod_m1_closure.o mod_m1_tests.o mod_m1_collisional.o mod_m1_eas.o mod_m1_eas_param.o mod_m1_eas_test.o mod_eos.o mod_m1_backreaction.o
mod_m1_collisional.o : mod_rootfinding.o mod_m1_closure.o mod_m1_internal.o mod_m1_eas.o #mod_metric_interface.o
mod_m1_closure.o: mod_rootfinding.o mod_m1_metric_interface.o mod_m1_closure_functions.o mod_m1_tests.o 
mod_m1_eas.o: amr_ghostcells.o mod_Weakhub_reader.o mod_m1_fermi.o mod_m1_tests.o mod_m1_eas_microphysical_gray.o mod_m1_eas_param.o mod_eos_tabulated.o mod_m1_beta_eq.o mod_m1_constants.o mod_eos.o
mod_m1_eas_microphysical_gray.o: mod_m1_eas_param.o mod_interpolate.o
mod_Weakhub_reader.o: mod_weakhub_parameters.o 
mod_m1_fermi.o: mod_m1_constants.o mod_m1_tests.o mod_m1_internal.o #mod_m1_closure.o 
mod_m1_beta_eq.o : mod_m1_fermi.o mod_rootfinding.o mod_eos_tabulated.o mod_m1_constants.o mod_lu.o mod_lnsrch.o mod_eos.o
mod_m1_tau_optical.o : mod_m1_eas.o mod_m1_eas_microphysical_gray.o mod_m1_eas_param.o mod_m1_metric_interface.o mod_m1_internal.o mod_m1_constants.o
mod_m1_outflow.o : mod_m1_metric_interface.o mod_m1_tests.o mod_m1_eas_param.o mod_m1_constants.o
mod_m1_metric_interface.o: mod_m1_tests.o
mod_m1_backreaction.o: mod_m1_tests.o mod_m1_eas.o mod_m1_eas_test.o mod_eos.o mod_m1_eas_param.o mod_m1_internal.o
mod_m1_diagnostics.o : mod_m1_closure.o mod_m1_fermi.o mod_rootfinding.o mod_eos_tabulated.o mod_eos.o mod_m1_constants.o mod_m1_eas_param.o mod_m1_internal.o

