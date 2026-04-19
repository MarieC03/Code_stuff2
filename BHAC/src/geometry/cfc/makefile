#########################################################################
FOBJECTS += mod_cfc$F mod_cfc_parameters$F mod_cfc_alp$F mod_cfc_beta$F mod_cfc_psi$F mod_gw_br_Iij_Qij$F mod_gw_br$F mod_gw_br_U_br$F mod_gw_br_U_bri$F mod_gw_br_R_br$F mod_gw_br_R2_br$F
OBJECTS += mod_cfc.o mod_cfc_parameters.o mod_cfc_alp.o mod_cfc_beta.o mod_cfc_psi.o mod_gw_br_Iij_Qij.o mod_gw_br.o mod_gw_br_U_br.o mod_gw_br_U_bri.o mod_gw_br_R_br.o mod_gw_br_R2_br.o

#########################################################################
# Dependencies:
#########################################################################
mod_cfc.o : mod_cfc_parameters.o mod_cfc_alp.o mod_cfc_beta.o mod_cfc_psi.o mod_multigrid_coupling.o mod_interpolate.o mod_imhd_intermediate.o mod_small_values.o mod_gw_br.o
mod_cfc_alp.o : mod_imhd_intermediate.o mod_multigrid_coupling.o mod_metric.o mod_eos.o
mod_cfc_beta.o : mod_multigrid_coupling.o mod_metric.o
mod_cfc_psi.o : mod_multigrid_coupling.o mod_metric.o
mod_gw_br_Iij_Qij.o : mod_imhd_intermediate.o mod_eos.o mod_cfc_parameters.o
mod_gw_br.o : mod_imhd_intermediate.o mod_eos.o mod_full_gr_source.o mod_cfc_parameters.o mod_gw_br_U_br.o mod_gw_br_U_bri.o mod_gw_br_Iij_Qij.o mod_gw_br_R_br.o mod_gw_br_R2_br.o mod_full_gr_source.o
mod_gw_br_U_br.o : mod_cfc_parameters.o
mod_gw_br_U_bri.o : mod_cfc_parameters.o mod_imhd_intermediate.o
mod_gw_br_R_br.o : mod_cfc_parameters.o mod_imhd_intermediate.o
mod_gw_br_R2_br.o : mod_cfc_parameters.o


