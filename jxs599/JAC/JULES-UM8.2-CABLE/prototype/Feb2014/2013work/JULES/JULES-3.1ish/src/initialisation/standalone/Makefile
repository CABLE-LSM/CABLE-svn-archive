# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(init_drive.o) $(ARC)(init.o) $(ARC)(init_time.o) $(ARC)(init_output.o) $(ARC)(init_urban.o) $(ARC)(init_parms.o) $(ARC)(init_imogen.o) $(ARC)(allocate_jules_arrays.o) $(ARC)(init_model_levels.o) $(ARC)(init_switches.o) $(ARC)(init_prescribed_data.o) $(ARC)(init_vars_tmp.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(init_drive.o): $(ARC)(switches.o) $(ARC)(logging_mod.o) $(ARC)(templating_mod.o) $(ARC)(datetime_mod.o) $(ARC)(string_utils_mod.o) $(ARC)(update_mod.o) $(ARC)(ancil_info.o) $(ARC)(io_constants.o) $(ARC)(model_interface_mod.o) $(ARC)(time_varying_input_mod.o)
$(ARC)(init.o): $(ARC)(logging_mod.o) $(ARC)(init_grid_mod.o) $(ARC)(init_params_mod.o) $(ARC)(dump_mod.o) $(ARC)(spinup_mod.o) $(ARC)(model_time_mod.o) $(ARC)(init_ancillaries_mod.o) $(ARC)(time_varying_input_mod.o) $(ARC)(initial_conditions_mod.o)
$(ARC)(init_time.o): $(ARC)(switches.o) $(ARC)(logging_mod.o) $(ARC)(datetime_mod.o) $(ARC)(string_utils_mod.o) $(ARC)(io_constants.o) $(ARC)(spinup_mod.o) $(ARC)(model_interface_mod.o) $(ARC)(trifctl.o) $(ARC)(model_time_mod.o)
$(ARC)(init_output.o): $(ARC)(output_mod.o) $(ARC)(logging_mod.o) $(ARC)(datetime_mod.o) $(ARC)(string_utils_mod.o) $(ARC)(theta_field_sizes.o) $(ARC)(ancil_info.o) $(ARC)(io_constants.o) $(ARC)(model_interface_mod.o) $(ARC)(model_time_mod.o)
$(ARC)(init_urban.o): $(ARC)(c_z0h_z0m.o) $(ARC)(switches.o) $(ARC)(logging_mod.o) $(ARC)(nstypes.o) $(ARC)(templating_mod.o) $(ARC)(urban_param.o) $(ARC)(nvegparm.o) $(ARC)(switches_urban.o) $(ARC)(string_utils_mod.o) $(ARC)(ancil_info.o) $(ARC)(io_constants.o) $(ARC)(model_interface_mod.o) $(ARC)(input_mod.o)
$(ARC)(init_parms.o): $(ARC)(switches.o) $(ARC)(u_v_grid.o) $(ARC)(logging_mod.o) $(ARC)(fluxes.o) $(ARC)(p_s_parms.o) $(ARC)(ancil_info.o) $(ARC)(coastal.o) $(ARC)(theta_field_sizes.o) $(ARC)(prognostics.o)
$(ARC)(init_imogen.o): $(ARC)(switches.o) $(ARC)(imogen_progs.o) $(ARC)(string_utils_mod.o) $(ARC)(imogen_run.o) $(ARC)(imogen_constants.o) $(ARC)(imogen_time.o) $(ARC)(io_constants.o) $(ARC)(dump_mod.o) $(ARC)(trifctl.o) $(ARC)(model_time_mod.o) $(ARC)(imogen_io_vars.o) $(ARC)(imogen_drive_vars.o) $(ARC)(logging_mod.o) $(ARC)(imogen_map.o) $(ARC)(datetime_mod.o) $(ARC)(imogen_clim.o) $(ARC)(ancil_info.o) $(ARC)(imogen_anlg_vals.o)
$(ARC)(allocate_jules_arrays.o): $(ARC)(c_z0h_z0m.o) $(ARC)(ozone_vars.o) $(ARC)(soil_param.o) $(ARC)(switches.o) $(ARC)(forcing.o) $(ARC)(pftparm.o) $(ARC)(u_v_grid.o) $(ARC)(c_elevate.o) $(ARC)(fluxes.o) $(ARC)(aero.o) $(ARC)(snow_param.o) $(ARC)(trif.o) $(ARC)(jules_mod.o) $(ARC)(dust_parameters_mod_jls.o) $(ARC)(trifctl.o) $(ARC)(coastal.o) $(ARC)(theta_field_sizes.o) $(ARC)(logging_mod.o) $(ARC)(nstypes.o) $(ARC)(screen.o) $(ARC)(nvegparm.o) $(ARC)(surf_param.o) $(ARC)(orog.o) $(ARC)(top_pdm.o) $(ARC)(prognostics.o) $(ARC)(p_s_parms.o) $(ARC)(ancil_info.o)
$(ARC)(init_model_levels.o): $(ARC)(switches.o) $(ARC)(logging_mod.o) $(ARC)(nstypes.o) $(ARC)(ancil_info.o) $(ARC)(io_constants.o) $(ARC)(model_interface_mod.o) $(ARC)(string_utils_mod.o)
$(ARC)(init_switches.o): $(ARC)(switches.o) $(ARC)(logging_mod.o) $(ARC)(datetime_mod.o) $(ARC)(string_utils_mod.o) $(ARC)(io_constants.o)
$(ARC)(init_prescribed_data.o): $(ARC)(switches.o) $(ARC)(logging_mod.o) $(ARC)(templating_mod.o) $(ARC)(datetime_mod.o) $(ARC)(string_utils_mod.o) $(ARC)(update_mod.o) $(ARC)(io_constants.o) $(ARC)(model_interface_mod.o) $(ARC)(time_varying_input_mod.o)
$(ARC)(init_vars_tmp.o): $(ARC)(switches.o) $(ARC)(forcing.o) $(ARC)(c_elevate.o) $(ARC)(c_rough.o) $(ARC)(aero.o) $(ARC)(trifctl.o) $(ARC)(coastal.o) $(ARC)(logging_mod.o) $(ARC)(sea_ice.o) $(ARC)(c_kappai.o) $(ARC)(c_densty.o) $(ARC)(orog.o) $(ARC)(top_pdm.o) $(ARC)(prognostics.o) $(ARC)(p_s_parms.o) $(ARC)(ancil_info.o)


# End of file.
