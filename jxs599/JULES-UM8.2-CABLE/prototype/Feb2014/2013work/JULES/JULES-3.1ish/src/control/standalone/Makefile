# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(next_time.o) $(ARC)(zenith_mod.o) $(ARC)(jules.o) $(ARC)(control.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(next_time.o): $(ARC)(logging_mod.o) $(ARC)(datetime_mod.o) $(ARC)(string_utils_mod.o) $(ARC)(dump_mod.o) $(ARC)(spinup_mod.o) $(ARC)(time_varying_input_mod.o) $(ARC)(trifctl.o) $(ARC)(model_time_mod.o)
$(ARC)(zenith_mod.o): $(ARC)(theta_field_sizes.o) $(ARC)(datetime_mod.o) $(ARC)(c_pi.o) $(ARC)(latlon_mod.o) $(ARC)(model_time_mod.o)
$(ARC)(jules.o): $(ARC)(update_mod.o) $(ARC)(output_mod.o) $(ARC)(time_varying_input_mod.o) $(ARC)(switches.o) $(ARC)(model_time_mod.o)
$(ARC)(control.o): $(ARC)(ozone_vars.o) $(ARC)(p_s_parms.o) $(ARC)(switches.o) $(ARC)(forcing.o) $(ARC)(c_elevate.o) $(ARC)(csigma.o) $(ARC)(fluxes.o) $(ARC)(c_epslon.o) $(ARC)(atm_fields_bounds_mod.o) $(ARC)(snow_param.o) $(ARC)(trifctl.o) $(ARC)(coastal.o) $(ARC)(theta_field_sizes.o) $(ARC)(u_v_grid.o) $(ARC)(sea_ice.o) $(ARC)(zenith_mod.o) $(ARC)(screen.o) $(ARC)(diag_swchs.o) $(ARC)(datetime_mod.o) $(ARC)(surf_param.o) $(ARC)(orog.o) $(ARC)(top_pdm.o) $(ARC)(update_mod.o) $(ARC)(prognostics.o) $(ARC)(bl_option_mod.o) $(ARC)(aero.o) $(ARC)(c_gamma.o) $(ARC)(model_time_mod.o) $(ARC)(ancil_info.o)


# End of file.
