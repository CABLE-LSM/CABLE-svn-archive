# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(delta_temp.o) $(ARC)(gcm_anlg.o) $(ARC)(imogen_check.o) $(ARC)(ocean_co2.o) $(ARC)(invert.o) $(ARC)(imogen_confirmed_run.o) $(ARC)(solang.o) $(ARC)(drdat.o) $(ARC)(diffcarb_land.o) $(ARC)(response.o) $(ARC)(redis.o) $(ARC)(radf_non_co2.o) $(ARC)(sunny.o) $(ARC)(clim_calc.o) $(ARC)(rndm.o) $(ARC)(imogen_update_clim.o) $(ARC)(radf_co2.o) $(ARC)(day_calc.o) $(ARC)(solpos.o) $(ARC)(imogen_update_carb.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(delta_temp.o): 
$(ARC)(gcm_anlg.o): $(ARC)(imogen_constants.o) $(ARC)(imogen_map.o) $(ARC)(io_constants.o)
$(ARC)(imogen_check.o): $(ARC)(logging_mod.o)
$(ARC)(ocean_co2.o): 
$(ARC)(invert.o): $(ARC)(logging_mod.o)
$(ARC)(imogen_confirmed_run.o): $(ARC)(logging_mod.o)
$(ARC)(solang.o): $(ARC)(c_pi.o)
$(ARC)(drdat.o): $(ARC)(io_constants.o)
$(ARC)(diffcarb_land.o): $(ARC)(c_pi.o) $(ARC)(earth_constants_mod_jls.o)
$(ARC)(response.o): 
$(ARC)(redis.o): 
$(ARC)(radf_non_co2.o): $(ARC)(logging_mod.o) $(ARC)(io_constants.o)
$(ARC)(sunny.o): $(ARC)(c_pi.o) $(ARC)(datetime_mod.o)
$(ARC)(clim_calc.o): 
$(ARC)(rndm.o): 
$(ARC)(imogen_update_clim.o): $(ARC)(imogen_drive_vars.o) $(ARC)(logging_mod.o) $(ARC)(imogen_progs.o) $(ARC)(imogen_time.o) $(ARC)(datetime_mod.o) $(ARC)(imogen_run.o) $(ARC)(trifctl.o) $(ARC)(imogen_clim.o) $(ARC)(aero.o) $(ARC)(ancil_info.o) $(ARC)(io_constants.o) $(ARC)(imogen_constants.o) $(ARC)(imogen_anlg_vals.o) $(ARC)(prognostics.o)
$(ARC)(radf_co2.o): 
$(ARC)(day_calc.o): $(ARC)(c_pi.o)
$(ARC)(solpos.o): $(ARC)(c_pi.o)
$(ARC)(imogen_update_carb.o): $(ARC)(imogen_io_vars.o) $(ARC)(logging_mod.o) $(ARC)(imogen_progs.o) $(ARC)(imogen_run.o) $(ARC)(trifctl.o) $(ARC)(imogen_clim.o) $(ARC)(imogen_time.o) $(ARC)(ancil_info.o) $(ARC)(imogen_constants.o) $(ARC)(imogen_anlg_vals.o) $(ARC)(prognostics.o)


# End of file.
