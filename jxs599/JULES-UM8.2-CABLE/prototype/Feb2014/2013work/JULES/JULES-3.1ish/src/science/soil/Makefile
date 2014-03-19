# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(gauss_jls.o) $(ARC)(darcy_ch.o) $(ARC)(surf_hyd_jls.o) $(ARC)(ice_htc_jls.o) $(ARC)(hydrol_jls.o) $(ARC)(calc_baseflow_jls.o) $(ARC)(hyd_con_vg_jls.o) $(ARC)(soilmc_jls.o) $(ARC)(soil_htc_jls.o) $(ARC)(calc_baseflow_jules.o) $(ARC)(heat_con_jls.o) $(ARC)(ch4_wetl_jls.o) $(ARC)(sieve_jls.o) $(ARC)(calc_fsat.o) $(ARC)(soilt_jls.o) $(ARC)(hyd_con_ch.o) $(ARC)(soil_hyd_jls.o) $(ARC)(darcy_vg_jls.o) $(ARC)(calc_zw_jls.o) $(ARC)(hyd_con_ic.o) $(ARC)(darcy_ic.o) $(ARC)(frunoff_jls.o) $(ARC)(pdm_jls.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(gauss_jls.o): 
$(ARC)(darcy_ch.o): 
$(ARC)(surf_hyd_jls.o): $(ARC)(switches.o) $(ARC)(soil_param.o)
$(ARC)(ice_htc_jls.o): $(ARC)(snow_param.o)
$(ARC)(hydrol_jls.o): $(ARC)(c_topog.o) $(ARC)(soil_param.o)
$(ARC)(calc_baseflow_jls.o): $(ARC)(c_topog.o) $(ARC)(soil_param.o)
$(ARC)(hyd_con_vg_jls.o): 
$(ARC)(soilmc_jls.o): $(ARC)(c_densty.o) $(ARC)(soil_param.o)
$(ARC)(soil_htc_jls.o): $(ARC)(c_0_dg_c.o) $(ARC)(c_lheat.o) $(ARC)(snow_param.o) $(ARC)(soil_param.o) $(ARC)(c_densty.o) $(ARC)(c_perma.o)
$(ARC)(calc_baseflow_jules.o): $(ARC)(switches.o) $(ARC)(c_topog.o) $(ARC)(soil_param.o)
$(ARC)(heat_con_jls.o): $(ARC)(switches.o) $(ARC)(soil_param.o) $(ARC)(snow_param.o)
$(ARC)(ch4_wetl_jls.o): $(ARC)(c_0_dg_c.o) $(ARC)(surf_param.o) $(ARC)(c_ch4.o)
$(ARC)(sieve_jls.o): 
$(ARC)(calc_fsat.o): $(ARC)(c_topog.o)
$(ARC)(soilt_jls.o): $(ARC)(soil_param.o)
$(ARC)(hyd_con_ch.o): 
$(ARC)(soil_hyd_jls.o): $(ARC)(c_topog.o) $(ARC)(c_densty.o) $(ARC)(soil_param.o)
$(ARC)(darcy_vg_jls.o): 
$(ARC)(calc_zw_jls.o): $(ARC)(c_topog.o) $(ARC)(c_densty.o)
$(ARC)(hyd_con_ic.o): $(ARC)(switches.o)
$(ARC)(darcy_ic.o): $(ARC)(switches.o)
$(ARC)(frunoff_jls.o): 
$(ARC)(pdm_jls.o): $(ARC)(switches.o) $(ARC)(c_densty.o) $(ARC)(soil_param.o)


# End of file.
