# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(dpm_rpm_jls.o) $(ARC)(tilepts_jls.o) $(ARC)(lotka_jls.o) $(ARC)(growth_jls.o) $(ARC)(veg-veg1a_jls.o) $(ARC)(phenol_jls.o) $(ARC)(veg-veg2a_jls.o) $(ARC)(compete_jls.o) $(ARC)(vegcarb_jls.o) $(ARC)(decay_jls.o) $(ARC)(triffid_jls.o) $(ARC)(soilcarb_jls.o) $(ARC)(sparm_jls.o) $(ARC)(pft_sparm_jls.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(dpm_rpm_jls.o): $(ARC)(nstypes.o)
$(ARC)(tilepts_jls.o): $(ARC)(switches.o) $(ARC)(nstypes.o)
$(ARC)(lotka_jls.o): $(ARC)(trif.o) $(ARC)(nstypes.o) $(ARC)(pftparm.o) $(ARC)(sigm.o) $(ARC)(seed.o)
$(ARC)(growth_jls.o): $(ARC)(trif.o) $(ARC)(descent.o) $(ARC)(pftparm.o)
$(ARC)(veg-veg1a_jls.o): $(ARC)(switches.o) $(ARC)(nstypes.o)
$(ARC)(phenol_jls.o): $(ARC)(trif.o) $(ARC)(pftparm.o)
$(ARC)(veg-veg2a_jls.o): $(ARC)(switches.o) $(ARC)(nstypes.o) $(ARC)(veg_param.o) $(ARC)(theta_field_sizes.o) $(ARC)(c_mdi.o) $(ARC)(atm_fields_bounds_mod.o) $(ARC)(seed.o) $(ARC)(descent.o)
$(ARC)(compete_jls.o): $(ARC)(seed.o) $(ARC)(descent.o) $(ARC)(nstypes.o)
$(ARC)(vegcarb_jls.o): $(ARC)(trif.o) $(ARC)(pftparm.o)
$(ARC)(decay_jls.o): $(ARC)(descent.o) $(ARC)(csmin.o)
$(ARC)(triffid_jls.o): $(ARC)(trif.o) $(ARC)(nstypes.o) $(ARC)(pftparm.o) $(ARC)(switches.o)
$(ARC)(soilcarb_jls.o): $(ARC)(nstypes.o)
$(ARC)(sparm_jls.o): $(ARC)(dust_param.o) $(ARC)(nstypes.o) $(ARC)(urban_param.o) $(ARC)(nvegparm.o) $(ARC)(switches_urban.o) $(ARC)(pftparm.o) $(ARC)(blend_h.o) $(ARC)(snow_param.o)
$(ARC)(pft_sparm_jls.o): $(ARC)(pftparm.o)


# End of file.
