# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(urbanemis.o) $(ARC)(tile_albedo_jls.o) $(ARC)(canyonalb.o) $(ARC)(albsnow_jls.o) $(ARC)(matinv.o) $(ARC)(ftsa_jls.o) $(ARC)(albpft_jls.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(urbanemis.o): 
$(ARC)(tile_albedo_jls.o): $(ARC)(switches.o) $(ARC)(nstypes.o) $(ARC)(urban_param.o) $(ARC)(nvegparm.o) $(ARC)(switches_urban.o) $(ARC)(pftparm.o) $(ARC)(c_0_dg_c.o) $(ARC)(rad_param.o) $(ARC)(p_s_parms.o) $(ARC)(lake_mod.o) $(ARC)(snow_param.o) $(ARC)(c_mdi.o) $(ARC)(jules_mod.o) $(ARC)(prognostics.o)
$(ARC)(canyonalb.o): 
$(ARC)(albsnow_jls.o): $(ARC)(rad_param.o) $(ARC)(nstypes.o)
$(ARC)(matinv.o): $(ARC)(logging_mod.o)
$(ARC)(ftsa_jls.o): $(ARC)(fluxes.o) $(ARC)(c_0_dg_c.o) $(ARC)(ancil_info.o) $(ARC)(c_kappai.o)
$(ARC)(albpft_jls.o): $(ARC)(nstypes.o) $(ARC)(jules_mod.o) $(ARC)(pftparm.o) $(ARC)(switches.o)


# End of file.
