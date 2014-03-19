# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(tridag.o) $(ARC)(snowtherm.o) $(ARC)(snow.o) $(ARC)(layersnow.o) $(ARC)(compactsnow.o) $(ARC)(relayersnow.o) $(ARC)(canopysnow.o) $(ARC)(snowpack.o) $(ARC)(snowgrain.o) $(ARC)(snow_intctl.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(tridag.o): 
$(ARC)(snowtherm.o): $(ARC)(ancil_info.o) $(ARC)(c_rhowat.o) $(ARC)(c_perma.o)
$(ARC)(snow.o): $(ARC)(switches.o) $(ARC)(c_lheat.o) $(ARC)(ancil_info.o) $(ARC)(snow_param.o)
$(ARC)(layersnow.o): $(ARC)(ancil_info.o) $(ARC)(snow_param.o)
$(ARC)(compactsnow.o): $(ARC)(c_g.o) $(ARC)(ancil_info.o)
$(ARC)(relayersnow.o): $(ARC)(switches.o) $(ARC)(rad_param.o) $(ARC)(ancil_info.o) $(ARC)(snow_param.o) $(ARC)(c_0_dg_c.o) $(ARC)(c_perma.o)
$(ARC)(canopysnow.o): $(ARC)(snow_param.o)
$(ARC)(snowpack.o): $(ARC)(c_0_dg_c.o) $(ARC)(snow_param.o) $(ARC)(c_lheat.o) $(ARC)(ancil_info.o) $(ARC)(soil_param.o) $(ARC)(c_densty.o) $(ARC)(c_perma.o)
$(ARC)(snowgrain.o): $(ARC)(c_0_dg_c.o) $(ARC)(rad_param.o) $(ARC)(c_pi.o) $(ARC)(ancil_info.o) $(ARC)(c_rmol.o)
$(ARC)(snow_intctl.o): 


# End of file.
