include ../../Makefile.common.mk
THISFILE = ./Makefile

##############################
## Dependencies             ##
##############################
.PHONY : all clean
all : parkind1 yomhook $(THISFILE)

parkind1 : parkind1.F90 $(THISFILE)
	$(FC) -c $(FFLAGS) $(FPP) $(FPP_DEFS) parkind1.F90
        
yomhook : yomhook.F90 $(THISFILE)
	$(FC) -c $(FFLAGS) $(FPP) $(FPP_DEFS) yomhook.F90

clean : 
	@$(RM) $(wildcard ./*.o) $(wildcard ./*$(MOD_FSUF))
