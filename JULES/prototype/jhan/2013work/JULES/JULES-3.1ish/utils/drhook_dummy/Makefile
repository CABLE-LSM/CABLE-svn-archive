include ../../Makefile.common.mk
THISFILE = ./Makefile

##############################
## Dependencies             ##
##############################
.PHONY : all clean
all : parkind1 yomhook $(THISFILE)

parkind1 :
	$(FC) -c $(FFLAGS) $(FPP) $(FPP_DEFS) parkind1.F90
        
yomhook :
	$(FC) -c $(FFLAGS) $(FPP) $(FPP_DEFS) yomhook.F90

clean : 
	@$(RM) $(wildcard ./*.o) $(wildcard ./*$(MOD_FSUF))
