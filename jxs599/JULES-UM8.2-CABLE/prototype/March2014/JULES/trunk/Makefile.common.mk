.SUFFIXES: .a .F90

.F90.a:
	$(FC) -c $(FFLAGS) $(FPP) $(FPP_DEFS) $< $(MOD_PATH) $(FPP_INC_PATH)
	@$(AR) $@ $(%F:.F90=.o)
	@$(RM) $(%F:.F90=.o)

######################################
## Architecture specific variables. ##
######################################
AR=ar r
EXTRACT=ar x
RM=rm -f
MOD_FSUF=.mod
