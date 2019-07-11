include ../../Makefile.common.mk
THISFILE = ./Makefile

OBJ_CDF_DUMMY=jules_netcdf_dummy.o

##############################
## Dependencies             ##
##############################
.PHONY : all clean
all : $(OBJ_CDF_DUMMY) $(THISFILE)

$(OBJ_CDF_DUMMY) : jules_netcdf_dummy.F90 $(THISFILE)
	$(FC) -c $(FFLAGS) $(FPP) $(FPP_DEFS) jules_netcdf_dummy.F90
clean : 
	@$(RM) $(OBJ_CDF_DUMMY) $(wildcard ./*$(MOD_FSUF))
