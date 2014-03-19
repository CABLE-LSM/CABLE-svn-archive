# This is a generated file - do not modify
        
include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(flake_init.o) $(ARC)(flake_init_0a.o) $(ARC)(jules_init.o) $(ARC)(allocate_jules_arrays.o) $(ARC)(init_urban.o) $(ARC)(total_snow_um.o) $(ARC)(read_jules_namelists.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(flake_init.o): $(THISFILE)
$(ARC)(flake_init_0a.o): $(THISFILE)
$(ARC)(jules_init.o): $(THISFILE)
$(ARC)(allocate_jules_arrays.o): $(THISFILE)
$(ARC)(init_urban.o): $(ARC)(allocate_arrays.o) $(ARC)(misc_utils.o) $(ARC)(readwrite_mod.o) $(ARC)(file_utils.o) $(ARC)(inout.o) $(ARC)(jules_netcdf.o) $(ARC)(switches.o) $(ARC)(nvegparm.o) $(ARC)(c_z0h_z0m.o) $(ARC)(urban_param.o) $(ARC)(switches_urban.o) $(ARC)(ancil_info.o) $(ARC)(nstypes.o) $(THISFILE)
$(ARC)(total_snow_um.o): $(THISFILE)
$(ARC)(read_jules_namelists.o): $(THISFILE)


.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

# End of file.
