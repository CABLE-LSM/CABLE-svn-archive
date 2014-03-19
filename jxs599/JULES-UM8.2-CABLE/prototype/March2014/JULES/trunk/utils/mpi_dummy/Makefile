include ../../Makefile.common.mk
THISFILE = ./Makefile

##############################
## Dependencies             ##
##############################
.PHONY : all clean
all : mpi.o $(THISFILE)

mpi.o : mpi.F90 $(THISFILE)
	$(FC) -c $(FFLAGS) $(FPP) $(FPP_DEFS) mpi.F90
        
clean : 
	@$(RM) $(wildcard ./*.o) $(wildcard ./*$(MOD_FSUF))
