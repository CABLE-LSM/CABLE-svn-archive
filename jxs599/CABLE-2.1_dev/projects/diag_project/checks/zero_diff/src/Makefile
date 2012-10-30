
#FC = gfortran 
#FC = /apps/intel-fc/11.0.074/bin/intel64/ifort
FC = ifort
FFLAGS = -O0 -g  -i8 -r8 
AUX= subr/
OBJ= debug_common.o debug_read.o 

diff_main: diff_main.f90
	$(FC) $(FFLAGS) -c $(AUX)debug_common.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_read.f90	
	$(FC) $(FFLAGS) -o diff_main diff_main.f90 $(OBJ) 



clean:
	rm -f *.o *~ *.mod diff_main






