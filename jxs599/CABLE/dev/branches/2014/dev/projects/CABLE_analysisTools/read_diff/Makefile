
#FC = gfortran 
#FC = /apps/intel-fc/11.0.074/bin/intel64/ifort
FC = ifort
FLINK = -c 
FFLAGS = -O2
AUX= subr/
OBJ= debug_common.o debug_read.o 

diff_main: diff_main.f90
	$(FC) $(FFLAGS) -c $(AUX)debug_common.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_read.f90	
	$(FC) $(FFLAGS) -o diff_main diff_main.f90 $(OBJ) 



clean:
	rm -f *.o *~ *.mod diff_main






