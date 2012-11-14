
FC = gfortran 
#FC = /apps/intel-fc/11.0.074/bin/intel64/ifort
#FC = ifort
#FFLAGS = -O0 -g  -i8 -r8 

diff_main: diff_main.f90
	$(FC) $(FFLAGS) -o diff_main diff_main.f90 

clean:
	rm -f *.o *~ *.mod diff_main






