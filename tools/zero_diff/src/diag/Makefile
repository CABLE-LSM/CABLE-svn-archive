
#FC = gfortran 
#FC = /apps/intel-fc/11.0.074/bin/intel64/ifort
FC = ifort
FFLAGS = -O0 -g -traceback  #-i8 -r8 
#FFLAGS = -O2 -x f95-cpp-input 

diff_main_r4: diff_main.f90
	$(FC) $(FFLAGS) -o diff_main_r4 diff_main.f90 

diff_main_r8: diff_main.f90
	#$(FC) $(FFLAGS) -i8 -r8 -o diff_main_r8 diff_main.f90 
	$(FC) $(FFLAGS) -fdefault-real-8 -fdefault-integer-8 -o diff_main_r8 diff_main.f90 

clean:
	rm -f *.o *~ *.mod diff_main






