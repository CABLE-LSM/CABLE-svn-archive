
#FC = gfortran 
FC = ifort
FLINK = -c 
FFLAGS = -O0 -fpp 
#FFLAGS = -x f95-cpp-input 
AUX= subr/
#shine
#FPGPLOTLIBS = -L/home/srb001/Public/pgplot/ifort -lpgplot -L/usr/lib64 -lpng -L/usr/lib64 -lz -L/usr/X11R6/lib64/ -lX11  
#cherax
FPGPLOTLIBS = -L/tools/pgplot/5.2.2 -lpgplot -L/usr/lib/ -lgfortran -L/usr/lib -lpng -L/usr/lib -lz -L/usr/X11R6/lib -lX11  
#xe.nci.org.au
#FPGPLOTLIBS = -L/home/599/jxs599/Public/pgplot/ -lpgplot -L/usr/lib64 -lpng -L/usr/lib64 -lz -L/usr/X11R6/lib64/ -lX11  
OBJ= debug_common.o debug_read.o debug_write.o debug_comp.o debug_plot.o

diag_main: diag_main.f90
	$(FC) $(FFLAGS) -c $(AUX)debug_common.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_read.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_write.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_comp.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_plot.f90	
	$(FC) $(FFLAGS) -o diag_main diag_main.f90 $(OBJ) $(FPGPLOTLIBS)




mock_cable: mock_cable.f90
	$(FC) $(FFLAGS) -c cable_debug.f90	
	$(FC) $(FFLAGS) -o mock_cable mock_cable.f90	cable_debug.o




clean:
	rm -f *.o *~ *.mod






