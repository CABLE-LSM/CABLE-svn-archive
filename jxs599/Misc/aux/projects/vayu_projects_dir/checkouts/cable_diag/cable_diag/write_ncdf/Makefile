
#FC = gfortran 
#FC = /apps/intel-fc/11.0.074/bin/intel64/ifort
FC = ifort
FLINK = -c 
FFLAGS = -O2 -fpp -traceback -g -I $(NCMOD) 
AUX= subr/

#vayu.nci.org.au
NCDIR = -L /apps/netcdf/3.6.3/lib/Intel/ -lnetcdf
NCMOD = /apps/netcdf/3.6.3/include/Intel/netcdf.mod

FPGPLOTLIBS = -L/home/599/jxs599/Public/pgplot/ -lpgplot -L/usr/lib64 -lpng -L/usr/lib64 -lz -L/usr/X11R6/lib64/ -lX11 
#OBJ= debug_common.o debug_read.o debug_write.o debug_comp.o debug_plot.o
OBJ= debug_read.o debug_write.o

ncdf_main: ncdf_main.f90
	$(FC) $(FFLAGS) -c $(AUX)debug_read.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_write.f90	
#	$(FC) $(FFLAGS) -c $(AUX)debug_comp.f90	
#	$(FC) $(FFLAGS) -c $(AUX)debug_plot.f90	
	$(FC) $(FFLAGS) -o ncdf_main ncdf_main.f90 $(OBJ) $(NCDIR)  



clean:
	rm -f *.o *~ *.mod 






