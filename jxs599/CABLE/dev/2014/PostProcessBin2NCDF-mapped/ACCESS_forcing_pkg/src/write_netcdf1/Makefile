
AUX= subr/

FLINK = -c 

#@#jiggle
#@FC = gfortran 
#@FFLAGS = -O0 -x f95-cpp-input -ffree-form -ffree-line-length-0    
#@NCDIR=/usr/local/lib
#@NCMOD=/usr/local/include

FC = ifort
FFLAGS = -O0 -fpp -traceback -g -I $(NCMOD) 
#@NCDF_ROOT=/apps/netcdf/4.2.1
NCDIR=$(NCDF_ROOT)/lib/Intel
NCMOD=$(NCDF_ROOT)/include/Intel

#cherax
#NCDIR = -L /apps/netcdf/4.1.3/lib -lnetcdf -lnetcdff
#NCMOD = /apps/netcdf/4.1.3/include

OBJ= debug_read.o debug_write.o

ncdf_main: ncdf_main.f90
	$(FC) $(FFLAGS) -c $(AUX)debug_read.f90
	$(FC) $(FFLAGS) -c $(AUX)debug_write.f90 -I$(NCMOD)
	$(FC) $(FFLAGS) -c ncdf_main.f90 
	$(FC) -o ncdf_main ncdf_main.o $(OBJ) -L $(NCDIR) -lnetcdf -lnetcdff 

clean:
	rm -f *.o *~ *.mod 






